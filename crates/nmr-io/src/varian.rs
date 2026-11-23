use crate::binread::{BinReader, Endian};
use crate::error::{NmrIoError, Result};
use crate::virtual_fs::{LoadedDataset, VirtualTree};
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D, Signal2D};
use num_complex::Complex64;
use std::collections::BTreeMap;

/// Status flags in Varian `fid` header.
const S_32: i16 = 0x04;
const S_FLOAT: i16 = 0x08;
const S_SPEC: i16 = 0x02;

/// Parse Varian/Agilent `procpar` text file into a key-value map.
pub fn parse_procpar(bytes: &[u8]) -> BTreeMap<String, String> {
    let text = String::from_utf8_lossy(bytes);
    let lines: Vec<&str> = text.lines().collect();
    let mut map = BTreeMap::new();
    let mut i = 0usize;

    while i < lines.len() {
        let header_parts: Vec<&str> = lines[i].split_whitespace().collect();
        // Parameter header lines start with parameter name followed by >= 10 numeric fields
        if header_parts.len() >= 11 && header_parts[1].parse::<i32>().is_ok() {
            let name = header_parts[0].to_string();
            let basictype: i32 = header_parts[2].parse().unwrap_or(1);
            if i + 1 < lines.len() {
                let val_line = lines[i + 1].trim();
                let mut parsed_val = String::new();
                if basictype == 2 {
                    // String parameter: first token is count, rest may contain quoted strings across lines
                    if let Some((_, rest)) = val_line.split_once(char::is_whitespace) {
                        parsed_val = rest.trim().trim_matches('"').to_string();
                    }
                } else {
                    let toks: Vec<&str> = val_line.split_whitespace().collect();
                    if toks.len() >= 2 {
                        parsed_val = toks[1..].join(" ");
                    }
                }
                map.insert(name, parsed_val);
                // Skip header line, value line(s), and enum line
                i += 3;
                continue;
            }
        }
        i += 1;
    }
    map
}

/// Load a Varian/Agilent `.fid` directory (`fid` + `procpar`).
pub fn load_varian_tree(tree: &VirtualTree) -> Result<LoadedDataset> {
    let fid_bytes = tree
        .get("fid")
        .ok_or_else(|| NmrIoError::MissingFile("fid".into()))?;
    let procpar = tree
        .get("procpar")
        .map(parse_procpar)
        .unwrap_or_default();

    let mut r = BinReader::new(fid_bytes, Endian::Big, "Varian");
    let nblocks = r.read_i32()? as usize;
    let ntraces = r.read_i32()? as usize;
    let np = r.read_i32()? as usize;
    let _ebytes = r.read_i32()?;
    let _tbytes = r.read_i32()?;
    let _bbytes = r.read_i32()?;
    let _vers_id = r.read_i16()?;
    let status = r.read_i16()?;
    let nbheaders = (r.read_i32()? as usize).max(1);

    let complex_pts = (np / 2).max(1);
    let total_traces = nblocks.saturating_mul(ntraces);
    let mut all_data = Vec::with_capacity(total_traces * complex_pts);

    for _ in 0..nblocks {
        // Skip block headers (28 bytes each)
        r.skip(nbheaders * 28)?;
        for _ in 0..ntraces {
            for _ in 0..complex_pts {
                let (re, im) = if (status & S_FLOAT) != 0 {
                    (r.read_f32()? as f64, r.read_f32()? as f64)
                } else if (status & S_32) != 0 {
                    (r.read_i32()? as f64, r.read_i32()? as f64)
                } else {
                    (r.read_i16()? as f64, r.read_i16()? as f64)
                };
                all_data.push(Complex64::new(re, im));
            }
        }
    }

    let domain = if (status & S_SPEC) != 0 {
        Domain::Frequency
    } else {
        Domain::Time
    };

    let sf = procpar
        .get("sfrq")
        .and_then(|s| s.split_whitespace().next())
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|&v| v > 0.0)
        .unwrap_or(500.0);
    let sw = procpar
        .get("sw")
        .and_then(|s| s.split_whitespace().next())
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(6000.0);
    let rfp = procpar
        .get("rfp")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);
    let rfl = procpar
        .get("rfl")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);
    let car_ppm = ((sw / 2.0) - rfl + rfp) / sf;
    let nucleus = procpar.get("tn").cloned().unwrap_or_else(|| "1H".to_string());
    let num_scans = procpar
        .get("nt")
        .and_then(|s| s.split_whitespace().next())
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1);

    let meta = AcquisitionParams {
        title: tree.root_name.clone(),
        format: NmrFormat::Varian,
        axis: AxisMeta {
            nucleus,
            spectrometer_freq_mhz: sf,
            spectral_width_hz: sw,
            offset_ppm: None,
            carrier_ppm: if car_ppm.is_finite() { car_ppm } else { 5.0 },
        },
        indirect_axis: None,
        group_delay: None,
        num_scans,
        vendor_params: procpar,
    };

    if total_traces > 1 {
        Ok(LoadedDataset::D2(Signal2D::new(
            total_traces,
            complex_pts,
            all_data,
            domain,
            Domain::Time,
            meta,
        )))
    } else {
        Ok(LoadedDataset::D1(Signal1D::new(all_data, domain, meta)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_varian_big_endian_float_fid() {
        let mut tree = VirtualTree::new("varian_test.fid");
        let procpar = "\
sfrq 1 1 1e9 0 0 2 1 0 1 64
1 400.13
0
sw 1 1 1e9 0 0 2 1 0 1 64
1 4800.0
0
tn 2 2 0 0 0 2 1 0 1 64
1 \"H1\"
0
";
        let mut fid = Vec::new();
        // 32-byte big-endian file header: nblocks=1, ntraces=1, np=8 (4 complex), ebytes=4, tbytes=32, bbytes=60, vers=1, status=S_FLOAT|0x11, nbheaders=1
        for v in [1i32, 1, 8, 4, 32, 60] {
            fid.extend_from_slice(&v.to_be_bytes());
        }
        fid.extend_from_slice(&1i16.to_be_bytes());
        fid.extend_from_slice(&(S_FLOAT | 0x11).to_be_bytes());
        fid.extend_from_slice(&1i32.to_be_bytes());
        // 28-byte block header
        fid.extend_from_slice(&[0u8; 28]);
        // 4 complex floats
        for i in 0..4 {
            fid.extend_from_slice(&((i as f32) + 1.0).to_be_bytes());
            fid.extend_from_slice(&(-((i as f32) + 0.5)).to_be_bytes());
        }

        tree.insert("procpar", procpar.as_bytes().to_vec());
        tree.insert("fid", fid);

        let ds = load_varian_tree(&tree).unwrap().into_1d().unwrap();
        assert_eq!(ds.len(), 4);
        assert!((ds.meta.axis.spectrometer_freq_mhz - 400.13).abs() < 1e-6);
        assert!((ds.data[2].re - 3.0).abs() < 1e-6);
        assert!((ds.data[2].im + 2.5).abs() < 1e-6);
    }
}
