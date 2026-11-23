use crate::binread::{BinReader, Endian};
use crate::error::{NmrIoError, Result};
use crate::virtual_fs::{LoadedDataset, VirtualTree};
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D, Signal2D};
use num_complex::Complex64;
use std::collections::BTreeMap;

/// Parse Bruker JCAMP-DX style parameter files (`acqus`, `acqu2s`, `procs`, `proc2s`).
pub fn parse_bruker_params(bytes: &[u8]) -> BTreeMap<String, String> {
    let text = String::from_utf8_lossy(bytes);
    let mut map = BTreeMap::new();
    let mut current_key: Option<String> = None;
    let mut current_val = String::new();
    let mut in_array = false;

    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.starts_with("$$") || line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix("##") {
            if let Some(prev) = current_key.take() {
                map.insert(prev, clean_bruker_value(&current_val));
                current_val.clear();
            }
            in_array = false;

            if rest.eq_ignore_ascii_case("END=") {
                break;
            }
            if let Some((k_part, v_part)) = rest.split_once('=') {
                let key = k_part.trim().trim_start_matches('$').to_string();
                let val = v_part.trim();
                if val.starts_with('(') && val.contains("..") {
                    // Array header e.g. "(0..31)" -> values follow on subsequent lines
                    current_key = Some(key);
                    in_array = true;
                } else if val.starts_with('<') && !val.ends_with('>') {
                    current_key = Some(key);
                    current_val.push_str(val);
                } else {
                    map.insert(key, clean_bruker_value(val));
                }
            }
        } else if let Some(ref _k) = current_key {
            if !current_val.is_empty() {
                current_val.push(' ');
            }
            current_val.push_str(line);
            if !in_array && line.ends_with('>') {
                let prev = current_key.take().unwrap();
                map.insert(prev, clean_bruker_value(&current_val));
                current_val.clear();
            }
        }
    }

    if let Some(prev) = current_key {
        map.insert(prev, clean_bruker_value(&current_val));
    }
    map
}

fn clean_bruker_value(v: &str) -> String {
    let s = v.trim();
    if s.starts_with('<') && s.ends_with('>') && s.len() >= 2 {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

/// Lookup standard Bruker digital filter group delay from `DSPFVS` and `DECIM`
/// when `GRPDLY` is absent or negative.
pub fn lookup_bruker_group_delay(dspfvs: i32, decim: i32) -> Option<f64> {
    // Standard Bruker Avance DSP firmware group-delay tables (points).
    let val = match (dspfvs, decim) {
        (10, 2) => 44.7500,
        (10, 3) => 33.5000,
        (10, 4) => 66.6250,
        (10, 6) => 59.0833,
        (10, 8) => 68.5625,
        (10, 12) => 60.3750,
        (10, 16) => 69.5313,
        (10, 24) => 61.0208,
        (10, 32) => 70.0156,
        (10, 48) => 61.3438,
        (10, 64) => 70.2578,
        (10, 96) => 61.5052,
        (10, 128) => 70.3789,
        (10, 192) => 61.5859,
        (10, 256) => 70.4395,
        (10, 384) => 61.6263,
        (10, 512) => 70.4697,
        (10, 768) => 61.6465,
        (10, 1024) => 70.4849,
        (10, 1536) => 61.6566,
        (10, 2048) => 70.4924,
        (11, 2) => 46.0000,
        (11, 3) => 36.5000,
        (11, 4) => 48.0000,
        (11, 6) => 50.1667,
        (11, 8) => 53.2500,
        (11, 12) => 69.5000,
        (11, 16) => 72.2500,
        (11, 24) => 70.1667,
        (11, 32) => 72.7500,
        (11, 48) => 70.5000,
        (11, 64) => 73.0000,
        (11, 96) => 70.6667,
        (11, 128) => 72.5000,
        (11, 192) => 71.3333,
        (11, 256) => 72.2500,
        (11, 384) => 71.6667,
        (11, 512) => 72.1250,
        (11, 768) => 71.8333,
        (11, 1024) => 72.0625,
        (11, 1536) => 71.9167,
        (11, 2048) => 72.0313,
        (12, 2) => 46.311,
        (12, 3) => 36.530,
        (12, 4) => 47.870,
        (12, 6) => 50.229,
        (12, 8) => 53.289,
        (12, 12) => 69.551,
        (12, 16) => 71.600,
        (12, 24) => 70.184,
        (12, 32) => 72.138,
        (12, 48) => 70.528,
        (12, 64) => 72.348,
        (12, 96) => 70.700,
        (12, 128) => 72.524,
        (12, 192) => 71.333,
        (12, 256) => 72.250,
        (12, 384) => 71.667,
        (12, 512) => 72.125,
        (12, 768) => 71.833,
        (12, 1024) => 72.063,
        (12, 1536) => 71.917,
        (12, 2048) => 72.031,
        (13, 2) => 2.750,
        (13, 3) => 2.833,
        (13, 4) => 2.875,
        (13, 6) => 2.917,
        (13, 8) => 2.938,
        (13, 12) => 2.958,
        (13, 16) => 2.969,
        (13, 24) => 2.979,
        (13, 32) => 2.984,
        (13, 48) => 2.989,
        (13, 64) => 2.992,
        (13, 96) => 2.995,
        _ => return None,
    };
    Some(val)
}

/// Load a Bruker TopSpin dataset (1D `fid`, 2D `ser`, or processed `pdata/1/1r`+`1i`).
pub fn load_bruker_tree(tree: &VirtualTree) -> Result<LoadedDataset> {
    let acqus = tree
        .get("acqus")
        .map(parse_bruker_params)
        .unwrap_or_default();
    let acqu2s = tree.get("acqu2s").map(parse_bruker_params);
    let procs = tree.get("procs").or_else(|| tree.get("pdata/1/procs")).map(parse_bruker_params);
    let proc2s = tree.get("proc2s").or_else(|| tree.get("pdata/1/proc2s")).map(parse_bruker_params);

    let title = tree
        .get("pdata/1/title")
        .map(|b| String::from_utf8_lossy(b).trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| tree.root_name.clone());

    let axis = build_bruker_axis(&acqus, procs.as_ref());
    let indirect_axis = acqu2s
        .as_ref()
        .map(|a2| build_bruker_axis(a2, proc2s.as_ref()));

    let group_delay = acqus
        .get("GRPDLY")
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|&v| v > 0.0)
        .or_else(|| {
            let dspfvs = acqus.get("DSPFVS")?.parse::<i32>().ok()?;
            let decim = acqus.get("DECIM")?.parse::<i32>().ok()?;
            lookup_bruker_group_delay(dspfvs, decim)
        });

    let num_scans = acqus
        .get("NS")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1);

    let mut vendor_params = acqus.clone();
    if let Some(p) = &procs {
        for (k, v) in p {
            vendor_params.insert(format!("procs.{k}"), v.clone());
        }
    }

    let meta = AcquisitionParams {
        title,
        format: NmrFormat::Bruker,
        axis,
        indirect_axis,
        group_delay,
        num_scans,
        vendor_params,
    };

    // Prefer 2D `ser` when acqu2s exists with > 1 row
    if let Some(ser_bytes) = tree.get("ser") {
        let td_direct = acqus
            .get("TD")
            .and_then(|s| s.parse::<usize>().ok())
            .ok_or_else(|| NmrIoError::MissingParameter("TD in acqus".into()))?;
        let td_indirect = acqu2s
            .as_ref()
            .and_then(|m| m.get("TD"))
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(1);
        let endian = bruker_endian(&acqus, "BYTORDP");
        let dtypa = bruker_dtypa(&acqus, "DTYPA");
        let reals = read_bruker_numbers(ser_bytes, endian, dtypa, 0)?;
        let cols = (td_direct / 2).max(1);
        let rows = td_indirect.min(reals.len() / (2 * cols)).max(1);
        let mut data = Vec::with_capacity(rows * cols);
        for idx in 0..(rows * cols) {
            data.push(Complex64::new(reals[2 * idx], reals[2 * idx + 1]));
        }
        return Ok(LoadedDataset::D2(Signal2D::new(
            rows,
            cols,
            data,
            Domain::Time,
            Domain::Time,
            meta,
        )));
    }

    // 1D raw FID
    if let Some(fid_bytes) = tree.get("fid") {
        let endian = bruker_endian(&acqus, "BYTORDP");
        let dtypa = bruker_dtypa(&acqus, "DTYPA");
        let reals = read_bruker_numbers(fid_bytes, endian, dtypa, 0)?;
        let td = acqus
            .get("TD")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(reals.len());
        let complex_len = (td / 2).min(reals.len() / 2);
        let data: Vec<Complex64> = (0..complex_len)
            .map(|i| Complex64::new(reals[2 * i], reals[2 * i + 1]))
            .collect();
        return Ok(LoadedDataset::D1(Signal1D::new(data, Domain::Time, meta)));
    }

    // 1D processed spectrum (`pdata/1/1r` + optional `1i`)
    if let Some(r_bytes) = tree.get("pdata/1/1r") {
        let p = procs.as_ref().ok_or_else(|| NmrIoError::MissingFile("pdata/1/procs".into()))?;
        let endian = bruker_endian(p, "BYTORDP");
        let dtypa = bruker_dtypa(p, "DTYPA");
        let nc_proc = p
            .get("NC_proc")
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0);
        let real_vals = read_bruker_numbers(r_bytes, endian, dtypa, nc_proc)?;
        let imag_vals = if let Some(i_bytes) = tree.get("pdata/1/1i") {
            read_bruker_numbers(i_bytes, endian, dtypa, nc_proc)?
        } else {
            vec![0.0; real_vals.len()]
        };
        let data: Vec<Complex64> = real_vals
            .into_iter()
            .zip(imag_vals.into_iter().chain(std::iter::repeat(0.0)))
            .map(|(r, im)| Complex64::new(r, im))
            .collect();
        let mut sig = Signal1D::new(data, Domain::Frequency, meta);
        sig.digital_filter_removed = true;
        return Ok(LoadedDataset::D1(sig));
    }

    Err(NmrIoError::MissingFile("fid, ser, or pdata/1/1r".into()))
}

fn build_bruker_axis(acqu: &BTreeMap<String, String>, procs: Option<&BTreeMap<String, String>>) -> AxisMeta {
    let nucleus = acqu
        .get("NUC1")
        .cloned()
        .unwrap_or_else(|| "1H".to_string());
    let sf = procs
        .and_then(|p| p.get("SF"))
        .or_else(|| acqu.get("BF1"))
        .or_else(|| acqu.get("SFO1"))
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|&v| v > 0.0)
        .unwrap_or(500.0);
    let sw_hz = procs
        .and_then(|p| p.get("SW_p"))
        .and_then(|s| s.parse::<f64>().ok())
        .or_else(|| acqu.get("SW_h").and_then(|s| s.parse::<f64>().ok()))
        .or_else(|| acqu.get("SW").and_then(|s| s.parse::<f64>().ok()).map(|ppm| ppm * sf))
        .unwrap_or(6000.0);
    let offset_ppm = procs
        .and_then(|p| p.get("OFFSET"))
        .and_then(|s| s.parse::<f64>().ok());
    let carrier_ppm = acqu
        .get("O1")
        .and_then(|s| s.parse::<f64>().ok())
        .map(|o1_hz| o1_hz / sf)
        .unwrap_or(5.0);
    AxisMeta {
        nucleus,
        spectrometer_freq_mhz: sf,
        spectral_width_hz: sw_hz,
        offset_ppm,
        carrier_ppm,
    }
}

fn bruker_endian(map: &BTreeMap<String, String>, key: &str) -> Endian {
    match map.get(key).and_then(|s| s.parse::<i32>().ok()).unwrap_or(0) {
        1 => Endian::Big,
        _ => Endian::Little,
    }
}

fn bruker_dtypa(map: &BTreeMap<String, String>, key: &str) -> i32 {
    map.get(key).and_then(|s| s.parse::<i32>().ok()).unwrap_or(0)
}

fn read_bruker_numbers(bytes: &[u8], endian: Endian, dtypa: i32, nc_proc: i32) -> Result<Vec<f64>> {
    let scale = 2.0_f64.powi(nc_proc);
    let mut reader = BinReader::new(bytes, endian, "Bruker");
    let mut out = Vec::new();
    match dtypa {
        2 => {
            let n = bytes.len() / 8;
            out.reserve(n);
            for _ in 0..n {
                out.push(reader.read_f64()? * scale);
            }
        }
        1 => {
            let n = bytes.len() / 4;
            out.reserve(n);
            for _ in 0..n {
                out.push(reader.read_f32()? as f64 * scale);
            }
        }
        _ => {
            let n = bytes.len() / 4;
            out.reserve(n);
            for _ in 0..n {
                out.push(reader.read_i32()? as f64 * scale);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_synthetic_bruker_1d_fid_and_group_delay() {
        let mut tree = VirtualTree::new("ethanol_exp");
        let acqus = "\
##TITLE= Ethanol 1H
##$TD= 16
##$BYTORDP= 0
##$DTYPA= 0
##$BF1= 500.0
##$SFO1= 500.0025
##$O1= 2500.0
##$SW_h= 6000.0
##$GRPDLY= 67.98
##$NS= 16
##$NUC1= <1H>
##END=
";
        let mut fid = Vec::new();
        for i in 0..8 {
            fid.extend_from_slice(&((i as i32 + 1) * 100).to_le_bytes());
            fid.extend_from_slice(&((-((i as i32) + 1)) * 50).to_le_bytes());
        }
        tree.insert("acqus", acqus.as_bytes().to_vec());
        tree.insert("fid", fid);

        let ds = load_bruker_tree(&tree).unwrap().into_1d().unwrap();
        assert_eq!(ds.len(), 8);
        assert_eq!(ds.domain, Domain::Time);
        assert_eq!(ds.meta.num_scans, 16);
        assert!((ds.meta.group_delay.unwrap() - 67.98).abs() < 1e-6);
        assert!((ds.data[0].re - 100.0).abs() < 1e-9);
        assert!((ds.data[0].im + 50.0).abs() < 1e-9);
    }
}
