use crate::binread::{BinReader, Endian};
use crate::error::{NmrIoError, Result};
use crate::jcamp;
use crate::virtual_fs::{LoadedDataset, VirtualTree};
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D};
use num_complex::Complex64;
use std::collections::BTreeMap;

/// Parse Magritek Spinsolve `acqu.par` / `proc.par` key-value text file.
pub fn parse_spinsolve_par(bytes: &[u8]) -> BTreeMap<String, String> {
    let text = String::from_utf8_lossy(bytes);
    let mut map = BTreeMap::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some((k, v)) = trimmed.split_once('=') {
            let key = k.trim().to_string();
            let val = v.trim().trim_matches('"').to_string();
            map.insert(key, val);
        }
    }
    map
}

/// Load a Magritek Spinsolve experiment directory (`acqu.par` + `nmr_fid.dx` / `data.1d` / `spectrum.1d`).
pub fn load_spinsolve_tree(tree: &VirtualTree) -> Result<LoadedDataset> {
    let acqu = tree
        .get("acqu.par")
        .map(parse_spinsolve_par)
        .unwrap_or_default();

    if let Some(dx_bytes) = tree.get("nmr_fid.dx") {
        let text = String::from_utf8_lossy(dx_bytes);
        let mut sig = jcamp::load_jcamp_str("nmr_fid.dx", &text)?;
        apply_spinsolve_acqu(&mut sig.meta, &acqu);
        return Ok(LoadedDataset::D1(sig));
    }

    let (bin_name, bin_bytes, domain) = if let Some(b) = tree.get("data.1d") {
        ("data.1d", b, Domain::Time)
    } else if let Some(b) = tree.get("fid.1d") {
        ("fid.1d", b, Domain::Time)
    } else if let Some(b) = tree.get("spectrum.1d") {
        ("spectrum.1d", b, Domain::Frequency)
    } else if let Some(b) = tree.get("spectrum_processed.1d") {
        ("spectrum_processed.1d", b, Domain::Frequency)
    } else {
        return Err(NmrIoError::MissingFile(
            "nmr_fid.dx, data.1d, fid.1d, or spectrum.1d".into(),
        ));
    };

    if bin_bytes.len() < 32 {
        return Err(NmrIoError::Malformed {
            format: "Spinsolve .1d",
            detail: format!("{bin_name} shorter than 32-byte header"),
        });
    }

    // 32-byte little-endian header: 8 * u32 (owner, format, version, dataType, xDim, yDim, zDim, qDim)
    let mut r = BinReader::new(bin_bytes, Endian::Little, "Spinsolve .1d");
    r.skip(32)?;
    let n_f32 = (bin_bytes.len() - 32) / 4;
    let mut floats = Vec::with_capacity(n_f32);
    for _ in 0..n_f32 {
        floats.push(r.read_f32()? as f64);
    }

    // First 1/3 is x-axis (time in s or chemical shift in ppm), remaining 2/3 is interleaved [re, im]
    let split = floats.len() / 3;
    let payload = &floats[split..];
    let n_pts = payload.len() / 2;
    let data: Vec<Complex64> = (0..n_pts)
        .map(|i| Complex64::new(payload[2 * i], payload[2 * i + 1]))
        .collect();

    let mut meta = AcquisitionParams {
        title: tree.root_name.clone(),
        format: NmrFormat::Spinsolve,
        ..Default::default()
    };
    apply_spinsolve_acqu(&mut meta, &acqu);

    Ok(LoadedDataset::D1(Signal1D::new(data, domain, meta)))
}

fn apply_spinsolve_acqu(meta: &mut AcquisitionParams, acqu: &BTreeMap<String, String>) {
    if acqu.is_empty() {
        return;
    }
    meta.format = NmrFormat::Spinsolve;
    if let Some(bw_khz) = acqu.get("bandwidth").and_then(|s| s.parse::<f64>().ok()) {
        meta.axis.spectral_width_hz = bw_khz * 1000.0;
    }
    if let Some(sf) = acqu.get("b1Freq").and_then(|s| s.parse::<f64>().ok()) {
        meta.axis.spectrometer_freq_mhz = sf;
    }
    if let Some(low_hz) = acqu.get("lowestFrequency").and_then(|s| s.parse::<f64>().ok()) {
        let car_hz = low_hz + 0.5 * meta.axis.spectral_width_hz;
        if meta.axis.spectrometer_freq_mhz > 0.0 {
            meta.axis.carrier_ppm = car_hz / meta.axis.spectrometer_freq_mhz;
        }
    }
    if let Some(nuc) = acqu.get("rxChannel") {
        meta.axis.nucleus = nuc.clone();
    }
    if let Some(ns) = acqu.get("nrScans").and_then(|s| s.parse::<usize>().ok()) {
        meta.num_scans = ns;
    }
    meta.axis = AxisMeta {
        nucleus: meta.axis.nucleus.clone(),
        spectrometer_freq_mhz: meta.axis.spectrometer_freq_mhz,
        spectral_width_hz: meta.axis.spectral_width_hz,
        offset_ppm: meta.axis.offset_ppm,
        carrier_ppm: meta.axis.carrier_ppm,
    };
    for (k, v) in acqu {
        meta.vendor_params.insert(k.clone(), v.clone());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_spinsolve_1d_binary() {
        let mut tree = VirtualTree::new("spinsolve_sample");
        let acqu = "bandwidth = 5.0\nb1Freq = 43.4\nlowestFrequency = -500.0\nrxChannel = \"1H\"\nnrScans = 4\n";
        let mut bin = vec![0u8; 32];
        // 4 points -> 4 xaxis floats + 4*2 complex floats = 12 floats
        for x in 0..4 {
            bin.extend_from_slice(&(x as f32 * 0.01).to_le_bytes());
        }
        for i in 0..4 {
            bin.extend_from_slice(&(i as f32 + 1.0).to_le_bytes());
            bin.extend_from_slice(&(-(i as f32)).to_le_bytes());
        }
        tree.insert("acqu.par", acqu.as_bytes().to_vec());
        tree.insert("data.1d", bin);

        let sig = load_spinsolve_tree(&tree).unwrap().into_1d().unwrap();
        assert_eq!(sig.len(), 4);
        assert_eq!(sig.meta.num_scans, 4);
        assert!((sig.meta.axis.spectral_width_hz - 5000.0).abs() < 1e-6);
        assert!((sig.data[2].re - 3.0).abs() < 1e-6);
        assert!((sig.data[2].im + 2.0).abs() < 1e-6);
    }
}
