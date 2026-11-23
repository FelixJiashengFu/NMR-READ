use crate::error::{NmrIoError, Result};
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, PeakList, Signal1D};
use num_complex::Complex64;
use std::fmt::Write;

/// Export a 1D NMR signal (with optional baseline column) to CSV.
pub fn export_csv(sig: &Signal1D) -> String {
    let mut out = String::new();
    let x_name = match sig.domain {
        Domain::Frequency => "ppm",
        Domain::Time => "time_s",
    };
    let has_baseline = sig.baseline.as_ref().is_some_and(|b| b.len() == sig.len());
    if has_baseline {
        let _ = writeln!(out, "{x_name},real,imag,magnitude,baseline");
    } else {
        let _ = writeln!(out, "{x_name},real,imag,magnitude");
    }

    let x_vals = match sig.domain {
        Domain::Frequency => sig.x_axis_ppm(),
        Domain::Time => sig.x_axis_time_s(),
    };

    for (i, c) in sig.data.iter().enumerate() {
        let x = x_vals.get(i).copied().unwrap_or(i as f64);
        if let Some(bl) = sig.baseline.as_ref().filter(|b| b.len() == sig.len()) {
            let _ = writeln!(out, "{x:.8},{:.8},{:.8},{:.8},{:.8}", c.re, c.im, c.norm(), bl[i]);
        } else {
            let _ = writeln!(out, "{x:.8},{:.8},{:.8},{:.8}", c.re, c.im, c.norm());
        }
    }
    out
}

/// Export detected peaks to CSV.
pub fn export_peaks_csv(peaks: &PeakList) -> String {
    let mut out = String::from("index,ppm,hz,intensity,fwhm_hz,snr\n");
    for p in &peaks.peaks {
        let fwhm = p.fwhm_hz.map(|v| format!("{v:.4}")).unwrap_or_default();
        let _ = writeln!(
            out,
            "{},{:.6},{:.4},{:.6},{},{:.2}",
            p.index, p.ppm, p.hz, p.intensity, fwhm, p.snr
        );
    }
    out
}

/// Load a CSV/TSV spectrum or FID file (`x, real [, imag]`).
pub fn load_csv_str(filename: &str, text: &str) -> Result<Signal1D> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed
            .split(|c| c == ',' || c == '\t' || c == ';')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();
        if cols.len() < 2 {
            continue;
        }
        let Ok(x) = cols[0].parse::<f64>() else {
            continue; // header row
        };
        let Ok(re) = cols[1].parse::<f64>() else {
            continue;
        };
        let im = cols.get(2).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        xs.push(x);
        ys.push(Complex64::new(re, im));
    }
    if ys.is_empty() {
        return Err(NmrIoError::Malformed {
            format: "CSV/TSV",
            detail: "No numeric rows parsed".to_string(),
        });
    }

    let first = xs[0];
    let last = *xs.last().unwrap();
    let domain = if first > last {
        Domain::Frequency
    } else {
        Domain::Time
    };
    let sf = 500.0;
    let axis = if domain == Domain::Frequency {
        let sw_ppm = (first - last).abs().max(1e-6);
        AxisMeta {
            nucleus: "1H".into(),
            spectrometer_freq_mhz: sf,
            spectral_width_hz: sw_ppm * sf,
            offset_ppm: Some(first),
            carrier_ppm: first - 0.5 * sw_ppm,
        }
    } else {
        let dt = ((last - first).abs() / (xs.len().saturating_sub(1).max(1)) as f64).max(1e-6);
        AxisMeta {
            nucleus: "1H".into(),
            spectrometer_freq_mhz: sf,
            spectral_width_hz: 1.0 / dt,
            offset_ppm: None,
            carrier_ppm: 5.0,
        }
    };

    let meta = AcquisitionParams {
        title: filename.to_string(),
        format: NmrFormat::Csv,
        axis,
        ..Default::default()
    };
    Ok(Signal1D::new(ys, domain, meta))
}
