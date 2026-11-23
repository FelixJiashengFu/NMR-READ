use nmr_core::{Domain, Signal1D};
use num_complex::Complex64;
use rustfft::FftPlanner;
use std::f64::consts::TAU;

/// Remove Bruker/vendor digital filter group delay (`grpdly` points) via exact
/// Fourier linear phase shift followed by truncating the wrapped leading tail.
pub fn remove_digital_filter(sig: &Signal1D, group_delay_override: Option<f64>) -> Signal1D {
    if sig.domain != Domain::Time || sig.is_empty() {
        return sig.clone();
    }
    let delay = group_delay_override
        .or(sig.meta.group_delay)
        .unwrap_or(0.0);
    if delay <= 0.0 {
        return sig.clone();
    }

    let n = sig.len();
    let mut buf = sig.data.clone();
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);
    let ifft = planner.plan_fft_inverse(n);

    fft.process(&mut buf);
    // Shift left by `delay` samples: X[k] * exp(+i * 2*pi * k_signed * delay / N)
    let half = n / 2;
    for (k, val) in buf.iter_mut().enumerate() {
        let k_signed = if k <= half {
            k as f64
        } else {
            (k as isize - n as isize) as f64
        };
        let angle = TAU * k_signed * delay / (n as f64);
        let rot = Complex64::from_polar(1.0, angle);
        *val *= rot;
    }
    ifft.process(&mut buf);

    let inv_n = 1.0 / n as f64;
    let drop_tail = (delay.ceil() as usize).min(n.saturating_sub(1));
    let keep = n - drop_tail;
    let shifted: Vec<Complex64> = buf.into_iter().take(keep).map(|c| c * inv_n).collect();

    let mut out = sig.clone();
    out.data = shifted;
    out.digital_filter_removed = true;
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use nmr_core::AcquisitionParams;

    #[test]
    fn removes_integer_group_delay_cleanly() {
        let n = 128;
        let delay = 10.0;
        let mut fid = vec![Complex64::new(0.0, 0.0); n];
        for (i, slot) in fid.iter_mut().enumerate() {
            if i >= 10 {
                let t = (i - 10) as f64 * 0.05;
                *slot = Complex64::from_polar((-0.2 * t).exp(), 1.3 * t);
            }
        }
        let mut meta = AcquisitionParams::default();
        meta.group_delay = Some(delay);
        let sig = Signal1D::new(fid, Domain::Time, meta);
        let clean = remove_digital_filter(&sig, None);
        assert_eq!(clean.len(), n - 10);
        // Sample 0 should now match t = 0 (1.0 + 0i) closely
        assert!((clean.data[0].re - 1.0).abs() < 0.08);
        assert!(clean.data[0].im.abs() < 0.08);
    }
}
