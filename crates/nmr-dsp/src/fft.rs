use nmr_core::{Domain, Signal1D, Signal2D};
use num_complex::Complex64;
use rustfft::FftPlanner;

/// Subtract complex DC offset estimated from the trailing fraction (default last 15%) of the FID.
pub fn dc_offset_correction(sig: &Signal1D, tail_fraction: f64) -> Signal1D {
    if sig.is_empty() {
        return sig.clone();
    }
    let frac = tail_fraction.clamp(0.01, 0.5);
    let tail_len = ((sig.len() as f64) * frac).round().max(1.0) as usize;
    let start = sig.len().saturating_sub(tail_len);
    let sum: Complex64 = sig.data[start..].iter().copied().sum();
    let mean = sum / (sig.len() - start) as f64;
    let mut out = sig.clone();
    for v in &mut out.data {
        *v -= mean;
    }
    out
}

/// Zero-fill a 1D signal to `target_len` (or next power of two >= `len * multiplier` if `target_len == 0`).
pub fn zero_fill(sig: &Signal1D, target_len: usize, multiplier: usize) -> Signal1D {
    let desired = if target_len > sig.len() {
        target_len
    } else {
        (sig.len() * multiplier.max(1)).next_power_of_two()
    };
    let mut out = sig.clone();
    if desired > out.data.len() {
        out.data.resize(desired, Complex64::new(0.0, 0.0));
    }
    out
}

/// Forward 1D NMR Fourier transform (Time -> Frequency), with optional half-dwell first-point
/// scaling (`fid[0] *= 0.5`) and standard downfield-left-to-upfield-right spectral ordering.
pub fn fft_1d(sig: &Signal1D, half_first_point: bool) -> Signal1D {
    let n = sig.len();
    if n == 0 {
        return sig.clone();
    }
    let mut buf = sig.data.clone();
    if half_first_point && sig.domain == Domain::Time {
        buf[0] *= 0.5;
    }

    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);
    fft.process(&mut buf);

    let norm = 1.0 / (n as f64).sqrt();
    let half = n / 2;
    let mut ordered = vec![Complex64::new(0.0, 0.0); n];
    for (i, slot) in ordered.iter_mut().enumerate() {
        // Map index 0 -> +SW/2 (downfield left), half -> 0 Hz (carrier), n-1 -> -SW/2 (upfield right)
        let k = (half + n - i) % n;
        *slot = buf[k] * norm;
    }

    let mut out = sig.clone();
    out.data = ordered;
    out.domain = Domain::Frequency;
    out
}

/// Inverse 1D NMR Fourier transform (Frequency -> Time), inverting [`fft_1d`].
pub fn ifft_1d(sig: &Signal1D) -> Signal1D {
    let n = sig.len();
    if n == 0 {
        return sig.clone();
    }
    let half = n / 2;
    let mut buf = vec![Complex64::new(0.0, 0.0); n];
    for k in 0..n {
        let i = (half + n - k) % n;
        buf[k] = sig.data[i];
    }

    let mut planner = FftPlanner::<f64>::new();
    let ifft = planner.plan_fft_inverse(n);
    ifft.process(&mut buf);

    let norm = 1.0 / (n as f64).sqrt();
    for v in &mut buf {
        *v *= norm;
    }

    let mut out = sig.clone();
    out.data = buf;
    out.domain = Domain::Time;
    out
}

/// Perform 2D NMR FFT along direct F2 (rows) and optionally indirect F1 (columns).
pub fn fft_2d(sig: &Signal2D, transform_f1: bool) -> Signal2D {
    let rows = sig.rows;
    let cols = sig.cols;
    if rows == 0 || cols == 0 {
        return sig.clone();
    }

    let mut planner = FftPlanner::<f64>::new();
    let fft_row = planner.plan_fft_forward(cols);
    let norm_c = 1.0 / (cols as f64).sqrt();
    let half_c = cols / 2;

    let mut out = sig.clone();
    let mut scratch = vec![Complex64::new(0.0, 0.0); cols];
    for r in 0..rows {
        let row = out.row_mut(r);
        scratch.copy_from_slice(row);
        scratch[0] *= 0.5;
        fft_row.process(&mut scratch);
        for (i, slot) in row.iter_mut().enumerate() {
            let k = (half_c + cols - i) % cols;
            *slot = scratch[k] * norm_c;
        }
    }
    out.f2_domain = Domain::Frequency;

    if transform_f1 && rows > 1 {
        let fft_col = planner.plan_fft_forward(rows);
        let norm_r = 1.0 / (rows as f64).sqrt();
        let half_r = rows / 2;
        let mut col_buf = vec![Complex64::new(0.0, 0.0); rows];
        for c in 0..cols {
            for r in 0..rows {
                col_buf[r] = out.data[r * cols + c];
            }
            col_buf[0] *= 0.5;
            fft_col.process(&mut col_buf);
            for r in 0..rows {
                let k = (half_r + rows - r) % rows;
                out.data[r * cols + c] = col_buf[k] * norm_r;
            }
        }
        out.f1_domain = Domain::Frequency;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use nmr_core::AcquisitionParams;

    #[test]
    fn fft_ifft_roundtrip_exact() {
        let data: Vec<Complex64> = (0..64)
            .map(|i| Complex64::new((i as f64 * 0.2).cos(), (i as f64 * 0.3).sin()))
            .collect();
        let sig = Signal1D::new(data.clone(), Domain::Time, AcquisitionParams::default());
        let spec = fft_1d(&sig, false);
        assert_eq!(spec.domain, Domain::Frequency);
        let rec = ifft_1d(&spec);
        assert_eq!(rec.domain, Domain::Time);
        for (a, b) in data.iter().zip(rec.data.iter()) {
            assert!((*a - *b).norm() < 1e-11);
        }
    }
}
