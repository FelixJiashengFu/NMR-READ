use faer::linalg::solvers::Solve;
use nmr_core::{Domain, Signal1D};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WindowType {
    Exponential,
    LorentzToGauss,
    SineBell,
    Hanning,
    Hamming,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ApodizationConfig {
    pub window: WindowType,
    /// Line broadening in Hz (positive = exponential decay, negative = Lorentz sharpening in GM).
    pub lb_hz: f64,
    /// Gaussian fraction [0.0..1.0] for Lorentz-to-Gauss window (maximum position as fraction of AQ).
    pub gb_frac: f64,
    /// Sine-bell shift (`0.0` = unshifted sin, `0.5` = cos-bell) and exponent power.
    pub sine_shift: f64,
    pub sine_power: f64,
}

impl Default for ApodizationConfig {
    fn default() -> Self {
        Self {
            window: WindowType::Exponential,
            lb_hz: 0.5,
            gb_frac: 0.15,
            sine_shift: 0.5,
            sine_power: 2.0,
        }
    }
}

/// Apply a time-domain apodization window to a FID, returning the modified signal
/// and the window envelope values (`Option<Vec<f64>>` stored for overlay plotting).
pub fn apply_apodization(sig: &Signal1D, cfg: &ApodizationConfig) -> (Signal1D, Vec<f64>) {
    let n = sig.len();
    if n == 0 {
        return (sig.clone(), Vec::new());
    }
    let dt = sig.meta.axis.dwell_time_s();
    let aq = ((n - 1) as f64 * dt).max(1e-9);
    let mut envelope = Vec::with_capacity(n);
    let mut out = sig.clone();

    for (i, c) in out.data.iter_mut().enumerate() {
        let t = i as f64 * dt;
        let frac = if n > 1 { i as f64 / (n - 1) as f64 } else { 0.0 };
        let w = match cfg.window {
            WindowType::Exponential => (-PI * cfg.lb_hz * t).exp(),
            WindowType::LorentzToGauss => {
                let t_max = cfg.gb_frac.clamp(0.01, 1.0) * aq;
                let a = PI * cfg.lb_hz.abs();
                let b = a / (2.0 * t_max);
                (a * t - b * t * t).exp()
            }
            WindowType::SineBell => {
                let s = cfg.sine_shift.clamp(0.0, 0.99);
                let arg = PI * (s + (1.0 - s) * frac);
                arg.sin().max(0.0).powf(cfg.sine_power.max(0.5))
            }
            WindowType::Hanning => 0.5 * (1.0 + (PI * frac).cos()),
            WindowType::Hamming => 0.54 + 0.46 * (PI * frac).cos(),
        };
        envelope.push(w);
        *c *= w;
    }

    (out, envelope)
}

/// Savitzky-Golay polynomial smoothing filter (for frequency-domain spectra or FIDs).
pub fn savitzky_golay_smooth(sig: &Signal1D, half_window: usize, poly_order: usize) -> Signal1D {
    let n = sig.len();
    if n < 5 || half_window == 0 {
        return sig.clone();
    }
    let hw = half_window.min((n - 1) / 2);
    let order = poly_order.min(hw * 2).min(5);
    let coeffs = compute_sg_coeffs(hw, order);

    let mut out = sig.clone();
    for i in 0..n {
        let mut acc = num_complex::Complex64::new(0.0, 0.0);
        for (k, &c) in coeffs.iter().enumerate() {
            let offset = k as isize - hw as isize;
            let idx = (i as isize + offset).clamp(0, (n - 1) as isize) as usize;
            acc += sig.data[idx] * c;
        }
        out.data[i] = acc;
    }
    if sig.domain == Domain::Frequency {
        out.domain = Domain::Frequency;
    }
    out
}

/// Compute 0th-derivative Savitzky-Golay FIR weights using normal equations via `faer`.
fn compute_sg_coeffs(hw: usize, order: usize) -> Vec<f64> {
    let m = 2 * hw + 1;
    let cols = order + 1;
    let mut a = faer::Mat::<f64>::zeros(m, cols);
    for r in 0..m {
        let x = (r as isize - hw as isize) as f64;
        let mut p = 1.0;
        for c in 0..cols {
            a[(r, c)] = p;
            p *= x;
        }
    }
    // Solve A^T A c = e_0 (first column of (A^T A)^{-1} A^T)
    let ata = a.transpose() * &a;
    let qr = ata.qr();
    let mut e0 = faer::Mat::<f64>::zeros(cols, 1);
    e0[(0, 0)] = 1.0;
    let sol = qr.solve(&e0);
    // c_j = sum_k sol[k] * A[j, k]
    (0..m)
        .map(|r| (0..cols).map(|k| a[(r, k)] * sol[(k, 0)]).sum())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use nmr_core::AcquisitionParams;
    use num_complex::Complex64;

    #[test]
    fn sg_preserves_quadratic_polynomial() {
        let data: Vec<Complex64> = (0..32)
            .map(|i| {
                let x = i as f64;
                Complex64::new(2.0 * x * x - 3.0 * x + 5.0, 0.0)
            })
            .collect();
        let sig = Signal1D::new(data, Domain::Frequency, AcquisitionParams::default());
        let sm = savitzky_golay_smooth(&sig, 3, 2);
        for i in 4..28 {
            assert!((sm.data[i].re - sig.data[i].re).abs() < 1e-8);
        }
    }
}
