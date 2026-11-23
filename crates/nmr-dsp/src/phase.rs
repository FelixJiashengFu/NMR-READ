use nmr_core::Signal1D;
use num_complex::Complex64;
use std::f64::consts::PI;

/// Apply 0th-order (`phi0_deg`) and 1st-order (`phi1_deg`) phase correction around `pivot_frac` (0..1).
pub fn phase_correct(sig: &Signal1D, phi0_deg: f64, phi1_deg: f64, pivot_frac: f64) -> Signal1D {
    let n = sig.len();
    if n == 0 {
        return sig.clone();
    }
    let p0 = phi0_deg * PI / 180.0;
    let p1 = phi1_deg * PI / 180.0;
    let piv = pivot_frac.clamp(0.0, 1.0);
    let denom = (n - 1).max(1) as f64;

    let mut out = sig.clone();
    for (i, c) in out.data.iter_mut().enumerate() {
        let frac = (i as f64 / denom) - piv;
        let angle = -(p0 + p1 * frac);
        *c *= Complex64::from_polar(1.0, angle);
    }
    out
}

/// Automatic 0th & 1st order phase correction via minimum-entropy & negative-lobe penalty
/// (Chen et al., *J. Magn. Reson.* 2002), returning `(corrected_signal, phi0_deg, phi1_deg)`.
pub fn auto_phase_entropy(sig: &Signal1D, allow_first_order: bool) -> (Signal1D, f64, f64) {
    let n = sig.len();
    if n < 8 {
        return (sig.clone(), 0.0, 0.0);
    }

    let score_fn = |p0: f64, p1: f64| -> f64 {
        let p0_rad = p0 * PI / 180.0;
        let p1_rad = p1 * PI / 180.0;
        let denom = (n - 1) as f64;

        let mut reals = Vec::with_capacity(n);
        for (i, c) in sig.data.iter().enumerate() {
            let angle = -(p0_rad + p1_rad * (i as f64 / denom - 0.5));
            let rot = c * Complex64::from_polar(1.0, angle);
            reals.push(rot.re);
        }

        // First derivative absolute sum for entropy
        let mut abs_d = Vec::with_capacity(n - 1);
        let mut sum_d = 0.0_f64;
        for w in reals.windows(2) {
            let d = (w[1] - w[0]).abs();
            abs_d.push(d);
            sum_d += d;
        }
        let sum_d = sum_d.max(1e-18);

        let mut entropy = 0.0_f64;
        for d in abs_d {
            let h = d / sum_d;
            if h > 1e-15 {
                entropy -= h * h.ln();
            }
        }

        // Penalize negative real lobes relative to total signal energy (scale- & length-invariant)
        let total_energy = reals.iter().map(|&v| v * v).sum::<f64>().max(1e-18);
        let neg_energy: f64 = reals.iter().filter(|&&v| v < 0.0).map(|&v| v * v).sum();
        let neg_ratio = neg_energy / total_energy;

        entropy + 1000.0 * neg_ratio
    };

    // Coarse grid search
    let mut best_p0 = 0.0;
    let mut best_p1 = 0.0;
    let mut best_score = f64::INFINITY;

    let p1_candidates: Vec<f64> = if allow_first_order {
        (-180..=180).step_by(30).map(|v| v as f64).collect()
    } else {
        vec![0.0]
    };

    for p0_i in (-180..180).step_by(15) {
        let p0 = p0_i as f64;
        for &p1 in &p1_candidates {
            let s = score_fn(p0, p1);
            if s < best_score {
                best_score = s;
                best_p0 = p0;
                best_p1 = p1;
            }
        }
    }

    // Multi-scale coordinate refinement
    for step in [8.0, 3.0, 1.0, 0.25] {
        let mut improved = true;
        let mut iters = 0;
        while improved && iters < 20 {
            improved = false;
            iters += 1;
            let dp1_opts: &[f64] = if allow_first_order { &[-step, 0.0, step] } else { &[0.0] };
            for dp0 in [-step, 0.0, step] {
                for &dp1 in dp1_opts {
                    if dp0 == 0.0 && dp1 == 0.0 {
                        continue;
                    }
                    let s = score_fn(best_p0 + dp0, best_p1 + dp1);
                    if s + 1e-9 < best_score {
                        best_score = s;
                        best_p0 += dp0;
                        best_p1 += dp1;
                        improved = true;
                    }
                }
            }
        }
    }

    // Normalize phi0 into [-180, 180)
    best_p0 = ((best_p0 + 180.0).rem_euclid(360.0)) - 180.0;
    let corrected = phase_correct(sig, best_p0, best_p1, 0.5);
    (corrected, best_p0, best_p1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nmr_core::{AcquisitionParams, Domain};

    #[test]
    fn auto_phase_recovers_zero_order_distortion() {
        let n = 256;
        // Pure absorptive Lorentzian peak at center rotated by +63 degrees
        let pure: Vec<Complex64> = (0..n)
            .map(|i| {
                let x = (i as f64 - 128.0) / 4.0;
                // Absorptive real + dispersive imag
                Complex64::new(1.0 / (1.0 + x * x), -x / (1.0 + x * x))
            })
            .collect();
        let sig = Signal1D::new(pure, Domain::Frequency, AcquisitionParams::default());
        let distorted = phase_correct(&sig, 63.0, 0.0, 0.5);
        let (fixed, p0, _p1) = auto_phase_entropy(&distorted, false);
        assert!((p0 + 63.0).abs() < 3.0, "Recovered p0={p0}");
        let min_re = fixed.data.iter().map(|c| c.re).fold(f64::INFINITY, f64::min);
        let max_re = fixed.data.iter().map(|c| c.re).fold(f64::NEG_INFINITY, f64::max);
        assert!(max_re > 0.9 && min_re > -0.1);
    }
}
