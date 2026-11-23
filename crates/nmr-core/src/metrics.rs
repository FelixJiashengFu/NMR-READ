use serde::{Deserialize, Serialize};

/// Quantitative signal quality & denoising comparison metrics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
    /// Estimated Signal-to-Noise Ratio in dB (based on robust noise sigma).
    pub snr_db: f64,
    /// Linear SNR (max signal / noise sigma).
    pub snr_linear: f64,
    /// Estimated noise standard deviation.
    pub noise_sigma: f64,
    /// Root Mean Square Error against reference (if provided).
    pub rmse: Option<f64>,
    /// Peak Signal-to-Noise Ratio in dB against reference (if provided).
    pub psnr_db: Option<f64>,
    /// 1D Structural Similarity Index (SSIM) against reference in [-1, 1] (if provided).
    pub ssim: Option<f64>,
}

/// Robustly estimate additive noise standard deviation using the Median Absolute Deviation (MAD)
/// of first-order differences (`diff / sqrt(2)`), which is resilient to sharp NMR peaks.
pub fn estimate_noise_sigma(values: &[f64]) -> f64 {
    if values.len() < 4 {
        return 0.0;
    }
    let mut diffs: Vec<f64> = values.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    diffs.sort_by(|a, b| a.total_cmp(b));
    let median_abs_diff = diffs[diffs.len() / 2];
    // For i.i.d. N(0, sigma^2), y[i+1] - y[i] ~ N(0, 2*sigma^2),
    // so MAD / 0.67448975 = sqrt(2) * sigma.
    median_abs_diff / (0.67448975 * std::f64::consts::SQRT_2)
}

/// Compute 1D Structural Similarity Index (SSIM) between `test` and `reference` signals.
///
/// Uses a sliding window of width `window_size` (clamped to array length).
pub fn ssim_1d(test: &[f64], reference: &[f64], window_size: usize) -> Option<f64> {
    if test.len() != reference.len() || test.is_empty() {
        return None;
    }
    let n = test.len();
    let max_abs = reference
        .iter()
        .chain(test.iter())
        .fold(0.0_f64, |acc, &v| acc.max(v.abs()))
        .max(1e-12);
    let c1 = (0.01 * max_abs).powi(2);
    let c2 = (0.03 * max_abs).powi(2);

    let win = window_size.max(3).min(n);
    let mut total = 0.0;
    let mut count = 0usize;

    for start in 0..=(n - win) {
        let x = &test[start..start + win];
        let y = &reference[start..start + win];
        let w = win as f64;
        let mean_x = x.iter().sum::<f64>() / w;
        let mean_y = y.iter().sum::<f64>() / w;
        let mut var_x = 0.0;
        let mut var_y = 0.0;
        let mut cov_xy = 0.0;
        for i in 0..win {
            let dx = x[i] - mean_x;
            let dy = y[i] - mean_y;
            var_x += dx * dx;
            var_y += dy * dy;
            cov_xy += dx * dy;
        }
        var_x /= w;
        var_y /= w;
        cov_xy /= w;

        let num = (2.0 * mean_x * mean_y + c1) * (2.0 * cov_xy + c2);
        let den = (mean_x * mean_x + mean_y * mean_y + c1) * (var_x + var_y + c2);
        if den > 0.0 {
            total += num / den;
            count += 1;
        }
    }

    if count == 0 {
        None
    } else {
        Some(total / count as f64)
    }
}

/// Compute single-signal SNR and optional reference-based metrics (RMSE, PSNR, SSIM).
pub fn compute_metrics(signal: &[f64], reference: Option<&[f64]>) -> QualityMetrics {
    let noise_sigma = estimate_noise_sigma(signal);
    let peak_amp = signal.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
    let snr_linear = if noise_sigma > 1e-15 {
        peak_amp / noise_sigma
    } else {
        f64::INFINITY
    };
    let snr_db = if snr_linear.is_finite() && snr_linear > 0.0 {
        20.0 * snr_linear.log10()
    } else {
        f64::INFINITY
    };

    let (rmse, psnr_db, ssim) = if let Some(ref_sig) = reference {
        if ref_sig.len() == signal.len() && !signal.is_empty() {
            let mse = signal
                .iter()
                .zip(ref_sig.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                / signal.len() as f64;
            let r = mse.sqrt();
            let ref_max = ref_sig.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
            let psnr = if r > 1e-15 && ref_max > 0.0 {
                Some(20.0 * (ref_max / r).log10())
            } else {
                None
            };
            let s = ssim_1d(signal, ref_sig, 31);
            (Some(r), psnr, s)
        } else {
            (None, None, None)
        }
    } else {
        (None, None, None)
    };

    QualityMetrics {
        snr_db,
        snr_linear,
        noise_sigma,
        rmse,
        psnr_db,
        ssim,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_signals_have_ssim_one_and_zero_rmse() {
        let a: Vec<f64> = (0..128).map(|i| (i as f64 * 0.1).sin()).collect();
        let m = compute_metrics(&a, Some(&a));
        assert!((m.rmse.unwrap() - 0.0).abs() < 1e-12);
        assert!((m.ssim.unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mad_noise_estimator_accurate_on_deterministic_sawtooth_noise() {
        // Constant zero signal + alternating +/- 0.5 has abs diff = 1.0
        let sig: Vec<f64> = (0..200).map(|i| if i % 2 == 0 { 0.5 } else { -0.5 }).collect();
        let sigma = estimate_noise_sigma(&sig);
        assert!(sigma > 0.5 && sigma < 1.5);
    }
}
