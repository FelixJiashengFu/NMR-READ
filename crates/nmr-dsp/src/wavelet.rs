use nmr_core::{Signal1D, estimate_noise_sigma};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::FRAC_1_SQRT_2;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WaveletFamily {
    Haar,
    Db2,
    Db4,
    Db6,
    Sym4,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ThresholdMode {
    Soft,
    Hard,
    Garrote,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaveletConfig {
    pub family: WaveletFamily,
    pub levels: usize,
    pub mode: ThresholdMode,
    /// Multiplier on universal VisuShrink threshold `sigma * sqrt(2 * ln(N))`.
    pub threshold_scale: f64,
}

impl Default for WaveletConfig {
    fn default() -> Self {
        Self {
            family: WaveletFamily::Sym4,
            levels: 4,
            mode: ThresholdMode::Soft,
            threshold_scale: 1.0,
        }
    }
}

/// Return orthonormal analysis low-pass filter `h_0` (decomposition, `||h_0||_2 = 1`).
fn low_pass_filter(family: WaveletFamily) -> &'static [f64] {
    match family {
        WaveletFamily::Haar => &[FRAC_1_SQRT_2, FRAC_1_SQRT_2],
        WaveletFamily::Db2 => &[
            -0.12940952255092145,
            0.22414386804185735,
            0.836516303737469,
            0.48296291314469025,
        ],
        WaveletFamily::Db4 => &[
            -0.010597401784997278,
            0.032883011666982945,
            0.030841381835986965,
            -0.18703481171888114,
            -0.02798376941698385,
            0.6308807679295904,
            0.7148465705525415,
            0.23037781330885523,
        ],
        WaveletFamily::Db6 => &[
            -0.00107730108499558,
            0.004777257510910411,
            0.0005538422009938016,
            -0.031582039318031156,
            0.02752286553001629,
            0.09750160558707936,
            -0.12976686756709563,
            -0.22626469396516913,
            0.3152503517092432,
            0.7511339080215775,
            0.4946238903983854,
            0.11154074335008017,
        ],
        WaveletFamily::Sym4 => &[
            -0.07576571478927333,
            -0.02963552764599851,
            0.49761866763201545,
            0.8037387518059161,
            0.29785779560527736,
            -0.09921954357684722,
            -0.012603967262037833,
            0.0322231006040427,
        ],
    }
}

/// Derive orthogonal quadrature mirror high-pass filter `h_1[k] = (-1)^k h_0[L - 1 - k]`.
fn high_pass_from_low(lo: &[f64]) -> Vec<f64> {
    let l = lo.len();
    (0..l)
        .map(|k| {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            sign * lo[l - 1 - k]
        })
        .collect()
}

/// 1D Stationary Wavelet Transform (à trous, periodic boundaries) returning `(approx_final, details_per_level)`.
pub fn swt_1d(signal: &[f64], family: WaveletFamily, levels: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = signal.len();
    if n == 0 || levels == 0 {
        return (signal.to_vec(), Vec::new());
    }
    let lo = low_pass_filter(family);
    let hi = high_pass_from_low(lo);
    let m = lo.len();

    let mut approx = signal.to_vec();
    let mut details = Vec::with_capacity(levels);

    for level in 0..levels {
        let step = 1usize << level;
        let mut next_a = vec![0.0; n];
        let mut d = vec![0.0; n];
        for i in 0..n {
            let mut sa = 0.0;
            let mut sd = 0.0;
            for k in 0..m {
                // Periodic wrap indexing (works cleanly for any signal length N)
                let idx = (i + k * step) % n;
                let v = approx[idx];
                sa += lo[k] * v;
                sd += hi[k] * v;
            }
            next_a[i] = sa;
            d[i] = sd;
        }
        approx = next_a;
        details.push(d);
    }

    (approx, details)
}

/// Inverse 1D Stationary Wavelet Transform (exact dual-shift average reconstruction).
pub fn iswt_1d(approx_final: &[f64], details: &[Vec<f64>], family: WaveletFamily) -> Vec<f64> {
    let n = approx_final.len();
    if n == 0 || details.is_empty() {
        return approx_final.to_vec();
    }
    let lo = low_pass_filter(family);
    let hi = high_pass_from_low(lo);
    let m = lo.len();

    let mut approx = approx_final.to_vec();
    for level in (0..details.len()).rev() {
        let step = 1usize << level;
        let d = &details[level];
        let mut rec = vec![0.0; n];
        // Average over the two cosets (equivalent to inverse of redundant undecimated filterbank)
        for i in 0..n {
            let mut sum = 0.0;
            for k in 0..m {
                let idx = (i + n - ((k * step) % n)) % n;
                sum += lo[k] * approx[idx] + hi[k] * d[idx];
            }
            rec[i] = 0.5 * sum;
        }
        approx = rec;
    }
    approx
}

fn apply_threshold(val: f64, thresh: f64, mode: ThresholdMode) -> f64 {
    let abs = val.abs();
    if abs <= thresh {
        return 0.0;
    }
    match mode {
        ThresholdMode::Hard => val,
        ThresholdMode::Soft => val.signum() * (abs - thresh),
        ThresholdMode::Garrote => val - (thresh * thresh) / val,
    }
}

fn denoise_real_swt(vals: &[f64], cfg: &WaveletConfig) -> Vec<f64> {
    let n = vals.len();
    if n < 8 {
        return vals.to_vec();
    }
    let levels = cfg.levels.clamp(1, 8);
    // Symmetric mirror extension [x_0..x_{N-1}, x_{N-1}..x_0] eliminates boundary wrap
    // between high-amplitude t=0 FID origin and noisy tail while preserving exact ISWT.
    let mut padded = Vec::with_capacity(2 * n);
    padded.extend_from_slice(vals);
    padded.extend(vals.iter().rev().copied());

    let (approx, mut details) = swt_1d(&padded, cfg.family, levels);

    let sigma = estimate_noise_sigma(vals).max(1e-15);
    let univ = sigma * (2.0 * (n as f64).ln()).sqrt() * cfg.threshold_scale.max(0.0);

    for d_band in details.iter_mut() {
        for v in d_band.iter_mut() {
            *v = apply_threshold(*v, univ, cfg.mode);
        }
    }

    let mut rec = iswt_1d(&approx, &details, cfg.family);
    rec.truncate(n);
    rec
}

/// Denoise a complex 1D NMR signal using Stationary Wavelet Transform thresholding.
pub fn wavelet_denoise(sig: &Signal1D, cfg: &WaveletConfig) -> Signal1D {
    if sig.len() < 8 {
        return sig.clone();
    }
    let re = denoise_real_swt(&sig.real(), cfg);
    let im = denoise_real_swt(&sig.imag(), cfg);
    let mut out = sig.clone();
    out.data = re
        .into_iter()
        .zip(im)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn swt_iswt_perfect_reconstruction_across_families() {
        let sig: Vec<f64> = (0..128)
            .map(|i| (i as f64 * 0.13).sin() + 0.5 * (i as f64 * 0.41).cos())
            .collect();
        for fam in [
            WaveletFamily::Haar,
            WaveletFamily::Db2,
            WaveletFamily::Db4,
            WaveletFamily::Db6,
            WaveletFamily::Sym4,
        ] {
            let (a, d) = swt_1d(&sig, fam, 4);
            let rec = iswt_1d(&a, &d, fam);
            let max_err = sig
                .iter()
                .zip(rec.iter())
                .map(|(x, y)| (x - y).abs())
                .fold(0.0_f64, f64::max);
            assert!(max_err < 1e-10, "Family {fam:?} max_err={max_err}");
        }
    }
}
