use crate::wavelet::{WaveletConfig, wavelet_denoise};
use nmr_core::{Signal1D, SignalSet};
use num_complex::Complex64;
use std::f64::consts::TAU;

/// Deterministic lightweight PRNG (zero OS/Wasm entropy dependencies).
pub(crate) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub(crate) fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(0x9E37_79B9_7F4A_7C15),
        }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub(crate) fn next_unit_f64(&mut self) -> f64 {
        // Uniform in (0, 1]
        let v = (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64);
        v.max(1e-12)
    }

    pub(crate) fn next_gaussian_pair(&mut self) -> (f64, f64) {
        let u1 = self.next_unit_f64();
        let u2 = self.next_unit_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = TAU * u2;
        (r * theta.cos(), r * theta.sin())
    }
}

/// Add complex Gaussian noise `N(0, sigma^2)` independently to real and imaginary channels.
pub fn add_complex_noise(sig: &Signal1D, sigma: f64, seed: u64) -> Signal1D {
    let mut rng = SplitMix64::new(seed);
    let mut out = sig.clone();
    for c in &mut out.data {
        let (nr, ni) = rng.next_gaussian_pair();
        *c += Complex64::new(nr * sigma, ni * sigma);
    }
    out
}

/// Generate a multi-scan [`SignalSet`] by adding independent noise realizations to a base scan.
pub fn generate_noisy_scans(base: &Signal1D, num_scans: usize, sigma: f64, base_seed: u64) -> SignalSet {
    let scans = (0..num_scans.max(1))
        .map(|k| add_complex_noise(base, sigma, base_seed.wrapping_add(k as u64 * 104729)))
        .collect();
    SignalSet::new(scans)
}

/// Average all scans in a [`SignalSet`] into a single [`Signal1D`].
pub fn average_scans(set: &SignalSet) -> Option<Signal1D> {
    let first = set.scans.first()?;
    let n = first.len();
    let count = set.scans.len() as f64;
    let mut acc = vec![Complex64::new(0.0, 0.0); n];
    for scan in &set.scans {
        for (dst, src) in acc.iter_mut().zip(scan.data.iter()) {
            *dst += *src;
        }
    }
    for v in &mut acc {
        *v /= count;
    }
    let mut out = first.clone();
    out.data = acc;
    out.meta.num_scans = set.scans.len();
    Some(out)
}

/// Apply Stationary Wavelet Denoising independently to every scan in a [`SignalSet`].
pub fn wavelet_denoise_scans(set: &SignalSet, cfg: &WaveletConfig) -> SignalSet {
    let scans = set
        .scans
        .iter()
        .map(|s| wavelet_denoise(s, cfg))
        .collect();
    SignalSet::new(scans)
}
