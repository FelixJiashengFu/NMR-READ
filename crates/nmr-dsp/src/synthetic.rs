use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D, Signal2D};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::TAU;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SyntheticPreset {
    /// 1H Ethanol (TMS singlet + CH3 triplet + CH2 quartet + OH singlet) at 500 MHz.
    Ethanol,
    /// 1H Ethylbenzene (TMS + CH3 triplet + CH2 quartet + aromatic multiplets) at 500 MHz.
    Ethylbenzene,
    /// Simple 2-peak diagnostic Lorentzian pair.
    TwoPeaks,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultipletPeak {
    pub shift_ppm: f64,
    pub amplitude: f64,
    pub t2_seconds: f64,
    /// J-coupling splitting in Hz (0.0 for singlet).
    pub j_hz: f64,
    /// Number of equivalent coupled neighbors (`0` = singlet, `2` = triplet, `3` = quartet).
    pub multiplicity_n: usize,
}

pub fn preset_multiplets(preset: SyntheticPreset) -> Vec<MultipletPeak> {
    match preset {
        SyntheticPreset::Ethanol => vec![
            // TMS reference
            MultipletPeak {
                shift_ppm: 0.0,
                amplitude: 0.5,
                t2_seconds: 0.65,
                j_hz: 0.0,
                multiplicity_n: 0,
            },
            // CH3 triplet (3H, 1:2:1)
            MultipletPeak {
                shift_ppm: 1.18,
                amplitude: 3.0,
                t2_seconds: 0.45,
                j_hz: 7.0,
                multiplicity_n: 2,
            },
            // CH2 quartet (2H, 1:3:3:1)
            MultipletPeak {
                shift_ppm: 3.65,
                amplitude: 2.0,
                t2_seconds: 0.45,
                j_hz: 7.0,
                multiplicity_n: 3,
            },
            // OH broad singlet (1H)
            MultipletPeak {
                shift_ppm: 4.82,
                amplitude: 1.0,
                t2_seconds: 0.22,
                j_hz: 0.0,
                multiplicity_n: 0,
            },
        ],
        SyntheticPreset::Ethylbenzene => vec![
            MultipletPeak {
                shift_ppm: 0.0,
                amplitude: 0.4,
                t2_seconds: 0.65,
                j_hz: 0.0,
                multiplicity_n: 0,
            },
            MultipletPeak {
                shift_ppm: 1.24,
                amplitude: 3.0,
                t2_seconds: 0.45,
                j_hz: 7.6,
                multiplicity_n: 2,
            },
            MultipletPeak {
                shift_ppm: 2.65,
                amplitude: 2.0,
                t2_seconds: 0.45,
                j_hz: 7.6,
                multiplicity_n: 3,
            },
            MultipletPeak {
                shift_ppm: 7.18,
                amplitude: 2.5,
                t2_seconds: 0.35,
                j_hz: 7.5,
                multiplicity_n: 2,
            },
            MultipletPeak {
                shift_ppm: 7.28,
                amplitude: 2.5,
                t2_seconds: 0.35,
                j_hz: 7.5,
                multiplicity_n: 1,
            },
        ],
        SyntheticPreset::TwoPeaks => vec![
            MultipletPeak {
                shift_ppm: 2.15,
                amplitude: 1.8,
                t2_seconds: 0.4,
                j_hz: 0.0,
                multiplicity_n: 0,
            },
            MultipletPeak {
                shift_ppm: 7.26,
                amplitude: 1.2,
                t2_seconds: 0.5,
                j_hz: 0.0,
                multiplicity_n: 0,
            },
        ],
    }
}

/// Generate a realistic synthetic 1D complex NMR FID using Pascal-triangle J-coupling splitting.
pub fn generate_synthetic_fid(
    preset: SyntheticPreset,
    num_points: usize,
    sf_mhz: f64,
    sw_ppm: f64,
    carrier_ppm: f64,
    noise_sigma: f64,
    seed: u64,
) -> Signal1D {
    let n = num_points.max(64);
    let sw_hz = sw_ppm * sf_mhz;
    let dt = 1.0 / sw_hz;
    let multiplets = preset_multiplets(preset);

    // Precompute Pascal multiplet lines (offset_hz, rel_amp, t2)
    let mut lines = Vec::new();
    for m in &multiplets {
        let base_hz = (m.shift_ppm - carrier_ppm) * sf_mhz;
        if m.multiplicity_n == 0 || m.j_hz <= 0.0 {
            lines.push((base_hz, m.amplitude, m.t2_seconds));
        } else {
            let weights = pascal_row(m.multiplicity_n);
            let total_w: f64 = weights.iter().sum();
            let k_max = m.multiplicity_n as f64;
            for (idx, w) in weights.into_iter().enumerate() {
                let offset = ((idx as f64) - 0.5 * k_max) * m.j_hz;
                lines.push((base_hz + offset, m.amplitude * (w / total_w), m.t2_seconds));
            }
        }
    }

    let mut data = vec![Complex64::new(0.0, 0.0); n];
    for (i, slot) in data.iter_mut().enumerate() {
        let t = i as f64 * dt;
        let mut sum = Complex64::new(0.0, 0.0);
        for &(freq_hz, amp, t2) in &lines {
            let decay = (-t / t2.max(1e-4)).exp();
            sum += Complex64::from_polar(amp * decay, TAU * freq_hz * t);
        }
        *slot = sum;
    }

    let meta = AcquisitionParams {
        title: format!("Synthetic {preset:?}"),
        format: NmrFormat::Synthetic,
        axis: AxisMeta {
            nucleus: "1H".to_string(),
            spectrometer_freq_mhz: sf_mhz,
            spectral_width_hz: sw_hz,
            offset_ppm: None,
            carrier_ppm,
        },
        indirect_axis: None,
        group_delay: None,
        num_scans: 1,
        vendor_params: Default::default(),
    };

    let sig = Signal1D::new(data, Domain::Time, meta);
    if noise_sigma > 0.0 {
        crate::multiscan::add_complex_noise(&sig, noise_sigma, seed)
    } else {
        sig
    }
}

fn pascal_row(n: usize) -> Vec<f64> {
    let mut row = vec![1.0];
    for _ in 0..n {
        let mut next = vec![1.0; row.len() + 1];
        for i in 1..row.len() {
            next[i] = row[i - 1] + row[i];
        }
        row = next;
    }
    row
}

/// Generate a synthetic 2D COSY-like time-domain matrix (`rows × cols`) for 2D FFT and heatmap visualization.
pub fn generate_synthetic_2d_cosy(rows: usize, cols: usize, noise_sigma: f64, seed: u64) -> Signal2D {
    let r_n = rows.max(16);
    let c_n = cols.max(16);
    let sf = 500.0;
    let sw_ppm = 10.0;
    let sw_hz = sf * sw_ppm;
    let dt = 1.0 / sw_hz;
    let carrier = 5.0;

    // Coupled resonances at 1.18 ppm and 3.65 ppm (diagonal + cross peaks)
    let f1_a = (1.18 - carrier) * sf;
    let f1_b = (3.65 - carrier) * sf;
    let pairs = [
        (f1_a, f1_a, 1.5),
        (f1_b, f1_b, 1.5),
        (f1_a, f1_b, 0.85),
        (f1_b, f1_a, 0.85),
    ];

    let mut data = vec![Complex64::new(0.0, 0.0); r_n * c_n];
    let mut rng = crate::multiscan::SplitMix64::new(seed);
    for r in 0..r_n {
        let t1 = r as f64 * dt;
        for c in 0..c_n {
            let t2 = c as f64 * dt;
            let mut val = Complex64::new(0.0, 0.0);
            for &(fa, fb, amp) in &pairs {
                let env = amp * (-t1 / 0.12).exp() * (-t2 / 0.12).exp();
                val += Complex64::from_polar(env, TAU * (fa * t1 + fb * t2));
            }
            if noise_sigma > 0.0 {
                let (nr, ni) = rng.next_gaussian_pair();
                val += Complex64::new(nr * noise_sigma, ni * noise_sigma);
            }
            data[r * c_n + c] = val;
        }
    }

    let axis = AxisMeta {
        nucleus: "1H".into(),
        spectrometer_freq_mhz: sf,
        spectral_width_hz: sw_hz,
        offset_ppm: None,
        carrier_ppm: carrier,
    };
    let meta = AcquisitionParams {
        title: "Synthetic 2D COSY".into(),
        format: NmrFormat::Synthetic,
        axis: axis.clone(),
        indirect_axis: Some(axis),
        group_delay: None,
        num_scans: 1,
        vendor_params: Default::default(),
    };

    Signal2D::new(r_n, c_n, data, Domain::Time, Domain::Time, meta)
}
