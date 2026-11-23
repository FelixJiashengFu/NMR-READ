use nmr_core::{Domain, IntegralRegion, Peak, PeakList, Signal1D, estimate_noise_sigma};

/// Detect resonance peaks in a frequency-domain spectrum (or magnitude array) with
/// parabolic sub-bin interpolation and FWHM estimation.
pub fn pick_peaks(
    sig: &Signal1D,
    snr_threshold: f64,
    min_separation_pts: usize,
    use_magnitude: bool,
) -> PeakList {
    let vals = if use_magnitude {
        sig.magnitude()
    } else {
        sig.real()
    };
    let n = vals.len();
    if n < 3 {
        return PeakList::default();
    }

    let noise_sigma = estimate_noise_sigma(&vals).max(1e-15);
    let min_height = snr_threshold.max(1.0) * noise_sigma;
    let sep = min_separation_pts.max(1);

    let mut candidates = Vec::new();
    for i in 1..(n - 1) {
        let y0 = vals[i - 1];
        let y1 = vals[i];
        let y2 = vals[i + 1];
        if y1 >= min_height && y1 > y0 && y1 >= y2 {
            // Parabolic vertex offset in [-0.5, 0.5]
            let denom = 2.0 * (y0 - 2.0 * y1 + y2);
            let delta = if denom.abs() > 1e-18 {
                ((y0 - y2) / denom).clamp(-0.5, 0.5)
            } else {
                0.0
            };
            let interp_idx = i as f64 + delta;
            let ppm = if sig.domain == Domain::Frequency {
                let frac = interp_idx / (n - 1).max(1) as f64;
                sig.meta.axis.left_ppm() - frac * sig.meta.axis.sw_ppm()
            } else {
                interp_idx
            };
            let hz = if sig.domain == Domain::Frequency {
                let frac = interp_idx / (n - 1).max(1) as f64;
                (0.5 - frac) * sig.meta.axis.spectral_width_hz
            } else {
                0.0
            };

            let fwhm_hz = estimate_fwhm_hz(&vals, i, sig.meta.axis.spectral_width_hz);
            candidates.push(Peak {
                index: i,
                ppm,
                hz,
                intensity: y1,
                fwhm_hz,
                snr: y1 / noise_sigma,
            });
        }
    }

    // Non-maximum suppression within `sep` points (sort highest first)
    candidates.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
    let mut kept: Vec<Peak> = Vec::new();
    for c in candidates {
        if kept
            .iter()
            .all(|k| k.index.abs_diff(c.index) >= sep)
        {
            kept.push(c);
        }
    }
    kept.sort_by(|a, b| a.index.cmp(&b.index));

    PeakList {
        peaks: kept,
        integrals: Vec::new(),
        noise_sigma,
    }
}

fn estimate_fwhm_hz(vals: &[f64], apex: usize, sw_hz: f64) -> Option<f64> {
    let half_max = 0.5 * vals[apex];
    if half_max <= 0.0 || vals.len() < 3 {
        return None;
    }
    let mut left = apex as f64;
    for i in (0..apex).rev() {
        if vals[i] <= half_max {
            let denom = (vals[i + 1] - vals[i]).max(1e-15);
            let t = (half_max - vals[i]) / denom;
            left = i as f64 + t;
            break;
        }
    }
    let mut right = apex as f64;
    for i in (apex + 1)..vals.len() {
        if vals[i] <= half_max {
            let denom = (vals[i - 1] - vals[i]).max(1e-15);
            let t = (vals[i - 1] - half_max) / denom;
            right = (i - 1) as f64 + t;
            break;
        }
    }
    let width_pts = (right - left).max(0.5);
    let hz_per_pt = sw_hz / (vals.len() - 1).max(1) as f64;
    Some(width_pts * hz_per_pt)
}

/// Integrate real spectrum intensity within `[ppm_low, ppm_high]`.
pub fn integrate_region(
    sig: &Signal1D,
    ppm_low: f64,
    ppm_high: f64,
    reference_scale: f64,
) -> IntegralRegion {
    let lo = ppm_low.min(ppm_high);
    let hi = ppm_low.max(ppm_high);
    let x = sig.x_axis_ppm();
    let sum: f64 = x
        .iter()
        .zip(sig.data.iter())
        .filter(|(p, _)| **p >= lo && **p <= hi)
        .map(|(_, c)| c.re)
        .sum();
    IntegralRegion {
        ppm_high: hi,
        ppm_low: lo,
        raw_integral: sum,
        normalized_integral: if reference_scale.abs() > 1e-15 {
            sum / reference_scale
        } else {
            sum
        },
    }
}
