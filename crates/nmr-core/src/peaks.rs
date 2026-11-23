use serde::{Deserialize, Serialize};

/// Detected NMR resonance peak with sub-bin interpolated position and width.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Peak {
    /// Sample index of the peak apex.
    pub index: usize,
    /// Chemical shift in ppm (interpolated).
    pub ppm: f64,
    /// Frequency relative to carrier in Hz.
    pub hz: f64,
    /// Peak intensity (real or magnitude).
    pub intensity: f64,
    /// Estimated full width at half maximum in Hz.
    pub fwhm_hz: Option<f64>,
    /// Signal-to-noise ratio of this peak.
    pub snr: f64,
}

/// Integrated chemical shift region [ppm_start, ppm_end].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntegralRegion {
    pub ppm_high: f64,
    pub ppm_low: f64,
    /// Raw sum of intensities inside region.
    pub raw_integral: f64,
    /// Integral normalized relative to the configured reference peak/region.
    pub normalized_integral: f64,
}

/// Collection of picked peaks and optional integral regions.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PeakList {
    pub peaks: Vec<Peak>,
    pub integrals: Vec<IntegralRegion>,
    pub noise_sigma: f64,
}
