use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Physical domain of an NMR signal array.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Domain {
    /// Time-domain Free Induction Decay (FID).
    Time,
    /// Frequency-domain spectrum (after Fourier transform).
    Frequency,
}

/// Source format identifier for provenance and metadata tracking.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum NmrFormat {
    Synthetic,
    Bruker,
    Varian,
    Jeol,
    Spinsolve,
    NmrPipe,
    SparkyUcsf,
    JcampDx,
    NmrMl,
    Csv,
}

/// Spectral axis calibration parameters for one NMR dimension.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AxisMeta {
    /// Nucleus label (e.g. "1H", "13C", "15N").
    pub nucleus: String,
    /// Spectrometer observing frequency in MHz (SFO1 / BF1).
    pub spectrometer_freq_mhz: f64,
    /// Spectral width in Hz.
    pub spectral_width_hz: f64,
    /// Downfield (leftmost) chemical shift edge in ppm (Bruker OFFSET convention), if known.
    pub offset_ppm: Option<f64>,
    /// Carrier/center chemical shift in ppm (fallback when `offset_ppm` is absent).
    pub carrier_ppm: f64,
}

impl Default for AxisMeta {
    fn default() -> Self {
        Self {
            nucleus: "1H".to_string(),
            spectrometer_freq_mhz: 500.0,
            spectral_width_hz: 6000.0,
            offset_ppm: None,
            carrier_ppm: 5.0,
        }
    }
}

impl AxisMeta {
    /// Total spectral width in ppm.
    pub fn sw_ppm(&self) -> f64 {
        if self.spectrometer_freq_mhz > 0.0 {
            self.spectral_width_hz / self.spectrometer_freq_mhz
        } else {
            10.0
        }
    }

    /// Dwell time between consecutive time-domain points in seconds.
    pub fn dwell_time_s(&self) -> f64 {
        if self.spectral_width_hz > 0.0 {
            1.0 / self.spectral_width_hz
        } else {
            1e-4
        }
    }

    /// Leftmost (highest/downfield) ppm coordinate.
    pub fn left_ppm(&self) -> f64 {
        self.offset_ppm
            .unwrap_or_else(|| self.carrier_ppm + 0.5 * self.sw_ppm())
    }

    /// Rightmost (lowest/upfield) ppm coordinate.
    pub fn right_ppm(&self) -> f64 {
        self.left_ppm() - self.sw_ppm()
    }

    /// Chemical shift (ppm) at index `i` for a frequency spectrum of length `n`
    /// stored in standard high-to-low (downfield-to-upfield) plotting order.
    pub fn index_to_ppm(&self, i: usize, n: usize) -> f64 {
        if n <= 1 {
            return self.left_ppm();
        }
        let frac = i as f64 / (n - 1) as f64;
        self.left_ppm() - frac * self.sw_ppm()
    }

    /// Frequency (Hz, relative to carrier) at index `i` for a spectrum of length `n`.
    pub fn index_to_hz(&self, i: usize, n: usize) -> f64 {
        if n <= 1 {
            return 0.5 * self.spectral_width_hz;
        }
        let frac = i as f64 / (n - 1) as f64;
        (0.5 - frac) * self.spectral_width_hz
    }

    /// Time in seconds at FID sample index `i`.
    pub fn index_to_time_s(&self, i: usize) -> f64 {
        i as f64 * self.dwell_time_s()
    }

    /// Build the full chemical shift axis vector (ppm, descending left-to-right).
    pub fn ppm_axis(&self, n: usize) -> Vec<f64> {
        (0..n).map(|i| self.index_to_ppm(i, n)).collect()
    }

    /// Build the full time axis vector in seconds.
    pub fn time_axis_s(&self, n: usize) -> Vec<f64> {
        (0..n).map(|i| self.index_to_time_s(i)).collect()
    }
}

/// Acquisition and vendor metadata associated with an NMR dataset.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AcquisitionParams {
    pub title: String,
    pub format: NmrFormat,
    pub axis: AxisMeta,
    /// Indirect axis metadata (for 2D datasets).
    pub indirect_axis: Option<AxisMeta>,
    /// Bruker digital filter group delay (points) if present and unremoved.
    pub group_delay: Option<f64>,
    /// Number of scans accumulated.
    pub num_scans: usize,
    /// Raw vendor parameter dictionary for inspection.
    pub vendor_params: BTreeMap<String, String>,
}

impl Default for AcquisitionParams {
    fn default() -> Self {
        Self {
            title: "Untitled NMR".to_string(),
            format: NmrFormat::Synthetic,
            axis: AxisMeta::default(),
            indirect_axis: None,
            group_delay: None,
            num_scans: 1,
            vendor_params: BTreeMap::new(),
        }
    }
}

/// Complex 1D NMR signal in either time (FID) or frequency (spectrum) domain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Signal1D {
    pub data: Vec<Complex64>,
    pub domain: Domain,
    pub meta: AcquisitionParams,
    /// Whether Bruker/vendor digital filter group delay has already been removed.
    pub digital_filter_removed: bool,
    /// Estimated baseline array (same length as `data`), if baseline subtraction was applied or computed.
    pub baseline: Option<Vec<f64>>,
}

impl Signal1D {
    pub fn new(data: Vec<Complex64>, domain: Domain, meta: AcquisitionParams) -> Self {
        let digital_filter_removed = meta.group_delay.map_or(true, |d| d <= 0.0);
        Self {
            data,
            domain,
            meta,
            digital_filter_removed,
            baseline: None,
        }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn real(&self) -> Vec<f64> {
        self.data.iter().map(|c| c.re).collect()
    }

    pub fn imag(&self) -> Vec<f64> {
        self.data.iter().map(|c| c.im).collect()
    }

    pub fn magnitude(&self) -> Vec<f64> {
        self.data.iter().map(|c| c.norm()).collect()
    }

    pub fn x_axis_ppm(&self) -> Vec<f64> {
        self.meta.axis.ppm_axis(self.data.len())
    }

    pub fn x_axis_time_s(&self) -> Vec<f64> {
        self.meta.axis.time_axis_s(self.data.len())
    }
}

/// Collection of 1D scans (for multi-scan averaging, denoising, or pseudo-2D series).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SignalSet {
    pub scans: Vec<Signal1D>,
}

impl SignalSet {
    pub fn new(scans: Vec<Signal1D>) -> Self {
        Self { scans }
    }

    pub fn len(&self) -> usize {
        self.scans.len()
    }

    pub fn is_empty(&self) -> bool {
        self.scans.is_empty()
    }
}

/// Complex 2D NMR matrix stored in row-major order (`rows` indirect F1 × `cols` direct F2).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Signal2D {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Complex64>,
    pub f2_domain: Domain,
    pub f1_domain: Domain,
    pub meta: AcquisitionParams,
}

impl Signal2D {
    pub fn new(
        rows: usize,
        cols: usize,
        data: Vec<Complex64>,
        f2_domain: Domain,
        f1_domain: Domain,
        meta: AcquisitionParams,
    ) -> Self {
        assert_eq!(rows * cols, data.len(), "2D matrix dimension mismatch");
        Self {
            rows,
            cols,
            data,
            f2_domain,
            f1_domain,
            meta,
        }
    }

    pub fn row(&self, r: usize) -> &[Complex64] {
        let start = r * self.cols;
        &self.data[start..start + self.cols]
    }

    pub fn row_mut(&mut self, r: usize) -> &mut [Complex64] {
        let start = r * self.cols;
        &mut self.data[start..start + self.cols]
    }

    /// Real values in row-major order (`rows * cols`).
    pub fn real_matrix(&self) -> Vec<f64> {
        self.data.iter().map(|c| c.re).collect()
    }

    /// Magnitude values in row-major order (`rows * cols`).
    pub fn magnitude_matrix(&self) -> Vec<f64> {
        self.data.iter().map(|c| c.norm()).collect()
    }
}
