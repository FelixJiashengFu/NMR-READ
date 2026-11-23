//! Core NMR data structures, spectral coordinate axes, peak models, and quality metrics.

pub mod metrics;
pub mod peaks;
pub mod signal;

pub use metrics::{QualityMetrics, compute_metrics, estimate_noise_sigma};
pub use peaks::{IntegralRegion, Peak, PeakList};
pub use signal::{
    AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D, Signal2D, SignalSet,
};
