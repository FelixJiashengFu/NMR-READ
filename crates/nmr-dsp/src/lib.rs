//! Pure-Rust NMR signal processing, denoising, peak picking, and synthetic FID generation.
//!
//! Zero C/Fortran dependencies: runs identically on native Desktop (macOS/Linux/Windows)
//! and browser WebAssembly (`wasm32-unknown-unknown`).

pub mod apodization;
pub mod baseline;
pub mod cadzow;
pub mod digital_filter;
pub mod fft;
pub mod multiscan;
pub mod peaks;
pub mod phase;
pub mod synthetic;
pub mod wavelet;

pub use apodization::{ApodizationConfig, WindowType, apply_apodization, savitzky_golay_smooth};
pub use baseline::{BaselineMethod, arpls_baseline, polynomial_baseline, subtract_baseline};
pub use cadzow::{CadzowConfig, cadzow_denoise};
pub use digital_filter::remove_digital_filter;
pub use fft::{dc_offset_correction, fft_1d, fft_2d, ifft_1d, zero_fill};
pub use multiscan::{add_complex_noise, average_scans, wavelet_denoise_scans};
pub use peaks::{integrate_region, pick_peaks};
pub use phase::{auto_phase_entropy, phase_correct};
pub use synthetic::{MultipletPeak, SyntheticPreset, generate_synthetic_2d_cosy, generate_synthetic_fid};
pub use wavelet::{ThresholdMode, WaveletFamily, WaveletConfig, wavelet_denoise};
