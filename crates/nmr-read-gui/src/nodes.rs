use nmr_core::{PeakList, QualityMetrics, Signal1D, Signal2D};
use nmr_dsp::{
    ApodizationConfig, BaselineMethod, CadzowConfig, SyntheticPreset, ThresholdMode,
    WaveletConfig, WaveletFamily, WindowType,
};
use serde::{Deserialize, Serialize};

/// Semantic data type carried across node graph wires.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PinKind {
    Signal1D,
    Signal2D,
    PeakList,
}

impl PinKind {
    pub fn color(self) -> egui::Color32 {
        match self {
            PinKind::Signal1D => egui::Color32::from_rgb(0x4f, 0xc3, 0xf7), // Cyan-blue
            PinKind::Signal2D => egui::Color32::from_rgb(0xab, 0x47, 0xbc), // Purple
            PinKind::PeakList => egui::Color32::from_rgb(0xff, 0xb7, 0x4d), // Amber
        }
    }
}

/// Serializable node definitions for the NMR-READ processing graph.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum NmrNode {
    // --- Sources ---
    Synthetic1D {
        preset: SyntheticPreset,
        points: usize,
        sf_mhz: f64,
        sw_ppm: f64,
        carrier_ppm: f64,
        noise_sigma: f64,
        seed: u64,
    },
    Synthetic2DCosy {
        rows: usize,
        cols: usize,
        noise_sigma: f64,
        seed: u64,
    },
    LoadedDataset1D {
        dataset_key: String,
    },
    LoadedDataset2D {
        dataset_key: String,
    },

    // --- Conditioning & Time Domain ---
    RemoveDigitalFilter {
        use_override: bool,
        grpdly_override: f64,
    },
    DcOffset {
        tail_fraction: f64,
    },
    Apodization {
        config: ApodizationConfig,
    },
    ZeroFill {
        multiplier: usize,
    },

    // --- Denoising ---
    WaveletDenoise {
        config: WaveletConfig,
    },
    CadzowDenoise {
        config: CadzowConfig,
    },
    SavitzkyGolay {
        half_window: usize,
        poly_order: usize,
    },

    // --- Transform & Frequency Domain ---
    Fft1D {
        half_first_point: bool,
    },
    Ifft1D,
    PhaseCorrect {
        phi0_deg: f64,
        phi1_deg: f64,
        pivot_frac: f64,
        auto_entropy: bool,
        allow_first_order: bool,
    },
    BaselineCorrect {
        method: BaselineMethod,
        lambda_exp: f64, // 10^lambda_exp for arPLS, or degree for Poly
        poly_degree: usize,
    },
    Fft2D {
        transform_f1: bool,
    },

    // --- Analysis & Outputs ---
    PeakPick {
        snr_threshold: f64,
        min_separation_pts: usize,
        use_magnitude: bool,
    },
    CompareMetrics {
        label: String,
    },
    SpectrumPlot1D {
        label: String,
        visible: bool,
        show_real: bool,
        show_imag: bool,
        show_magnitude: bool,
        show_baseline: bool,
    },
    HeatmapPlot2D {
        label: String,
        visible: bool,
        use_magnitude: bool,
        floor_pct: f64,
    },
}

impl NmrNode {
    pub fn title(&self) -> &'static str {
        match self {
            NmrNode::Synthetic1D { .. } => "Synthetic 1D FID",
            NmrNode::Synthetic2DCosy { .. } => "Synthetic 2D COSY",
            NmrNode::LoadedDataset1D { .. } => "Loaded 1D Dataset",
            NmrNode::LoadedDataset2D { .. } => "Loaded 2D Dataset",
            NmrNode::RemoveDigitalFilter { .. } => "Digital Filter Removal",
            NmrNode::DcOffset { .. } => "DC Offset Removal",
            NmrNode::Apodization { .. } => "Apodization Window",
            NmrNode::ZeroFill { .. } => "Zero Fill",
            NmrNode::WaveletDenoise { .. } => "Wavelet Denoise (SWT)",
            NmrNode::CadzowDenoise { .. } => "Cadzow Hankel SVD",
            NmrNode::SavitzkyGolay { .. } => "Savitzky-Golay Smooth",
            NmrNode::Fft1D { .. } => "1D FFT",
            NmrNode::Ifft1D => "1D Inverse FFT",
            NmrNode::PhaseCorrect { .. } => "Phase Correction",
            NmrNode::BaselineCorrect { .. } => "Baseline Correction",
            NmrNode::Fft2D { .. } => "2D FFT",
            NmrNode::PeakPick { .. } => "Peak Picking",
            NmrNode::CompareMetrics { .. } => "Compare SNR / SSIM",
            NmrNode::SpectrumPlot1D { .. } => "1D Spectrum Plot",
            NmrNode::HeatmapPlot2D { .. } => "2D Heatmap Plot",
        }
    }

    pub fn input_pins(&self) -> &'static [(&'static str, PinKind)] {
        match self {
            NmrNode::Synthetic1D { .. }
            | NmrNode::Synthetic2DCosy { .. }
            | NmrNode::LoadedDataset1D { .. }
            | NmrNode::LoadedDataset2D { .. } => &[],
            NmrNode::RemoveDigitalFilter { .. }
            | NmrNode::DcOffset { .. }
            | NmrNode::Apodization { .. }
            | NmrNode::ZeroFill { .. }
            | NmrNode::WaveletDenoise { .. }
            | NmrNode::CadzowDenoise { .. }
            | NmrNode::SavitzkyGolay { .. }
            | NmrNode::Fft1D { .. }
            | NmrNode::Ifft1D
            | NmrNode::PhaseCorrect { .. }
            | NmrNode::BaselineCorrect { .. }
            | NmrNode::PeakPick { .. } => &[("Signal", PinKind::Signal1D)],
            NmrNode::Fft2D { .. } | NmrNode::HeatmapPlot2D { .. } => &[("2D Matrix", PinKind::Signal2D)],
            NmrNode::CompareMetrics { .. } => {
                &[("Test", PinKind::Signal1D), ("Reference", PinKind::Signal1D)]
            }
            NmrNode::SpectrumPlot1D { .. } => {
                &[("Signal", PinKind::Signal1D), ("Peaks (opt)", PinKind::PeakList)]
            }
        }
    }

    pub fn output_pins(&self) -> &'static [(&'static str, PinKind)] {
        match self {
            NmrNode::Synthetic1D { .. }
            | NmrNode::LoadedDataset1D { .. }
            | NmrNode::RemoveDigitalFilter { .. }
            | NmrNode::DcOffset { .. }
            | NmrNode::Apodization { .. }
            | NmrNode::ZeroFill { .. }
            | NmrNode::WaveletDenoise { .. }
            | NmrNode::CadzowDenoise { .. }
            | NmrNode::SavitzkyGolay { .. }
            | NmrNode::Fft1D { .. }
            | NmrNode::Ifft1D
            | NmrNode::PhaseCorrect { .. }
            | NmrNode::BaselineCorrect { .. } => &[("Signal", PinKind::Signal1D)],
            NmrNode::Synthetic2DCosy { .. }
            | NmrNode::LoadedDataset2D { .. }
            | NmrNode::Fft2D { .. } => &[("2D Matrix", PinKind::Signal2D)],
            NmrNode::PeakPick { .. } => &[("Peaks", PinKind::PeakList)],
            NmrNode::CompareMetrics { .. }
            | NmrNode::SpectrumPlot1D { .. }
            | NmrNode::HeatmapPlot2D { .. } => &[],
        }
    }

    pub fn default_synthetic_ethanol() -> Self {
        Self::Synthetic1D {
            preset: SyntheticPreset::Ethanol,
            points: 1024,
            sf_mhz: 500.0,
            sw_ppm: 10.0,
            carrier_ppm: 5.0,
            noise_sigma: 0.25,
            seed: 42,
        }
    }

    pub fn default_wavelet() -> Self {
        Self::WaveletDenoise {
            config: WaveletConfig {
                family: WaveletFamily::Sym4,
                levels: 4,
                mode: ThresholdMode::Soft,
                threshold_scale: 1.0,
            },
        }
    }

    pub fn default_cadzow() -> Self {
        Self::CadzowDenoise {
            config: CadzowConfig {
                rank: 16,
                iterations: 2,
                pencil_fraction: 0.4,
            },
        }
    }

    pub fn default_apodization() -> Self {
        Self::Apodization {
            config: ApodizationConfig {
                window: WindowType::Exponential,
                lb_hz: 0.8,
                ..Default::default()
            },
        }
    }
}

/// Runtime evaluated output value stored per graph node.
#[derive(Debug, Clone)]
pub enum NodeOutput {
    D1(Signal1D),
    D2(Signal2D),
    Peaks(PeakList),
    Metrics(QualityMetrics),
    Plot1D {
        label: String,
        signal: Signal1D,
        peaks: Option<PeakList>,
        visible: bool,
        show_real: bool,
        show_imag: bool,
        show_magnitude: bool,
        show_baseline: bool,
    },
    Plot2D {
        label: String,
        matrix: Signal2D,
        visible: bool,
        use_magnitude: bool,
        floor_pct: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkflowPreset {
    EthanolDenoiseCompare,
    PhaseAndBaseline,
    Cosy2D,
}

/// Build ready-to-run workflow presets for immediate exploration.
pub fn build_preset_snarl(preset: WorkflowPreset) -> egui_snarl::Snarl<NmrNode> {
    use egui::pos2;
    use egui_snarl::{InPinId, OutPinId, Snarl};

    let mut snarl = Snarl::new();
    match preset {
        WorkflowPreset::EthanolDenoiseCompare => {
            let clean_src = snarl.insert_node(
                pos2(40.0, 40.0),
                NmrNode::Synthetic1D {
                    preset: SyntheticPreset::Ethanol,
                    points: 1024,
                    sf_mhz: 500.0,
                    sw_ppm: 10.0,
                    carrier_ppm: 5.0,
                    noise_sigma: 0.0,
                    seed: 42,
                },
            );
            let noisy_src = snarl.insert_node(pos2(40.0, 260.0), NmrNode::default_synthetic_ethanol());

            let clean_fft = snarl.insert_node(pos2(300.0, 40.0), NmrNode::Fft1D { half_first_point: true });
            let noisy_fft = snarl.insert_node(pos2(300.0, 210.0), NmrNode::Fft1D { half_first_point: true });
            let swt = snarl.insert_node(pos2(540.0, 340.0), NmrNode::default_wavelet());

            let cadzow = snarl.insert_node(pos2(280.0, 490.0), NmrNode::default_cadzow());
            let cadzow_fft = snarl.insert_node(pos2(540.0, 490.0), NmrNode::Fft1D { half_first_point: true });

            let peaks = snarl.insert_node(
                pos2(780.0, 390.0),
                NmrNode::PeakPick {
                    snr_threshold: 6.0,
                    min_separation_pts: 4,
                    use_magnitude: false,
                },
            );
            let cmp = snarl.insert_node(
                pos2(780.0, 120.0),
                NmrNode::CompareMetrics {
                    label: "Cadzow vs Clean".into(),
                },
            );

            let plot_noisy = snarl.insert_node(
                pos2(780.0, 250.0),
                NmrNode::SpectrumPlot1D {
                    label: "Noisy Ethanol".into(),
                    visible: true,
                    show_real: true,
                    show_imag: false,
                    show_magnitude: false,
                    show_baseline: false,
                },
            );
            let plot_cadzow = snarl.insert_node(
                pos2(1020.0, 460.0),
                NmrNode::SpectrumPlot1D {
                    label: "Cadzow SVD".into(),
                    visible: true,
                    show_real: true,
                    show_imag: false,
                    show_magnitude: false,
                    show_baseline: false,
                },
            );
            let plot_swt = snarl.insert_node(
                pos2(780.0, 560.0),
                NmrNode::SpectrumPlot1D {
                    label: "Wavelet (Sym4)".into(),
                    visible: false,
                    show_real: true,
                    show_imag: false,
                    show_magnitude: false,
                    show_baseline: false,
                },
            );

            snarl.connect(OutPinId { node: clean_src, output: 0 }, InPinId { node: clean_fft, input: 0 });
            snarl.connect(OutPinId { node: noisy_src, output: 0 }, InPinId { node: noisy_fft, input: 0 });
            snarl.connect(OutPinId { node: noisy_fft, output: 0 }, InPinId { node: swt, input: 0 });
            snarl.connect(OutPinId { node: noisy_src, output: 0 }, InPinId { node: cadzow, input: 0 });
            snarl.connect(OutPinId { node: cadzow, output: 0 }, InPinId { node: cadzow_fft, input: 0 });

            snarl.connect(OutPinId { node: cadzow_fft, output: 0 }, InPinId { node: peaks, input: 0 });
            snarl.connect(OutPinId { node: cadzow_fft, output: 0 }, InPinId { node: cmp, input: 0 });
            snarl.connect(OutPinId { node: clean_fft, output: 0 }, InPinId { node: cmp, input: 1 });

            snarl.connect(OutPinId { node: noisy_fft, output: 0 }, InPinId { node: plot_noisy, input: 0 });
            snarl.connect(OutPinId { node: cadzow_fft, output: 0 }, InPinId { node: plot_cadzow, input: 0 });
            snarl.connect(OutPinId { node: peaks, output: 0 }, InPinId { node: plot_cadzow, input: 1 });
            snarl.connect(OutPinId { node: swt, output: 0 }, InPinId { node: plot_swt, input: 0 });
        }
        WorkflowPreset::PhaseAndBaseline => {
            let src = snarl.insert_node(
                pos2(40.0, 120.0),
                NmrNode::Synthetic1D {
                    preset: SyntheticPreset::Ethylbenzene,
                    points: 2048,
                    sf_mhz: 500.0,
                    sw_ppm: 10.0,
                    carrier_ppm: 5.0,
                    noise_sigma: 0.08,
                    seed: 7,
                },
            );
            let apod = snarl.insert_node(pos2(280.0, 120.0), NmrNode::default_apodization());
            let zf = snarl.insert_node(pos2(510.0, 120.0), NmrNode::ZeroFill { multiplier: 2 });
            let fft = snarl.insert_node(pos2(710.0, 120.0), NmrNode::Fft1D { half_first_point: true });
            let phase = snarl.insert_node(
                pos2(910.0, 120.0),
                NmrNode::PhaseCorrect {
                    phi0_deg: 0.0,
                    phi1_deg: 0.0,
                    pivot_frac: 0.5,
                    auto_entropy: true,
                    allow_first_order: false,
                },
            );
            let bl = snarl.insert_node(
                pos2(1140.0, 120.0),
                NmrNode::BaselineCorrect {
                    method: BaselineMethod::ArPls,
                    lambda_exp: 6.0,
                    poly_degree: 3,
                },
            );
            let peaks = snarl.insert_node(
                pos2(1140.0, 320.0),
                NmrNode::PeakPick {
                    snr_threshold: 8.0,
                    min_separation_pts: 4,
                    use_magnitude: false,
                },
            );
            let plot = snarl.insert_node(
                pos2(1380.0, 180.0),
                NmrNode::SpectrumPlot1D {
                    label: "Ethylbenzene (Phased + arPLS)".into(),
                    visible: true,
                    show_real: true,
                    show_imag: false,
                    show_magnitude: false,
                    show_baseline: true,
                },
            );
            snarl.connect(OutPinId { node: src, output: 0 }, InPinId { node: apod, input: 0 });
            snarl.connect(OutPinId { node: apod, output: 0 }, InPinId { node: zf, input: 0 });
            snarl.connect(OutPinId { node: zf, output: 0 }, InPinId { node: fft, input: 0 });
            snarl.connect(OutPinId { node: fft, output: 0 }, InPinId { node: phase, input: 0 });
            snarl.connect(OutPinId { node: phase, output: 0 }, InPinId { node: bl, input: 0 });
            snarl.connect(OutPinId { node: bl, output: 0 }, InPinId { node: peaks, input: 0 });
            snarl.connect(OutPinId { node: bl, output: 0 }, InPinId { node: plot, input: 0 });
            snarl.connect(OutPinId { node: peaks, output: 0 }, InPinId { node: plot, input: 1 });
        }
        WorkflowPreset::Cosy2D => {
            let src = snarl.insert_node(
                pos2(40.0, 140.0),
                NmrNode::Synthetic2DCosy {
                    rows: 128,
                    cols: 128,
                    noise_sigma: 0.05,
                    seed: 12,
                },
            );
            let fft2 = snarl.insert_node(pos2(300.0, 140.0), NmrNode::Fft2D { transform_f1: true });
            let heat = snarl.insert_node(
                pos2(560.0, 140.0),
                NmrNode::HeatmapPlot2D {
                    label: "2D COSY Magnitude".into(),
                    visible: true,
                    use_magnitude: true,
                    floor_pct: 5.0,
                },
            );
            snarl.connect(OutPinId { node: src, output: 0 }, InPinId { node: fft2, input: 0 });
            snarl.connect(OutPinId { node: fft2, output: 0 }, InPinId { node: heat, input: 0 });
        }
    }
    snarl
}
