use clap::{Parser, Subcommand, ValueEnum};
use nmr_core::{Domain, compute_metrics};
use nmr_dsp::{
    ApodizationConfig, BaselineMethod, CadzowConfig, SyntheticPreset, ThresholdMode, WaveletConfig,
    WaveletFamily, apply_apodization, auto_phase_entropy, cadzow_denoise, dc_offset_correction,
    fft_1d, generate_synthetic_fid, pick_peaks, remove_digital_filter, subtract_baseline,
    wavelet_denoise, zero_fill,
};
use nmr_io::export::{export_csv, export_peaks_csv};
use nmr_io::jcamp::write_jcamp_dx;
use nmr_io::pipe::write_nmrpipe_1d;
use nmr_io::{LoadedDataset, load_from_path};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "nmr-read",
    version,
    about = "NMR-READ CLI: Multi-vendor NMR format inspector, converter, and pure-Rust denoising pipeline"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Clone, Copy, ValueEnum)]
enum ExportFormat {
    Csv,
    Jdx,
    Pipe,
}

#[derive(Clone, Copy, ValueEnum)]
enum PresetArg {
    Ethanol,
    Ethylbenzene,
    TwoPeaks,
}

#[derive(Subcommand)]
enum Commands {
    /// Inspect an NMR file or experiment directory (Bruker, Varian, JEOL, Spinsolve, NMRPipe, JCAMP-DX, Sparky, nmrML, ZIP)
    Inspect { path: PathBuf },
    /// Process and denoise an NMR dataset with optional FFT, auto-phase, baseline, wavelet, and Cadzow
    Process {
        input: PathBuf,
        #[arg(long)]
        remove_df: bool,
        #[arg(long)]
        lb_hz: Option<f64>,
        #[arg(long, default_value_t = false)]
        fft: bool,
        #[arg(long, default_value_t = false)]
        auto_phase: bool,
        #[arg(long)]
        baseline_lambda: Option<f64>,
        #[arg(long)]
        wavelet_scale: Option<f64>,
        #[arg(long)]
        cadzow_rank: Option<usize>,
        #[arg(long, default_value_t = 3)]
        cadzow_iters: usize,
        #[arg(long)]
        out_csv: Option<PathBuf>,
        #[arg(long)]
        out_peaks: Option<PathBuf>,
    },
    /// Convert any supported 1D NMR dataset to CSV, JCAMP-DX (.jdx), or NMRPipe (.ft1/.fid)
    Convert {
        input: PathBuf,
        #[arg(long, value_enum)]
        to: ExportFormat,
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Generate a synthetic 1H NMR spectrum and benchmark Wavelet vs Cadzow denoising
    Benchmark {
        #[arg(long, value_enum, default_value_t = PresetArg::Ethanol)]
        preset: PresetArg,
        #[arg(long, default_value_t = 2048)]
        points: usize,
        #[arg(long, default_value_t = 0.25)]
        noise: f64,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Inspect { path } => {
            let ds = load_from_path(&path)?;
            match ds {
                LoadedDataset::D1(s) => {
                    println!("Dataset: {} (1D)", s.meta.title);
                    println!("Format:  {:?}", s.meta.format);
                    println!("Domain:  {:?}", s.domain);
                    println!("Points:  {}", s.len());
                    println!(
                        "Nucleus: {} | SF: {:.4} MHz | SW: {:.2} Hz ({:.3} ppm) | Carrier: {:.3} ppm",
                        s.meta.axis.nucleus,
                        s.meta.axis.spectrometer_freq_mhz,
                        s.meta.axis.spectral_width_hz,
                        s.meta.axis.sw_ppm(),
                        s.meta.axis.carrier_ppm
                    );
                    if let Some(d) = s.meta.group_delay {
                        println!("Group Delay (GRPDLY): {d:.4} pts");
                    }
                    let m = compute_metrics(&s.real(), None);
                    println!("SNR: {:.2} dB (sigma={:.4e})", m.snr_db, m.noise_sigma);
                }
                LoadedDataset::D2(s2) => {
                    println!("Dataset: {} (2D)", s2.meta.title);
                    println!("Format:  {:?}", s2.meta.format);
                    println!("Shape:   {} rows (F1) x {} cols (F2)", s2.rows, s2.cols);
                    println!("Domains: F2={:?}, F1={:?}", s2.f2_domain, s2.f1_domain);
                }
            }
        }
        Commands::Process {
            input,
            remove_df,
            lb_hz,
            fft,
            auto_phase,
            baseline_lambda,
            wavelet_scale,
            cadzow_rank,
            cadzow_iters,
            out_csv,
            out_peaks,
        } => {
            let mut sig = load_from_path(&input)?.into_1d()?;
            if sig.domain == Domain::Time {
                if remove_df {
                    sig = remove_digital_filter(&sig, None);
                }
                sig = dc_offset_correction(&sig, 0.15);
                if let Some(rank) = cadzow_rank {
                    sig = cadzow_denoise(
                        &sig,
                        &CadzowConfig {
                            rank,
                            iterations: cadzow_iters,
                            pencil_fraction: 0.4,
                        },
                    );
                }
                if let Some(scale) = wavelet_scale {
                    sig = wavelet_denoise(
                        &sig,
                        &WaveletConfig {
                            family: WaveletFamily::Sym4,
                            levels: 4,
                            mode: ThresholdMode::Soft,
                            threshold_scale: scale,
                        },
                    );
                }
                if let Some(lb) = lb_hz {
                    let (apod, _) = apply_apodization(
                        &sig,
                        &ApodizationConfig {
                            lb_hz: lb,
                            ..Default::default()
                        },
                    );
                    sig = apod;
                }
                if fft {
                    sig = zero_fill(&sig, 0, 2);
                    sig = fft_1d(&sig, true);
                }
            } else if let Some(scale) = wavelet_scale {
                sig = wavelet_denoise(
                    &sig,
                    &WaveletConfig {
                        family: WaveletFamily::Sym4,
                        levels: 4,
                        mode: ThresholdMode::Soft,
                        threshold_scale: scale,
                    },
                );
            }

            if sig.domain == Domain::Frequency && auto_phase {
                let (phased, p0, p1) = auto_phase_entropy(&sig, true);
                sig = phased;
                eprintln!("Auto-phase applied: phi0={p0:.2} deg, phi1={p1:.2} deg");
            }
            if sig.domain == Domain::Frequency {
                if let Some(lam) = baseline_lambda {
                    sig = subtract_baseline(&sig, BaselineMethod::ArPls, lam);
                }
            }

            let peaks = pick_peaks(&sig, 5.0, 4, false);
            let m = compute_metrics(&sig.real(), None);
            println!(
                "Processed {} points ({:?}) | SNR={:.2} dB | Peaks={}",
                sig.len(),
                sig.domain,
                m.snr_db,
                peaks.peaks.len()
            );

            if let Some(csv_path) = out_csv {
                std::fs::write(&csv_path, export_csv(&sig))?;
                println!("Wrote spectrum CSV -> {}", csv_path.display());
            }
            if let Some(pk_path) = out_peaks {
                std::fs::write(&pk_path, export_peaks_csv(&peaks))?;
                println!("Wrote peak list CSV -> {}", pk_path.display());
            }
        }
        Commands::Convert { input, to, output } => {
            let sig = load_from_path(&input)?.into_1d()?;
            match to {
                ExportFormat::Csv => std::fs::write(&output, export_csv(&sig))?,
                ExportFormat::Jdx => std::fs::write(&output, write_jcamp_dx(&sig))?,
                ExportFormat::Pipe => std::fs::write(&output, write_nmrpipe_1d(&sig))?,
            }
            println!("Converted {} -> {}", input.display(), output.display());
        }
        Commands::Benchmark {
            preset,
            points,
            noise,
        } => {
            let p = match preset {
                PresetArg::Ethanol => SyntheticPreset::Ethanol,
                PresetArg::Ethylbenzene => SyntheticPreset::Ethylbenzene,
                PresetArg::TwoPeaks => SyntheticPreset::TwoPeaks,
            };
            let clean_fid = generate_synthetic_fid(p, points, 500.0, 10.0, 5.0, 0.0, 42);
            let noisy_fid = generate_synthetic_fid(p, points, 500.0, 10.0, 5.0, noise, 42);

            let clean_spec_sig = fft_1d(&clean_fid, true);
            let clean_spec = clean_spec_sig.real();
            let noisy_spec_sig = fft_1d(&noisy_fid, true);
            let noisy_spec = noisy_spec_sig.real();
            let swt_spec = wavelet_denoise(&noisy_spec_sig, &WaveletConfig::default()).real();
            let cad_fid = cadzow_denoise(
                &noisy_fid,
                &CadzowConfig {
                    rank: 16,
                    iterations: 2,
                    pencil_fraction: 0.4,
                },
            );
            let cad_spec = fft_1d(&cad_fid, true).real();

            let m_noisy = compute_metrics(&noisy_spec, Some(&clean_spec));
            let m_swt = compute_metrics(&swt_spec, Some(&clean_spec));
            let m_cad = compute_metrics(&cad_spec, Some(&clean_spec));

            println!("=== NMR-READ Pure-Rust Denoising Benchmark ({p:?}, N={points}, noise={noise}) ===");
            println!(
                "Noisy:   SNR={:.2} dB | RMSE={:.4} | SSIM={:.4}",
                m_noisy.snr_db,
                m_noisy.rmse.unwrap_or(0.0),
                m_noisy.ssim.unwrap_or(0.0)
            );
            println!(
                "Wavelet: SNR={:.2} dB | RMSE={:.4} | SSIM={:.4}",
                m_swt.snr_db,
                m_swt.rmse.unwrap_or(0.0),
                m_swt.ssim.unwrap_or(0.0)
            );
            println!(
                "Cadzow:  SNR={:.2} dB | RMSE={:.4} | SSIM={:.4}",
                m_cad.snr_db,
                m_cad.rmse.unwrap_or(0.0),
                m_cad.ssim.unwrap_or(0.0)
            );
        }
    }
    Ok(())
}
