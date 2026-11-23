use crate::engine::EvalResults;
use crate::nodes::NodeOutput;
use egui::{Color32, ColorImage, TextureHandle, TextureOptions, Ui};
use egui_plot::{Legend, Line, Plot, PlotImage, PlotPoint, PlotPoints, Points};
use nmr_core::{Domain, PeakList, Signal2D};
use std::collections::HashMap;

const PALETTE: &[Color32] = &[
    Color32::from_rgb(0x4f, 0xc3, 0xf7),
    Color32::from_rgb(0xff, 0x8a, 0x65),
    Color32::from_rgb(0x81, 0xc7, 0x84),
    Color32::from_rgb(0xf0, 0x62, 0x92),
    Color32::from_rgb(0xff, 0xd5, 0x4f),
];

/// State cache for 2D heatmap GPU textures.
#[derive(Default)]
pub struct PlotViewState {
    pub textures_2d: HashMap<String, (usize, usize, TextureHandle)>,
    pub active_tab: PlotTab,
    pub selected_peak_list: Option<PeakList>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PlotTab {
    #[default]
    Spectrum1D,
    Matrix2D,
    PeaksAndMetrics,
}

/// Render the interactive 1D / 2D NMR visualization dock and peak/metrics inspection tables.
pub fn show_plot_dock(ui: &mut Ui, eval: &EvalResults, state: &mut PlotViewState) {
    let has_1d = eval
        .outputs
        .values()
        .any(|v| matches!(v, NodeOutput::Plot1D { visible: true, .. }));
    let has_2d = eval
        .outputs
        .values()
        .any(|v| matches!(v, NodeOutput::Plot2D { visible: true, .. }));

    if has_2d && !has_1d && state.active_tab == PlotTab::Spectrum1D {
        state.active_tab = PlotTab::Matrix2D;
    } else if has_1d && !has_2d && state.active_tab == PlotTab::Matrix2D {
        state.active_tab = PlotTab::Spectrum1D;
    }

    ui.horizontal(|ui| {
        ui.selectable_value(&mut state.active_tab, PlotTab::Spectrum1D, "📈 1D Spectrum / FID");
        ui.selectable_value(&mut state.active_tab, PlotTab::Matrix2D, "🟦 2D Heatmap");
        ui.selectable_value(
            &mut state.active_tab,
            PlotTab::PeaksAndMetrics,
            "📊 Peaks & SNR Metrics",
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.weak(format!("Eval: {:.2} ms", eval.total_ms));
        });
    });
    ui.separator();

    match state.active_tab {
        PlotTab::Spectrum1D => show_1d_plot(ui, eval),
        PlotTab::Matrix2D => show_2d_plot(ui, eval, state),
        PlotTab::PeaksAndMetrics => show_peaks_and_metrics(ui, eval),
    }
}

fn show_1d_plot(ui: &mut Ui, eval: &EvalResults) {
    let mut plots_1d: Vec<_> = eval
        .outputs
        .iter()
        .filter_map(|(id, v)| match v {
            NodeOutput::Plot1D {
                label,
                signal,
                peaks,
                visible: true,
                show_real,
                show_imag,
                show_magnitude,
                show_baseline,
            } => Some((
                id.0,
                label,
                signal,
                peaks,
                *show_real,
                *show_imag,
                *show_magnitude,
                *show_baseline,
            )),
            _ => None,
        })
        .collect();
    plots_1d.sort_by_key(|(id, ..)| *id);

    if plots_1d.is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("Connect a signal to a '1D Spectrum Plot' node (and enable Visible) to view spectra.");
        });
        return;
    }

    let is_freq = plots_1d
        .iter()
        .any(|(_, _, s, ..)| s.domain == Domain::Frequency);

    let plot = Plot::new("nmr_1d_plot")
        .legend(Legend::default())
        .x_axis_label(if is_freq { "Chemical Shift (ppm)" } else { "Time (s)" })
        .y_axis_label("Intensity")
        .x_axis_formatter(move |mark, _range| {
            if is_freq {
                format!("{:.2}", -mark.value)
            } else {
                format!("{:.3}", mark.value)
            }
        });

    plot.show(ui, |plot_ui| {
        for (idx, (_id, label, sig, peaks, show_re, show_im, show_mag, show_bl)) in
            plots_1d.iter().enumerate()
        {
            let base_color = PALETTE[idx % PALETTE.len()];
            let x_vals = if sig.domain == Domain::Frequency {
                sig.x_axis_ppm().into_iter().map(|p| -p).collect::<Vec<f64>>()
            } else {
                sig.x_axis_time_s()
            };

            if *show_re {
                let pts: PlotPoints = x_vals
                    .iter()
                    .zip(sig.data.iter())
                    .map(|(&x, c)| [x, c.re])
                    .collect();
                plot_ui.line(Line::new(format!("{label} (Re)"), pts).color(base_color).width(1.6));
            }
            if *show_im {
                let pts: PlotPoints = x_vals
                    .iter()
                    .zip(sig.data.iter())
                    .map(|(&x, c)| [x, c.im])
                    .collect();
                plot_ui.line(
                    Line::new(format!("{label} (Im)"), pts)
                        .color(base_color.gamma_multiply(0.6))
                        .width(1.0),
                );
            }
            if *show_mag {
                let pts: PlotPoints = x_vals
                    .iter()
                    .zip(sig.data.iter())
                    .map(|(&x, c)| [x, c.norm()])
                    .collect();
                plot_ui.line(
                    Line::new(format!("{label} (|Mag|)"), pts)
                        .color(Color32::from_rgb(0xba, 0x68, 0xc8))
                        .width(1.4),
                );
            }
            if *show_bl {
                if let Some(bl) = &sig.baseline {
                    if bl.len() == x_vals.len() {
                        let pts: PlotPoints = x_vals
                            .iter()
                            .zip(bl.iter())
                            .map(|(&x, &b)| [x, b])
                            .collect();
                        plot_ui.line(
                            Line::new(format!("{label} (Baseline)"), pts)
                                .color(Color32::from_rgb(0xef, 0x53, 0x50))
                                .width(1.2),
                        );
                    }
                }
            }
            if let Some(pk_list) = peaks {
                let pk_pts: PlotPoints = pk_list
                    .peaks
                    .iter()
                    .map(|p| {
                        let x = if sig.domain == Domain::Frequency {
                            -p.ppm
                        } else {
                            p.index as f64 * sig.meta.axis.dwell_time_s()
                        };
                        [x, p.intensity]
                    })
                    .collect();
                plot_ui.points(
                    Points::new(format!("{label} Peaks"), pk_pts)
                        .color(Color32::from_rgb(0xff, 0xd5, 0x4f))
                        .radius(4.0),
                );
            }
        }
    });
}

fn show_2d_plot(ui: &mut Ui, eval: &EvalResults, state: &mut PlotViewState) {
    let plot2d = eval.outputs.iter().find_map(|(id, v)| match v {
        NodeOutput::Plot2D {
            label,
            matrix,
            visible: true,
            use_magnitude,
            floor_pct,
        } => Some((id.0, label, matrix, *use_magnitude, *floor_pct)),
        _ => None,
    });

    let Some((id, label, matrix, use_mag, floor_pct)) = plot2d else {
        ui.centered_and_justified(|ui| {
            ui.label("Connect a 2D dataset or Synthetic 2D COSY -> 2D FFT -> 2D Heatmap Plot.");
        });
        return;
    };

    let tex_key = format!("{id}_{label}_{use_mag}_{floor_pct:.1}");
    let tex_handle = if let Some((r, c, handle)) = state.textures_2d.get(&tex_key) {
        if *r == matrix.rows && *c == matrix.cols {
            handle.clone()
        } else {
            let img = build_2d_color_image(matrix, use_mag, floor_pct);
            let h = ui
                .ctx()
                .load_texture(tex_key.clone(), img, TextureOptions::LINEAR);
            state
                .textures_2d
                .insert(tex_key, (matrix.rows, matrix.cols, h.clone()));
            h
        }
    } else {
        let img = build_2d_color_image(matrix, use_mag, floor_pct);
        let h = ui
            .ctx()
            .load_texture(tex_key.clone(), img, TextureOptions::LINEAR);
        state
            .textures_2d
            .insert(tex_key, (matrix.rows, matrix.cols, h.clone()));
        h
    };

    let f2_left = matrix.meta.axis.left_ppm();
    let f2_right = matrix.meta.axis.right_ppm();
    let f1_axis = matrix.meta.indirect_axis.as_ref().unwrap_or(&matrix.meta.axis);
    let f1_left = f1_axis.left_ppm();
    let f1_right = f1_axis.right_ppm();

    let center = PlotPoint::new(
        -0.5 * (f2_left + f2_right),
        -0.5 * (f1_left + f1_right),
    );
    let size = egui::vec2(
        (f2_left - f2_right).abs().max(0.1) as f32,
        (f1_left - f1_right).abs().max(0.1) as f32,
    );

    Plot::new("nmr_2d_heatmap")
        .x_axis_label("F2 (ppm)")
        .y_axis_label("F1 (ppm)")
        .x_axis_formatter(|mark, _| format!("{:.2}", -mark.value))
        .y_axis_formatter(|mark, _| format!("{:.2}", -mark.value))
        .show(ui, |plot_ui| {
            plot_ui.image(PlotImage::new(label, tex_handle.id(), center, size));
        });
}

fn build_2d_color_image(matrix: &Signal2D, use_magnitude: bool, floor_pct: f64) -> ColorImage {
    let vals = if use_magnitude {
        matrix.magnitude_matrix()
    } else {
        matrix.real_matrix()
    };
    let max_val = vals.iter().fold(0.0_f64, |m, &v| m.max(v.abs())).max(1e-12);
    let floor = (floor_pct.clamp(0.0, 80.0) / 100.0) * max_val;

    let mut pixels = Vec::with_capacity(matrix.rows * matrix.cols);
    for v in vals {
        let norm = ((v.abs() - floor) / (max_val - floor).max(1e-12)).clamp(0.0, 1.0);
        // Inferno-inspired perceptual ramp on dark background
        let r = (255.0 * (norm.powf(0.55))) as u8;
        let g = (255.0 * (norm.powf(1.35) * 0.9)) as u8;
        let b = (255.0 * ((1.0 - (norm - 0.35).abs() * 1.6).max(0.06))) as u8;
        let alpha = if norm <= 0.005 { 30 } else { 255 };
        pixels.push(Color32::from_rgba_unmultiplied(r, g, b, alpha));
    }
    ColorImage {
        size: [matrix.cols, matrix.rows],
        source_size: egui::vec2(matrix.cols as f32, matrix.rows as f32),
        pixels,
    }
}

fn show_peaks_and_metrics(ui: &mut Ui, eval: &EvalResults) {
    egui::ScrollArea::vertical().show(ui, |ui| {
        ui.heading("Denoising & Quality Comparison");
        let mut any_metric = false;
        for (id, out) in &eval.outputs {
            if let NodeOutput::Metrics(m) = out {
                any_metric = true;
                ui.group(|ui| {
                    ui.strong(format!("Node #{} Metrics", id.0));
                    ui.label(format!(
                        "SNR: {:.2} dB (linear {:.1}) | Noise sigma: {:.4e}",
                        m.snr_db, m.snr_linear, m.noise_sigma
                    ));
                    if let Some(rmse) = m.rmse {
                        ui.label(format!(
                            "RMSE: {:.5} | PSNR: {:.2} dB | SSIM: {:.4}",
                            rmse,
                            m.psnr_db.unwrap_or(0.0),
                            m.ssim.unwrap_or(0.0)
                        ));
                    }
                });
            }
        }
        if !any_metric {
            ui.weak("Add or connect a 'Compare SNR / SSIM' node to evaluate denoising metrics.");
        }

        ui.add_space(10.0);
        ui.heading("Picked Resonance Peaks");
        for (id, out) in &eval.outputs {
            if let NodeOutput::Peaks(pk_list) = out {
                ui.strong(format!(
                    "Node #{} — {} peaks (noise sigma = {:.3e})",
                    id.0,
                    pk_list.peaks.len(),
                    pk_list.noise_sigma
                ));
                egui::Grid::new(format!("peaks_grid_{}", id.0))
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("#");
                        ui.strong("Shift (ppm)");
                        ui.strong("Freq (Hz)");
                        ui.strong("Intensity");
                        ui.strong("FWHM (Hz)");
                        ui.strong("SNR");
                        ui.end_row();
                        for (idx, p) in pk_list.peaks.iter().enumerate() {
                            ui.label(format!("{}", idx + 1));
                            ui.label(format!("{:.4}", p.ppm));
                            ui.label(format!("{:.2}", p.hz));
                            ui.label(format!("{:.4}", p.intensity));
                            ui.label(
                                p.fwhm_hz
                                    .map(|w| format!("{w:.2}"))
                                    .unwrap_or_else(|| "-".into()),
                            );
                            ui.label(format!("{:.1}", p.snr));
                            ui.end_row();
                        }
                    });
                ui.add_space(8.0);
            }
        }
    });
}
