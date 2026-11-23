use crate::nodes::{NmrNode, NodeOutput};
use egui_snarl::{NodeId, Snarl};
use nmr_core::{Signal1D, Signal2D, compute_metrics};
use nmr_dsp::{
    apply_apodization, auto_phase_entropy, cadzow_denoise, dc_offset_correction, fft_1d, fft_2d,
    generate_synthetic_2d_cosy, generate_synthetic_fid, ifft_1d, phase_correct, pick_peaks,
    remove_digital_filter, savitzky_golay_smooth, subtract_baseline, wavelet_denoise, zero_fill,
};
use nmr_io::LoadedDataset;
use std::collections::{BTreeMap, HashMap, HashSet};

/// Cross-platform elapsed timer compatible with both native and `wasm32-unknown-unknown`.
#[cfg(not(target_arch = "wasm32"))]
fn now_ms() -> f64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs_f64() * 1000.0)
        .unwrap_or(0.0)
}

#[cfg(target_arch = "wasm32")]
fn now_ms() -> f64 {
    web_sys::window()
        .and_then(|w| w.performance())
        .map(|p| p.now())
        .unwrap_or(0.0)
}

/// Evaluation results for the entire node graph.
#[derive(Default)]
pub struct EvalResults {
    pub outputs: HashMap<NodeId, NodeOutput>,
    pub timings_ms: HashMap<NodeId, f64>,
    pub errors: HashMap<NodeId, String>,
    pub total_ms: f64,
}

/// Evaluate all reachable nodes in topological order without recursion overflow or cycle panics.
pub fn evaluate_graph(
    snarl: &Snarl<NmrNode>,
    datasets_1d: &BTreeMap<String, Signal1D>,
    datasets_2d: &BTreeMap<String, Signal2D>,
) -> EvalResults {
    let t_start = now_ms();
    let mut results = EvalResults::default();

    // Build incoming wire lookup: (target_node, input_idx) -> source_node
    let mut incoming: HashMap<(NodeId, usize), NodeId> = HashMap::new();
    for (out_pin, in_pin) in snarl.wires() {
        incoming.insert((in_pin.node, in_pin.input), out_pin.node);
    }

    let all_ids: Vec<NodeId> = snarl.node_ids().map(|(id, _)| id).collect();
    let mut visited = HashSet::new();
    let mut visiting = HashSet::new();
    let mut order = Vec::with_capacity(all_ids.len());

    fn dfs(
        node: NodeId,
        snarl: &Snarl<NmrNode>,
        incoming: &HashMap<(NodeId, usize), NodeId>,
        visited: &mut HashSet<NodeId>,
        visiting: &mut HashSet<NodeId>,
        order: &mut Vec<NodeId>,
        errors: &mut HashMap<NodeId, String>,
    ) {
        if visited.contains(&node) {
            return;
        }
        if !visiting.insert(node) {
            errors.insert(node, "Cycle detected in node graph".into());
            return;
        }
        if let Some(n_data) = snarl.get_node(node) {
            for idx in 0..n_data.input_pins().len() {
                if let Some(&src) = incoming.get(&(node, idx)) {
                    dfs(src, snarl, incoming, visited, visiting, order, errors);
                }
            }
        }
        visiting.remove(&node);
        visited.insert(node);
        order.push(node);
    }

    for id in all_ids {
        dfs(
            id,
            snarl,
            &incoming,
            &mut visited,
            &mut visiting,
            &mut order,
            &mut results.errors,
        );
    }

    for node_id in order {
        let Some(node) = snarl.get_node(node_id) else {
            continue;
        };
        let t0 = now_ms();

        let get_in_d1 = |pin_idx: usize, outs: &HashMap<NodeId, NodeOutput>| -> Option<Signal1D> {
            let src = incoming.get(&(node_id, pin_idx))?;
            match outs.get(src)? {
                NodeOutput::D1(s) => Some(s.clone()),
                _ => None,
            }
        };
        let get_in_d2 = |pin_idx: usize, outs: &HashMap<NodeId, NodeOutput>| -> Option<Signal2D> {
            let src = incoming.get(&(node_id, pin_idx))?;
            match outs.get(src)? {
                NodeOutput::D2(s) => Some(s.clone()),
                _ => None,
            }
        };
        let get_in_peaks = |pin_idx: usize, outs: &HashMap<NodeId, NodeOutput>| -> Option<nmr_core::PeakList> {
            let src = incoming.get(&(node_id, pin_idx))?;
            match outs.get(src)? {
                NodeOutput::Peaks(p) => Some(p.clone()),
                _ => None,
            }
        };

        let produced: Option<NodeOutput> = match node {
            NmrNode::Synthetic1D {
                preset,
                points,
                sf_mhz,
                sw_ppm,
                carrier_ppm,
                noise_sigma,
                seed,
            } => Some(NodeOutput::D1(generate_synthetic_fid(
                *preset,
                *points,
                *sf_mhz,
                *sw_ppm,
                *carrier_ppm,
                *noise_sigma,
                *seed,
            ))),
            NmrNode::Synthetic2DCosy {
                rows,
                cols,
                noise_sigma,
                seed,
            } => Some(NodeOutput::D2(generate_synthetic_2d_cosy(
                *rows,
                *cols,
                *noise_sigma,
                *seed,
            ))),
            NmrNode::LoadedDataset1D { dataset_key } => {
                if let Some(s) = datasets_1d.get(dataset_key) {
                    Some(NodeOutput::D1(s.clone()))
                } else {
                    results
                        .errors
                        .insert(node_id, format!("Dataset '{dataset_key}' not loaded"));
                    None
                }
            }
            NmrNode::LoadedDataset2D { dataset_key } => {
                if let Some(s) = datasets_2d.get(dataset_key) {
                    Some(NodeOutput::D2(s.clone()))
                } else {
                    results
                        .errors
                        .insert(node_id, format!("2D dataset '{dataset_key}' not loaded"));
                    None
                }
            }
            NmrNode::RemoveDigitalFilter {
                use_override,
                grpdly_override,
            } => get_in_d1(0, &results.outputs).map(|s| {
                NodeOutput::D1(remove_digital_filter(
                    &s,
                    if *use_override {
                        Some(*grpdly_override)
                    } else {
                        None
                    },
                ))
            }),
            NmrNode::DcOffset { tail_fraction } => {
                get_in_d1(0, &results.outputs).map(|s| NodeOutput::D1(dc_offset_correction(&s, *tail_fraction)))
            }
            NmrNode::Apodization { config } => get_in_d1(0, &results.outputs).map(|s| {
                let (out, _env) = apply_apodization(&s, config);
                NodeOutput::D1(out)
            }),
            NmrNode::ZeroFill { multiplier } => {
                get_in_d1(0, &results.outputs).map(|s| NodeOutput::D1(zero_fill(&s, 0, *multiplier)))
            }
            NmrNode::WaveletDenoise { config } => {
                get_in_d1(0, &results.outputs).map(|s| NodeOutput::D1(wavelet_denoise(&s, config)))
            }
            NmrNode::CadzowDenoise { config } => {
                get_in_d1(0, &results.outputs).map(|s| NodeOutput::D1(cadzow_denoise(&s, config)))
            }
            NmrNode::SavitzkyGolay {
                half_window,
                poly_order,
            } => get_in_d1(0, &results.outputs)
                .map(|s| NodeOutput::D1(savitzky_golay_smooth(&s, *half_window, *poly_order))),
            NmrNode::Fft1D { half_first_point } => {
                get_in_d1(0, &results.outputs).map(|s| NodeOutput::D1(fft_1d(&s, *half_first_point)))
            }
            NmrNode::Ifft1D => get_in_d1(0, &results.outputs).map(|s| NodeOutput::D1(ifft_1d(&s))),
            NmrNode::PhaseCorrect {
                phi0_deg,
                phi1_deg,
                pivot_frac,
                auto_entropy,
                allow_first_order,
            } => get_in_d1(0, &results.outputs).map(|s| {
                let out = if *auto_entropy {
                    auto_phase_entropy(&s, *allow_first_order).0
                } else {
                    phase_correct(&s, *phi0_deg, *phi1_deg, *pivot_frac)
                };
                NodeOutput::D1(out)
            }),
            NmrNode::BaselineCorrect {
                method,
                lambda_exp,
                poly_degree,
            } => get_in_d1(0, &results.outputs).map(|s| {
                let param = match method {
                    nmr_dsp::BaselineMethod::ArPls => 10.0_f64.powf(*lambda_exp),
                    nmr_dsp::BaselineMethod::Polynomial => *poly_degree as f64,
                };
                NodeOutput::D1(subtract_baseline(&s, *method, param))
            }),
            NmrNode::Fft2D { transform_f1 } => {
                get_in_d2(0, &results.outputs).map(|s| NodeOutput::D2(fft_2d(&s, *transform_f1)))
            }
            NmrNode::PeakPick {
                snr_threshold,
                min_separation_pts,
                use_magnitude,
            } => get_in_d1(0, &results.outputs)
                .map(|s| NodeOutput::Peaks(pick_peaks(&s, *snr_threshold, *min_separation_pts, *use_magnitude))),
            NmrNode::CompareMetrics { .. } => {
                let test = get_in_d1(0, &results.outputs);
                let reference = get_in_d1(1, &results.outputs);
                test.map(|t| {
                    let ref_vec = reference.as_ref().map(|r| r.real());
                    NodeOutput::Metrics(compute_metrics(&t.real(), ref_vec.as_deref()))
                })
            }
            NmrNode::SpectrumPlot1D {
                label,
                visible,
                show_real,
                show_imag,
                show_magnitude,
                show_baseline,
            } => get_in_d1(0, &results.outputs).map(|signal| {
                let peaks = get_in_peaks(1, &results.outputs);
                NodeOutput::Plot1D {
                    label: label.clone(),
                    signal,
                    peaks,
                    visible: *visible,
                    show_real: *show_real,
                    show_imag: *show_imag,
                    show_magnitude: *show_magnitude,
                    show_baseline: *show_baseline,
                }
            }),
            NmrNode::HeatmapPlot2D {
                label,
                visible,
                use_magnitude,
                floor_pct,
            } => get_in_d2(0, &results.outputs).map(|matrix| NodeOutput::Plot2D {
                label: label.clone(),
                matrix,
                visible: *visible,
                use_magnitude: *use_magnitude,
                floor_pct: *floor_pct,
            }),
        };

        if let Some(val) = produced {
            results.outputs.insert(node_id, val);
        }
        results.timings_ms.insert(node_id, (now_ms() - t0).max(0.0));
    }

    results.total_ms = (now_ms() - t_start).max(0.0);
    results
}

pub fn register_loaded_dataset(
    name: String,
    ds: LoadedDataset,
    datasets_1d: &mut BTreeMap<String, Signal1D>,
    datasets_2d: &mut BTreeMap<String, Signal2D>,
) -> (bool, String) {
    match ds {
        LoadedDataset::D1(s) => {
            datasets_1d.insert(name.clone(), s);
            (false, name)
        }
        LoadedDataset::D2(s2) => {
            datasets_2d.insert(name.clone(), s2);
            (true, name)
        }
    }
}
