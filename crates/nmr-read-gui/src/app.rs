use crate::engine::{EvalResults, evaluate_graph, register_loaded_dataset};
use crate::nodes::{NmrNode, NodeOutput, WorkflowPreset, build_preset_snarl};
use crate::plot::{PlotViewState, show_plot_dock};
use egui::{Color32, Ui};
use egui_snarl::{
    InPin, NodeId, OutPin, Snarl,
    ui::{PinInfo, SnarlPin, SnarlStyle, SnarlViewer, SnarlWidget},
};
use nmr_core::{Signal1D, Signal2D};
use nmr_dsp::{BaselineMethod, SyntheticPreset, ThresholdMode, WaveletFamily, WindowType};
use nmr_io::export::{export_csv, export_peaks_csv};
use nmr_io::jcamp::write_jcamp_dx;
use nmr_io::pipe::write_nmrpipe_1d;
use nmr_io::{LoadedDataset, VirtualTree, load_from_bytes, load_from_virtual_tree};
use std::collections::BTreeMap;
use std::sync::mpsc::{Receiver, Sender};

enum AsyncLoadMsg {
    SingleFile { name: String, bytes: Vec<u8> },
    VirtualDir(VirtualTree),
    Error(String),
}

pub struct NmrReadApp {
    pub snarl: Snarl<NmrNode>,
    pub datasets_1d: BTreeMap<String, Signal1D>,
    pub datasets_2d: BTreeMap<String, Signal2D>,
    pub eval: EvalResults,
    pub plot_state: PlotViewState,
    pub dirty: bool,
    pub status_msg: String,
    rx: Receiver<AsyncLoadMsg>,
    tx: Sender<AsyncLoadMsg>,
}

impl Default for NmrReadApp {
    fn default() -> Self {
        let snarl = build_preset_snarl(WorkflowPreset::EthanolDenoiseCompare);
        let datasets_1d = BTreeMap::new();
        let datasets_2d = BTreeMap::new();
        let eval = evaluate_graph(&snarl, &datasets_1d, &datasets_2d);
        let (tx, rx) = std::sync::mpsc::channel();
        Self {
            snarl,
            datasets_1d,
            datasets_2d,
            eval,
            plot_state: PlotViewState::default(),
            dirty: false,
            status_msg: "Ready — Pure-Rust NMR Studio (Desktop & Wasm)".into(),
            rx,
            tx,
        }
    }
}

impl NmrReadApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        cc.egui_ctx.set_visuals(egui::Visuals::dark());
        Self::default()
    }

    fn poll_async_loads(&mut self) {
        while let Ok(msg) = self.rx.try_recv() {
            match msg {
                AsyncLoadMsg::SingleFile { name, bytes } => {
                    if name.to_ascii_lowercase().ends_with(".nmrflow.json") {
                        match serde_json::from_slice::<Snarl<NmrNode>>(&bytes) {
                            Ok(loaded_snarl) => {
                                self.snarl = loaded_snarl;
                                self.dirty = true;
                                self.status_msg = format!("Loaded workflow: {name}");
                            }
                            Err(e) => self.status_msg = format!("Invalid workflow ({name}): {e}"),
                        }
                        continue;
                    }
                    match load_from_bytes(&name, &bytes) {
                        Ok(ds) => self.attach_loaded_dataset(name, ds),
                        Err(e) => self.status_msg = format!("Load error ({name}): {e}"),
                    }
                }
                AsyncLoadMsg::VirtualDir(tree) => {
                    let name = tree.root_name.clone();
                    match load_from_virtual_tree(&tree) {
                        Ok(ds) => self.attach_loaded_dataset(name, ds),
                        Err(e) => self.status_msg = format!("Folder load error ({name}): {e}"),
                    }
                }
                AsyncLoadMsg::Error(e) => {
                    self.status_msg = format!("Error: {e}");
                }
            }
        }
    }

    fn attach_loaded_dataset(&mut self, name: String, ds: LoadedDataset) {
        let (is_2d, key) = register_loaded_dataset(
            name,
            ds,
            &mut self.datasets_1d,
            &mut self.datasets_2d,
        );
        let pos = egui::pos2(60.0, 60.0 + (self.datasets_1d.len() + self.datasets_2d.len()) as f32 * 40.0);
        if is_2d {
            self.snarl
                .insert_node(pos, NmrNode::LoadedDataset2D { dataset_key: key.clone() });
            self.status_msg = format!("Loaded 2D dataset: {key}");
        } else {
            self.snarl
                .insert_node(pos, NmrNode::LoadedDataset1D { dataset_key: key.clone() });
            self.status_msg = format!("Loaded 1D dataset: {key}");
        }
        self.dirty = true;
    }

    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if dropped.is_empty() {
            return;
        }

        #[cfg(not(target_arch = "wasm32"))]
        {
            if dropped.len() > 1 {
                let mut tree = VirtualTree::new("dropped_dataset");
                let mut has_bytes = false;
                for f in &dropped {
                    let name = f
                        .path()
                        .file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or("file")
                        .to_string();
                    if let Ok(bytes) = f.bytes() {
                        tree.insert(&name, bytes);
                        has_bytes = true;
                    }
                }
                if has_bytes {
                    if let Ok(ds) = load_from_virtual_tree(&tree) {
                        self.attach_loaded_dataset("dropped_dataset".into(), ds);
                        return;
                    }
                }
            }

            for f in dropped {
                let path = f.path();
                let name = path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("dataset")
                    .to_string();
                if path.exists() {
                    match nmr_io::load_from_path(path) {
                        Ok(ds) => self.attach_loaded_dataset(name, ds),
                        Err(e) => self.status_msg = format!("Failed to load {name}: {e}"),
                    }
                } else if let Ok(bytes) = f.bytes() {
                    match load_from_bytes(&name, &bytes) {
                        Ok(ds) => self.attach_loaded_dataset(name, ds),
                        Err(e) => self.status_msg = format!("Failed to load {name}: {e}"),
                    }
                }
            }
        }

        #[cfg(target_arch = "wasm32")]
        {
            let tx = self.tx.clone();
            let ctx = ctx.clone();
            spawn_task(async move {
                if dropped.len() > 1 {
                    let mut tree = VirtualTree::new("dropped_dataset");
                    for f in &dropped {
                        let name = f
                            .path()
                            .file_name()
                            .and_then(|s| s.to_str())
                            .unwrap_or("file")
                            .to_string();
                        if let Ok(bytes) = f.bytes_async().await {
                            tree.insert(&name, bytes);
                        }
                    }
                    let _ = tx.send(AsyncLoadMsg::VirtualDir(tree));
                } else if let Some(f) = dropped.into_iter().next() {
                    let name = f
                        .path()
                        .file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or("dataset")
                        .to_string();
                    match f.bytes_async().await {
                        Ok(bytes) => {
                            let _ = tx.send(AsyncLoadMsg::SingleFile { name, bytes });
                        }
                        Err(e) => {
                            let _ = tx.send(AsyncLoadMsg::Error(e));
                        }
                    }
                }
                ctx.request_repaint();
            });
        }
    }

    fn spawn_open_file_dialog(&self, ctx: &egui::Context) {
        let tx = self.tx.clone();
        let ctx = ctx.clone();
        spawn_task(async move {
            let file = rfd::AsyncFileDialog::new()
                .set_title("Open NMR File or ZIP Archive")
                .add_filter(
                    "NMR Formats",
                    &["jdx", "dx", "jdf", "ft1", "ft2", "fid", "ucsf", "nmrML", "csv", "zip"],
                )
                .pick_file()
                .await;
            if let Some(f) = file {
                let name = f.file_name();
                let bytes = f.read().await;
                let _ = tx.send(AsyncLoadMsg::SingleFile { name, bytes });
                ctx.request_repaint();
            }
        });
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn spawn_open_folder_dialog(&self, ctx: &egui::Context) {
        let tx = self.tx.clone();
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            if let Some(path) = rfd::FileDialog::new()
                .set_title("Open Bruker / Varian / Spinsolve Directory")
                .pick_folder()
            {
                let name = path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("experiment")
                    .to_string();
                let mut tree = VirtualTree::new(name.clone());
                if let Ok(ds) = nmr_io::load_from_path(&path) {
                    // Quick wrapper via single virtual dataset isn't needed if we reuse load_from_path,
                    // here send tree or re-read:
                    drop(ds);
                }
                if let Err(e) = collect_folder(&path, &path, &mut tree, 3) {
                    let _ = tx.send(AsyncLoadMsg::Error(e.to_string()));
                } else {
                    let _ = tx.send(AsyncLoadMsg::VirtualDir(tree));
                }
                ctx.request_repaint();
            }
        });
    }

    fn first_visible_1d(&self) -> Option<&Signal1D> {
        self.eval.outputs.values().find_map(|v| match v {
            NodeOutput::Plot1D {
                signal,
                visible: true,
                ..
            } => Some(signal),
            _ => None,
        })
    }

    fn first_peaks(&self) -> Option<&nmr_core::PeakList> {
        self.eval.outputs.values().find_map(|v| match v {
            NodeOutput::Peaks(p) => Some(p),
            _ => None,
        })
    }

    fn add_node_centered(&mut self, node: NmrNode) {
        let offset = (self.snarl.node_ids().count() as f32 % 8.0) * 28.0;
        self.snarl.insert_node(egui::pos2(260.0 + offset, 180.0 + offset), node);
        self.dirty = true;
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn collect_folder(
    root: &std::path::Path,
    dir: &std::path::Path,
    tree: &mut VirtualTree,
    depth: usize,
) -> std::io::Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let p = entry.path();
        if p.is_dir() && depth > 0 {
            collect_folder(root, &p, tree, depth - 1)?;
        } else if p.is_file() {
            if let Ok(rel) = p.strip_prefix(root) {
                if let Ok(data) = std::fs::read(&p) {
                    tree.insert(rel.to_string_lossy().as_ref(), data);
                }
            }
        }
    }
    Ok(())
}

#[cfg(not(target_arch = "wasm32"))]
fn spawn_task<F: std::future::Future<Output = ()> + Send + 'static>(fut: F) {
    std::thread::spawn(move || {
        // Minimal single-future executor on native thread
        use std::task::{Context, Poll, Wake, Waker};
        struct ThreadWaker(std::thread::Thread);
        impl Wake for ThreadWaker {
            fn wake(self: std::sync::Arc<Self>) {
                self.0.unpark();
            }
        }
        let waker = Waker::from(std::sync::Arc::new(ThreadWaker(std::thread::current())));
        let mut cx = Context::from_waker(&waker);
        let mut fut = std::pin::pin!(fut);
        while let Poll::Pending = fut.as_mut().poll(&mut cx) {
            std::thread::park();
        }
    });
}

#[cfg(target_arch = "wasm32")]
fn spawn_task<F: std::future::Future<Output = ()> + 'static>(fut: F) {
    wasm_bindgen_futures::spawn_local(fut);
}

fn trigger_save_bytes(filename: String, bytes: Vec<u8>) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::thread::spawn(move || {
            if let Some(path) = rfd::FileDialog::new().set_file_name(&filename).save_file() {
                let _ = std::fs::write(path, bytes);
            }
        });
    }
    #[cfg(target_arch = "wasm32")]
    {
        use wasm_bindgen::JsCast;
        let Some(window) = web_sys::window() else { return };
        let Some(doc) = window.document() else { return };
        let arr = js_sys::Uint8Array::from(bytes.as_slice());
        let parts = js_sys::Array::new();
        parts.push(&arr.buffer());
        let opts = web_sys::BlobPropertyBag::new();
        opts.set_type("application/octet-stream");
        if let Ok(blob) = web_sys::Blob::new_with_u8_array_sequence_and_options(&parts, &opts) {
            if let Ok(url) = web_sys::Url::create_object_url_with_blob(&blob) {
                if let Ok(el) = doc.create_element("a") {
                    if let Ok(a) = el.dyn_into::<web_sys::HtmlAnchorElement>() {
                        a.set_href(&url);
                        a.set_download(&filename);
                        a.click();
                    }
                }
                let _ = web_sys::Url::revoke_object_url(&url);
            }
        }
    }
}

impl eframe::App for NmrReadApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.poll_async_loads();
        self.handle_dropped_files(&ctx);

        // Top menu bar
        egui::Panel::top("top_bar").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.strong("⚛ NMR-READ");
                ui.separator();

                ui.menu_button("File", |ui| {
                    if ui.button("Open File / ZIP...").clicked() {
                        self.spawn_open_file_dialog(&ctx);
                        ui.close();
                    }
                    #[cfg(not(target_arch = "wasm32"))]
                    if ui.button("Open Bruker / Varian / Spinsolve Folder...").clicked() {
                        self.spawn_open_folder_dialog(&ctx);
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("Save Workflow (.nmrflow.json)...").clicked() {
                        if let Ok(json) = serde_json::to_string_pretty(&self.snarl) {
                            trigger_save_bytes("workflow.nmrflow.json".into(), json.into_bytes());
                        }
                        ui.close();
                    }
                    if ui.button("Export Active Spectrum -> CSV").clicked() {
                        if let Some(sig) = self.first_visible_1d() {
                            trigger_save_bytes("spectrum.csv".into(), export_csv(sig).into_bytes());
                        }
                        ui.close();
                    }
                    if ui.button("Export Active Spectrum -> JCAMP-DX (.jdx)").clicked() {
                        if let Some(sig) = self.first_visible_1d() {
                            trigger_save_bytes("spectrum.jdx".into(), write_jcamp_dx(sig).into_bytes());
                        }
                        ui.close();
                    }
                    if ui.button("Export Active Spectrum -> NMRPipe (.ft1)").clicked() {
                        if let Some(sig) = self.first_visible_1d() {
                            trigger_save_bytes("spectrum.ft1".into(), write_nmrpipe_1d(sig));
                        }
                        ui.close();
                    }
                    if ui.button("Export Peaks -> CSV").clicked() {
                        if let Some(p) = self.first_peaks() {
                            trigger_save_bytes("peaks.csv".into(), export_peaks_csv(p).into_bytes());
                        }
                        ui.close();
                    }
                });

                ui.menu_button("Presets", |ui| {
                    if ui.button("1H Ethanol: Noisy vs Wavelet vs Cadzow").clicked() {
                        self.snarl = build_preset_snarl(WorkflowPreset::EthanolDenoiseCompare);
                        self.dirty = true;
                        ui.close();
                    }
                    if ui.button("1H Ethylbenzene: Apod + Auto-Phase + arPLS Baseline").clicked() {
                        self.snarl = build_preset_snarl(WorkflowPreset::PhaseAndBaseline);
                        self.dirty = true;
                        ui.close();
                    }
                    if ui.button("2D COSY FFT & Heatmap").clicked() {
                        self.snarl = build_preset_snarl(WorkflowPreset::Cosy2D);
                        self.dirty = true;
                        ui.close();
                    }
                });

                ui.separator();
                if ui.button("⟳ Re-run Graph").clicked() {
                    self.dirty = true;
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.weak(&self.status_msg);
                });
            });
        });

        // Left sidebar: Add Node Palette + Loaded Datasets
        egui::Panel::left("palette_panel")
            .resizable(true)
            .default_size(210.0)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Node Palette");
                    ui.weak("Click to add onto graph");
                    ui.separator();

                    ui.strong("Sources");
                    if ui.button("+ Synthetic 1D FID").clicked() {
                        self.add_node_centered(NmrNode::default_synthetic_ethanol());
                    }
                    if ui.button("+ Synthetic 2D COSY").clicked() {
                        self.add_node_centered(NmrNode::Synthetic2DCosy {
                            rows: 128,
                            cols: 128,
                            noise_sigma: 0.05,
                            seed: 12,
                        });
                    }
                    ui.add_space(6.0);

                    ui.strong("Conditioning");
                    if ui.button("+ Digital Filter Removal").clicked() {
                        self.add_node_centered(NmrNode::RemoveDigitalFilter {
                            use_override: false,
                            grpdly_override: 68.0,
                        });
                    }
                    if ui.button("+ DC Offset Removal").clicked() {
                        self.add_node_centered(NmrNode::DcOffset { tail_fraction: 0.15 });
                    }
                    if ui.button("+ Apodization Window").clicked() {
                        self.add_node_centered(NmrNode::default_apodization());
                    }
                    if ui.button("+ Zero Fill").clicked() {
                        self.add_node_centered(NmrNode::ZeroFill { multiplier: 2 });
                    }
                    ui.add_space(6.0);

                    ui.strong("Denoising");
                    if ui.button("+ Wavelet Denoise (SWT)").clicked() {
                        self.add_node_centered(NmrNode::default_wavelet());
                    }
                    if ui.button("+ Cadzow Hankel SVD").clicked() {
                        self.add_node_centered(NmrNode::default_cadzow());
                    }
                    if ui.button("+ Savitzky-Golay").clicked() {
                        self.add_node_centered(NmrNode::SavitzkyGolay {
                            half_window: 4,
                            poly_order: 2,
                        });
                    }
                    ui.add_space(6.0);

                    ui.strong("Transform & Spectral");
                    if ui.button("+ 1D FFT").clicked() {
                        self.add_node_centered(NmrNode::Fft1D {
                            half_first_point: true,
                        });
                    }
                    if ui.button("+ 1D IFFT").clicked() {
                        self.add_node_centered(NmrNode::Ifft1D);
                    }
                    if ui.button("+ Phase Correction").clicked() {
                        self.add_node_centered(NmrNode::PhaseCorrect {
                            phi0_deg: 0.0,
                            phi1_deg: 0.0,
                            pivot_frac: 0.5,
                            auto_entropy: true,
                            allow_first_order: false,
                        });
                    }
                    if ui.button("+ Baseline Correction").clicked() {
                        self.add_node_centered(NmrNode::BaselineCorrect {
                            method: BaselineMethod::ArPls,
                            lambda_exp: 6.0,
                            poly_degree: 3,
                        });
                    }
                    if ui.button("+ 2D FFT").clicked() {
                        self.add_node_centered(NmrNode::Fft2D { transform_f1: true });
                    }
                    ui.add_space(6.0);

                    ui.strong("Analysis & Output");
                    if ui.button("+ Peak Picking").clicked() {
                        self.add_node_centered(NmrNode::PeakPick {
                            snr_threshold: 6.0,
                            min_separation_pts: 4,
                            use_magnitude: false,
                        });
                    }
                    if ui.button("+ Compare SNR / SSIM").clicked() {
                        self.add_node_centered(NmrNode::CompareMetrics {
                            label: "Comparison".into(),
                        });
                    }
                    if ui.button("+ 1D Spectrum Plot").clicked() {
                        self.add_node_centered(NmrNode::SpectrumPlot1D {
                            label: "Spectrum".into(),
                            visible: true,
                            show_real: true,
                            show_imag: false,
                            show_magnitude: false,
                            show_baseline: false,
                        });
                    }
                    if ui.button("+ 2D Heatmap Plot").clicked() {
                        self.add_node_centered(NmrNode::HeatmapPlot2D {
                            label: "2D Map".into(),
                            visible: true,
                            use_magnitude: true,
                            floor_pct: 5.0,
                        });
                    }

                    ui.separator();
                    ui.strong("Supported Formats");
                    ui.small("Bruker (fid/ser/1r), Varian, JEOL (.jdf), Spinsolve (.1d/.dx), NMRPipe, Sparky (.ucsf), JCAMP-DX, nmrML, CSV, ZIP");
                    if !self.datasets_1d.is_empty() || !self.datasets_2d.is_empty() {
                        ui.separator();
                        ui.strong("Loaded Datasets");
                        for (k, s) in &self.datasets_1d {
                            ui.small(format!("1D: {k} ({} pts, {:?})", s.len(), s.meta.format));
                        }
                        for (k, s) in &self.datasets_2d {
                            ui.small(format!("2D: {k} ({}x{}, {:?})", s.rows, s.cols, s.meta.format));
                        }
                    }
                });
            });

        // Bottom dock: Interactive 1D/2D NMR Plots & Metrics
        egui::Panel::bottom("plot_dock")
            .resizable(true)
            .default_size(330.0)
            .min_size(180.0)
            .show(ui, |ui| {
                show_plot_dock(ui, &self.eval, &mut self.plot_state);
            });

        // Center: Interactive Node Graph Canvas
        egui::CentralPanel::default().show(ui, |ui| {
            let mut viewer = NmrSnarlViewer {
                dirty: &mut self.dirty,
                eval: &self.eval,
                dataset_keys_1d: self.datasets_1d.keys().cloned().collect(),
                dataset_keys_2d: self.datasets_2d.keys().cloned().collect(),
            };
            let style = SnarlStyle::new();
            SnarlWidget::new()
                .id_salt("nmr_read_snarl")
                .style(style)
                .show(&mut self.snarl, &mut viewer, ui);
        });

        if self.dirty {
            self.eval = evaluate_graph(&self.snarl, &self.datasets_1d, &self.datasets_2d);
            self.dirty = false;
        }
    }
}

struct NmrSnarlViewer<'a> {
    dirty: &'a mut bool,
    eval: &'a EvalResults,
    dataset_keys_1d: Vec<String>,
    dataset_keys_2d: Vec<String>,
}

impl SnarlViewer<NmrNode> for NmrSnarlViewer<'_> {
    fn title(&mut self, node: &NmrNode) -> String {
        node.title().to_string()
    }

    fn inputs(&mut self, node: &NmrNode) -> usize {
        node.input_pins().len()
    }

    fn outputs(&mut self, node: &NmrNode) -> usize {
        node.output_pins().len()
    }

    fn show_input(&mut self, pin: &InPin, ui: &mut Ui, snarl: &mut Snarl<NmrNode>) -> impl SnarlPin + 'static {
        let (label, kind) = snarl[pin.id.node].input_pins()[pin.id.input];
        ui.label(label);
        PinInfo::circle().with_fill(kind.color())
    }

    fn show_output(&mut self, pin: &OutPin, ui: &mut Ui, snarl: &mut Snarl<NmrNode>) -> impl SnarlPin + 'static {
        let (label, kind) = snarl[pin.id.node].output_pins()[pin.id.output];
        ui.label(label);
        PinInfo::circle().with_fill(kind.color())
    }

    fn has_body(&mut self, _node: &NmrNode) -> bool {
        true
    }

    fn show_body(
        &mut self,
        node_id: NodeId,
        _inputs: &[InPin],
        _outputs: &[OutPin],
        ui: &mut Ui,
        snarl: &mut Snarl<NmrNode>,
    ) {
        let mut changed = false;
        let node = &mut snarl[node_id];

        match node {
            NmrNode::Synthetic1D {
                preset,
                points,
                noise_sigma,
                ..
            } => {
                egui::ComboBox::from_id_salt(("syn_preset", node_id))
                    .selected_text(format!("{preset:?}"))
                    .show_ui(ui, |ui| {
                        for p in [
                            SyntheticPreset::Ethanol,
                            SyntheticPreset::Ethylbenzene,
                            SyntheticPreset::TwoPeaks,
                        ] {
                            if ui.selectable_value(preset, p, format!("{p:?}")).changed() {
                                changed = true;
                            }
                        }
                    });
                changed |= ui
                    .add(egui::Slider::new(points, 128..=8192).text("Points"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(noise_sigma, 0.0..=1.5).text("Noise"))
                    .changed();
            }
            NmrNode::Synthetic2DCosy {
                rows,
                cols,
                noise_sigma,
                ..
            } => {
                changed |= ui.add(egui::Slider::new(rows, 32..=256).text("Rows")).changed();
                changed |= ui.add(egui::Slider::new(cols, 32..=256).text("Cols")).changed();
                changed |= ui
                    .add(egui::Slider::new(noise_sigma, 0.0..=0.5).text("Noise"))
                    .changed();
            }
            NmrNode::LoadedDataset1D { dataset_key } => {
                egui::ComboBox::from_id_salt(("ds1d", node_id))
                    .selected_text(if dataset_key.is_empty() {
                        "Select..."
                    } else {
                        dataset_key.as_str()
                    })
                    .show_ui(ui, |ui| {
                        for k in &self.dataset_keys_1d {
                            if ui.selectable_value(dataset_key, k.clone(), k).changed() {
                                changed = true;
                            }
                        }
                    });
            }
            NmrNode::LoadedDataset2D { dataset_key } => {
                egui::ComboBox::from_id_salt(("ds2d", node_id))
                    .selected_text(if dataset_key.is_empty() {
                        "Select..."
                    } else {
                        dataset_key.as_str()
                    })
                    .show_ui(ui, |ui| {
                        for k in &self.dataset_keys_2d {
                            if ui.selectable_value(dataset_key, k.clone(), k).changed() {
                                changed = true;
                            }
                        }
                    });
            }
            NmrNode::RemoveDigitalFilter {
                use_override,
                grpdly_override,
            } => {
                changed |= ui.checkbox(use_override, "Override GRPDLY").changed();
                if *use_override {
                    changed |= ui
                        .add(egui::DragValue::new(grpdly_override).speed(0.1).prefix("pts: "))
                        .changed();
                }
            }
            NmrNode::DcOffset { tail_fraction } => {
                changed |= ui
                    .add(egui::Slider::new(tail_fraction, 0.02..=0.5).text("Tail"))
                    .changed();
            }
            NmrNode::Apodization { config } => {
                egui::ComboBox::from_id_salt(("apod", node_id))
                    .selected_text(format!("{:?}", config.window))
                    .show_ui(ui, |ui| {
                        for w in [
                            WindowType::Exponential,
                            WindowType::LorentzToGauss,
                            WindowType::SineBell,
                            WindowType::Hanning,
                            WindowType::Hamming,
                        ] {
                            if ui
                                .selectable_value(&mut config.window, w, format!("{w:?}"))
                                .changed()
                            {
                                changed = true;
                            }
                        }
                    });
                changed |= ui
                    .add(egui::Slider::new(&mut config.lb_hz, -5.0..=15.0).text("LB (Hz)"))
                    .changed();
            }
            NmrNode::ZeroFill { multiplier } => {
                changed |= ui
                    .add(egui::Slider::new(multiplier, 1..=8).text("x Fill"))
                    .changed();
            }
            NmrNode::WaveletDenoise { config } => {
                egui::ComboBox::from_id_salt(("wav_fam", node_id))
                    .selected_text(format!("{:?}", config.family))
                    .show_ui(ui, |ui| {
                        for f in [
                            WaveletFamily::Haar,
                            WaveletFamily::Db2,
                            WaveletFamily::Db4,
                            WaveletFamily::Db6,
                            WaveletFamily::Sym4,
                        ] {
                            if ui
                                .selectable_value(&mut config.family, f, format!("{f:?}"))
                                .changed()
                            {
                                changed = true;
                            }
                        }
                    });
                egui::ComboBox::from_id_salt(("wav_mode", node_id))
                    .selected_text(format!("{:?}", config.mode))
                    .show_ui(ui, |ui| {
                        for m in [ThresholdMode::Soft, ThresholdMode::Hard, ThresholdMode::Garrote] {
                            if ui
                                .selectable_value(&mut config.mode, m, format!("{m:?}"))
                                .changed()
                            {
                                changed = true;
                            }
                        }
                    });
                changed |= ui
                    .add(egui::Slider::new(&mut config.levels, 1..=7).text("Levels"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(&mut config.threshold_scale, 0.1..=3.0).text("Scale"))
                    .changed();
            }
            NmrNode::CadzowDenoise { config } => {
                changed |= ui
                    .add(egui::Slider::new(&mut config.rank, 1..=64).text("Rank"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(&mut config.iterations, 1..=8).text("Iters"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(&mut config.pencil_fraction, 0.15..=0.5).text("Pencil"))
                    .changed();
            }
            NmrNode::SavitzkyGolay {
                half_window,
                poly_order,
            } => {
                changed |= ui
                    .add(egui::Slider::new(half_window, 1..=20).text("HalfWin"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(poly_order, 1..=4).text("Degree"))
                    .changed();
            }
            NmrNode::Fft1D { half_first_point } => {
                changed |= ui.checkbox(half_first_point, "Half 1st point").changed();
            }
            NmrNode::Ifft1D => {}
            NmrNode::PhaseCorrect {
                phi0_deg,
                phi1_deg,
                auto_entropy,
                allow_first_order,
                ..
            } => {
                changed |= ui.checkbox(auto_entropy, "Auto (Entropy)").changed();
                if *auto_entropy {
                    changed |= ui.checkbox(allow_first_order, "Allow PH1").changed();
                } else {
                    changed |= ui
                        .add(egui::Slider::new(phi0_deg, -180.0..=180.0).text("PH0"))
                        .changed();
                    changed |= ui
                        .add(egui::Slider::new(phi1_deg, -180.0..=180.0).text("PH1"))
                        .changed();
                }
            }
            NmrNode::BaselineCorrect {
                method,
                lambda_exp,
                poly_degree,
            } => {
                egui::ComboBox::from_id_salt(("bl_method", node_id))
                    .selected_text(format!("{method:?}"))
                    .show_ui(ui, |ui| {
                        for m in [BaselineMethod::ArPls, BaselineMethod::Polynomial] {
                            if ui.selectable_value(method, m, format!("{m:?}")).changed() {
                                changed = true;
                            }
                        }
                    });
                match method {
                    BaselineMethod::ArPls => {
                        changed |= ui
                            .add(egui::Slider::new(lambda_exp, 2.0..=8.0).text("log10(λ)"))
                            .changed();
                    }
                    BaselineMethod::Polynomial => {
                        changed |= ui
                            .add(egui::Slider::new(poly_degree, 1..=7).text("Degree"))
                            .changed();
                    }
                }
            }
            NmrNode::Fft2D { transform_f1 } => {
                changed |= ui.checkbox(transform_f1, "Transform F1").changed();
            }
            NmrNode::PeakPick {
                snr_threshold,
                min_separation_pts,
                use_magnitude,
            } => {
                changed |= ui
                    .add(egui::Slider::new(snr_threshold, 2.0..=30.0).text("SNR Min"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(min_separation_pts, 1..=32).text("Sep pts"))
                    .changed();
                changed |= ui.checkbox(use_magnitude, "Magnitude").changed();
            }
            NmrNode::CompareMetrics { label } => {
                ui.text_edit_singleline(label);
                if let Some(NodeOutput::Metrics(m)) = self.eval.outputs.get(&node_id) {
                    ui.colored_label(
                        Color32::from_rgb(0x81, 0xc7, 0x84),
                        format!("SNR: {:.1} dB", m.snr_db),
                    );
                    if let Some(ssim) = m.ssim {
                        ui.small(format!("SSIM: {ssim:.4} | RMSE: {:.4}", m.rmse.unwrap_or(0.0)));
                    }
                }
            }
            NmrNode::SpectrumPlot1D {
                label,
                visible,
                show_real,
                show_imag,
                show_magnitude,
                show_baseline,
            } => {
                ui.text_edit_singleline(label);
                ui.horizontal(|ui| {
                    ui.checkbox(visible, "Vis");
                    ui.checkbox(show_real, "Re");
                    ui.checkbox(show_imag, "Im");
                    ui.checkbox(show_magnitude, "|M|");
                    ui.checkbox(show_baseline, "BL");
                });
            }
            NmrNode::HeatmapPlot2D {
                label,
                visible,
                use_magnitude,
                floor_pct,
            } => {
                ui.text_edit_singleline(label);
                ui.checkbox(visible, "Visible");
                ui.checkbox(use_magnitude, "Magnitude");
                ui.add(egui::Slider::new(floor_pct, 0.0..=50.0).text("Floor %"));
            }
        }

        if let Some(err) = self.eval.errors.get(&node_id) {
            ui.colored_label(Color32::LIGHT_RED, err);
        } else if let Some(ms) = self.eval.timings_ms.get(&node_id) {
            ui.weak(format!("{ms:.2} ms"));
        }

        if changed {
            *self.dirty = true;
        }
    }

    fn has_node_menu(&mut self, _node: &NmrNode) -> bool {
        true
    }

    fn show_node_menu(
        &mut self,
        node: NodeId,
        _inputs: &[InPin],
        _outputs: &[OutPin],
        ui: &mut Ui,
        snarl: &mut Snarl<NmrNode>,
    ) {
        if ui.button("Delete Node").clicked() {
            snarl.remove_node(node);
            *self.dirty = true;
            ui.close();
        }
    }

    fn connect(&mut self, from: &OutPin, to: &InPin, snarl: &mut Snarl<NmrNode>) {
        let out_kind = snarl[from.id.node].output_pins()[from.id.output].1;
        let in_kind = snarl[to.id.node].input_pins()[to.id.input].1;
        if out_kind != in_kind {
            return;
        }
        for &remote in &to.remotes {
            snarl.disconnect(remote, to.id);
        }
        snarl.connect(from.id, to.id);
        *self.dirty = true;
    }

    fn disconnect(&mut self, from: &OutPin, to: &InPin, snarl: &mut Snarl<NmrNode>) {
        snarl.disconnect(from.id, to.id);
        *self.dirty = true;
    }

    fn drop_inputs(&mut self, pin: &InPin, snarl: &mut Snarl<NmrNode>) {
        snarl.drop_inputs(pin.id);
        *self.dirty = true;
    }

    fn drop_outputs(&mut self, pin: &OutPin, snarl: &mut Snarl<NmrNode>) {
        snarl.drop_outputs(pin.id);
        *self.dirty = true;
    }
}
