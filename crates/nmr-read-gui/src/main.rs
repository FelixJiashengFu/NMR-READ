#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

#[cfg(not(target_arch = "wasm32"))]
use nmr_read_gui::NmrReadApp;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result<()> {
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1420.0, 900.0])
            .with_min_inner_size([960.0, 620.0])
            .with_title("NMR-READ — Pure-Rust NMR Denoising & Analysis Studio"),
        ..Default::default()
    };
    eframe::run_native(
        "NMR-READ",
        native_options,
        Box::new(|cc| Ok(Box::new(NmrReadApp::new(cc)))),
    )
}

#[cfg(target_arch = "wasm32")]
fn main() {
    let _ = nmr_read_gui::start_wasm();
}
