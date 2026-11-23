pub mod app;
pub mod engine;
pub mod nodes;
pub mod plot;

pub use app::NmrReadApp;

#[cfg(target_arch = "wasm32")]
pub fn start_wasm() -> Result<(), wasm_bindgen::JsValue> {
    use wasm_bindgen::JsCast;
    let web_options = eframe::WebOptions::default();
    let window = web_sys::window().ok_or_else(|| wasm_bindgen::JsValue::from_str("No window"))?;
    let document = window
        .document()
        .ok_or_else(|| wasm_bindgen::JsValue::from_str("No document"))?;
    let canvas = document
        .get_element_by_id("nmr_read_canvas")
        .ok_or_else(|| wasm_bindgen::JsValue::from_str("Missing #nmr_read_canvas"))?
        .dyn_into::<web_sys::HtmlCanvasElement>()?;

    wasm_bindgen_futures::spawn_local(async move {
        let runner = eframe::WebRunner::new();
        let res = runner
            .start(
                canvas,
                web_options,
                Box::new(|cc| Ok(Box::new(NmrReadApp::new(cc)))),
            )
            .await;
        if res.is_ok() {
            if let Some(win) = web_sys::window() {
                if let Some(doc) = win.document() {
                    if let Some(el) = doc.get_element_by_id("loading") {
                        el.remove();
                    }
                }
            }
        }
    });
    Ok(())
}
