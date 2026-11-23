pub mod app;
pub mod engine;
pub mod nodes;
pub mod plot;

pub use app::NmrReadApp;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub fn start_wasm() -> Result<(), JsValue> {
    use wasm_bindgen::JsCast;
    let web_options = eframe::WebOptions::default();
    let window = web_sys::window().ok_or_else(|| JsValue::from_str("No window"))?;
    let document = window.document().ok_or_else(|| JsValue::from_str("No document"))?;
    let canvas = document
        .get_element_by_id("nmr_read_canvas")
        .ok_or_else(|| JsValue::from_str("Missing #nmr_read_canvas"))?
        .dyn_into::<web_sys::HtmlCanvasElement>()?;

    wasm_bindgen_futures::spawn_local(async move {
        let runner = eframe::WebRunner::new();
        let _ = runner
            .start(
                canvas,
                web_options,
                Box::new(|cc| Ok(Box::new(NmrReadApp::new(cc)))),
            )
            .await;
    });
    Ok(())
}
