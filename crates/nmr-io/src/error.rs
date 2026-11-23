use thiserror::Error;

pub type Result<T> = std::result::Result<T, NmrIoError>;

#[derive(Debug, Error)]
pub enum NmrIoError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Unrecognized or unsupported NMR dataset format: {0}")]
    UnsupportedFormat(String),

    #[error("Missing required file '{0}' in dataset")]
    MissingFile(String),

    #[error("Missing required parameter '{0}'")]
    MissingParameter(String),

    #[error("Malformed header or data in {format}: {detail}")]
    Malformed {
        format: &'static str,
        detail: String,
    },

    #[error("ZIP archive error: {0}")]
    Zip(#[from] zip::result::ZipError),

    #[error("XML parse error: {0}")]
    Xml(String),

    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
}
