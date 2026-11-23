use crate::error::{NmrIoError, Result};
use crate::{bruker, jcamp, jeol, nmrml, pipe, sparky, spinsolve, varian};
use nmr_core::{Signal1D, Signal2D};
use std::collections::BTreeMap;
use std::io::Read;

/// Unified loaded NMR dataset (either 1D or 2D).
#[derive(Debug, Clone)]
pub enum LoadedDataset {
    D1(Signal1D),
    D2(Signal2D),
}

impl LoadedDataset {
    pub fn title(&self) -> &str {
        match self {
            LoadedDataset::D1(s) => &s.meta.title,
            LoadedDataset::D2(s) => &s.meta.title,
        }
    }

    pub fn into_1d(self) -> Result<Signal1D> {
        match self {
            LoadedDataset::D1(s) => Ok(s),
            LoadedDataset::D2(_) => Err(NmrIoError::Malformed {
                format: "Dataset",
                detail: "Expected 1D dataset, got 2D".to_string(),
            }),
        }
    }
}

/// Platform-independent in-memory directory tree (normalized `/` paths, lowercase lookup helper).
#[derive(Debug, Clone, Default)]
pub struct VirtualTree {
    pub root_name: String,
    pub files: BTreeMap<String, Vec<u8>>,
}

impl VirtualTree {
    pub fn new(root_name: impl Into<String>) -> Self {
        Self {
            root_name: root_name.into(),
            files: BTreeMap::new(),
        }
    }

    pub fn insert(&mut self, rel_path: impl AsRef<str>, bytes: Vec<u8>) {
        let norm = normalize_rel_path(rel_path.as_ref());
        if !norm.is_empty() {
            self.files.insert(norm, bytes);
        }
    }

    /// Find a file by exact relative path or common nested prefix (strips leading single wrapper folder in ZIPs).
    pub fn get(&self, rel: &str) -> Option<&[u8]> {
        let norm = normalize_rel_path(rel);
        if let Some(b) = self.files.get(&norm) {
            return Some(b.as_slice());
        }
        // Fallback: match `<wrapper>/<norm>` if user zipped a top-level folder
        let suffix = format!("/{norm}");
        self.files
            .iter()
            .find(|(k, _)| k.eq_ignore_ascii_case(&norm) || k.to_ascii_lowercase().ends_with(&suffix.to_ascii_lowercase()))
            .map(|(_, v)| v.as_slice())
    }

    pub fn contains(&self, rel: &str) -> bool {
        self.get(rel).is_some()
    }

    /// Unpack a ZIP archive into a [`VirtualTree`] entirely in memory.
    pub fn from_zip_bytes(name: impl Into<String>, zip_bytes: &[u8]) -> Result<Self> {
        let reader = std::io::Cursor::new(zip_bytes);
        let mut archive = zip::ZipArchive::new(reader)?;
        let mut tree = Self::new(name);
        for i in 0..archive.len() {
            let mut file = archive.by_index(i)?;
            if file.is_dir() {
                continue;
            }
            let raw_name = file.name()?.to_string();
            // Skip macOS AppleDouble / metadata noise
            if raw_name.contains("__MACOSX") || raw_name.rsplit('/').next().is_some_and(|f| f.starts_with("._")) {
                continue;
            }
            let mut buf = Vec::with_capacity(file.size() as usize);
            file.read_to_end(&mut buf)?;
            tree.insert(&raw_name, buf);
        }
        Ok(tree)
    }
}

fn normalize_rel_path(p: &str) -> String {
    p.replace('\\', "/")
        .trim_start_matches("./")
        .trim_matches('/')
        .to_string()
}

/// Automatically detect format and load from a [`VirtualTree`].
pub fn load_from_virtual_tree(tree: &VirtualTree) -> Result<LoadedDataset> {
    // 1. Bruker TopSpin directory
    if tree.contains("acqus") || tree.contains("pdata/1/procs") || tree.contains("pdata/1/1r") {
        return bruker::load_bruker_tree(tree);
    }
    // 2. Varian / Agilent directory
    if tree.contains("procpar") && tree.contains("fid") {
        return varian::load_varian_tree(tree);
    }
    // 3. Magritek Spinsolve directory
    if tree.contains("acqu.par")
        || tree.contains("data.1d")
        || tree.contains("spectrum.1d")
        || tree.contains("spectrum_processed.1d")
    {
        return spinsolve::load_spinsolve_tree(tree);
    }
    // 4. Single-file inside tree fallback
    if let Some((name, bytes)) = tree.files.iter().next() {
        return load_from_bytes(name, bytes);
    }
    Err(NmrIoError::UnsupportedFormat(tree.root_name.clone()))
}

/// Automatically detect format from filename + magic bytes and load from an in-memory buffer.
/// Works seamlessly on both Desktop and WebAssembly.
pub fn load_from_bytes(filename: &str, bytes: &[u8]) -> Result<LoadedDataset> {
    let lower = filename.to_ascii_lowercase();

    // ZIP bundle containing an experiment directory (Bruker / Varian / Spinsolve)
    if bytes.starts_with(b"PK\x03\x04") || lower.ends_with(".zip") {
        let tree = VirtualTree::from_zip_bytes(filename, bytes)?;
        return load_from_virtual_tree(&tree);
    }

    // JEOL Delta (.jdf): "JEOL.NMR" or "RMN.LOEJ"
    if bytes.len() >= 16 && (&bytes[0..8] == b"JEOL.NMR" || &bytes[0..8] == b"RMN.LOEJ") {
        return jeol::load_jdf_bytes(filename, bytes);
    }

    // Sparky UCSF (.ucsf): starts with "UCSF NMR"
    if bytes.len() >= 180 && &bytes[0..8] == b"UCSF NMR" {
        return sparky::load_ucsf_bytes(filename, bytes).map(LoadedDataset::D2);
    }

    // NMRPipe (2048-byte header)
    if pipe::looks_like_nmrpipe(bytes) || lower.ends_with(".ft1") || lower.ends_with(".ft2") || lower.ends_with(".fid") {
        if let Ok(ds) = pipe::load_nmrpipe_bytes(filename, bytes) {
            return Ok(ds);
        }
    }

    // Text-based formats: JCAMP-DX, nmrML, CSV/TSV
    if let Ok(text) = std::str::from_utf8(bytes) {
        let head = &text[..text.len().min(4096)];
        if head.contains("<nmrML") || lower.ends_with(".nmrml") {
            return nmrml::load_nmrml_str(filename, text).map(LoadedDataset::D1);
        }
        if head.contains("##JCAMP-DX") || head.contains("##XYDATA") || lower.ends_with(".dx") || lower.ends_with(".jdx") {
            return jcamp::load_jcamp_str(filename, text).map(LoadedDataset::D1);
        }
        if lower.ends_with(".csv") || lower.ends_with(".tsv") || lower.ends_with(".txt") {
            return crate::export::load_csv_str(filename, text).map(LoadedDataset::D1);
        }
    }

    Err(NmrIoError::UnsupportedFormat(filename.to_string()))
}

/// Load a file or experiment directory from the native filesystem.
#[cfg(not(target_arch = "wasm32"))]
pub fn load_from_path(path: impl AsRef<std::path::Path>) -> Result<LoadedDataset> {
    let path = path.as_ref();
    if path.is_dir() {
        let mut tree = VirtualTree::new(
            path.file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("dataset"),
        );
        collect_dir_files(path, path, &mut tree, 3)?;
        return load_from_virtual_tree(&tree);
    }
    let bytes = std::fs::read(path)?;
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("dataset");
    load_from_bytes(name, &bytes)
}

#[cfg(not(target_arch = "wasm32"))]
fn collect_dir_files(
    root: &std::path::Path,
    dir: &std::path::Path,
    tree: &mut VirtualTree,
    max_depth: usize,
) -> Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let p = entry.path();
        if p.is_dir() {
            if max_depth > 0 {
                collect_dir_files(root, &p, tree, max_depth - 1)?;
            }
        } else if p.is_file() {
            if let Ok(rel) = p.strip_prefix(root) {
                let rel_str = rel.to_string_lossy();
                // Skip giant unused sub-files unless relevant
                if let Ok(meta) = entry.metadata() {
                    if meta.len() <= 256 * 1024 * 1024 {
                        if let Ok(data) = std::fs::read(&p) {
                            tree.insert(rel_str.as_ref(), data);
                        }
                    }
                }
            }
        }
    }
    Ok(())
}
