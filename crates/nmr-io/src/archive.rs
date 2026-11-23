use crate::error::Result;
use crate::virtual_fs::{LoadedDataset, VirtualTree};

/// Convenience helper for loading a zipped NMR dataset buffer.
pub fn load_zip_bytes(name: &str, bytes: &[u8]) -> Result<LoadedDataset> {
    let tree = VirtualTree::from_zip_bytes(name, bytes)?;
    crate::virtual_fs::load_from_virtual_tree(&tree)
}
