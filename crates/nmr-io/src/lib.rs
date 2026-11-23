//! Multi-vendor NMR format readers and exporters in pure Rust (Desktop + WebAssembly compatible).
//!
//! Supported formats:
//! - **Bruker TopSpin** (`acqus`/`acqu2s`/`procs`, `fid`/`ser`, `pdata/1/1r`+`1i`, `GRPDLY`/`DECIM`+`DSPFVS`)
//! - **JCAMP-DX** (`.dx`/`.jdx`, including AFFN, ASDF SQZ/DIF/DUP decompression, and `NTUPLES`/`PAGE`)
//! - **Varian / Agilent** (`procpar` + big-endian `fid`)
//! - **JEOL Delta** (`.jdf` binary header, parameter block, submatrix reordering)
//! - **NMRPipe** (`.fid`/`.ft1`/`.ft2` 512×4-byte `fdata` header + 1D/2D real/complex arrays)
//! - **Magritek Spinsolve** (`acqu.par`/`proc.par` + `data.1d`/`spectrum.1d`/`nmr_fid.dx`)
//! - **Sparky UCSF** (`.ucsf` 2D tiled spectrum)
//! - **nmrML** (`.nmrML` XML + base64 arrays)
//! - **ZIP experiment archives** (works identically in browser Wasm and desktop)
//! - **CSV / TSV / JSON** export

pub mod archive;
pub(crate) mod binread;
pub mod bruker;
pub mod error;
pub mod export;
pub mod jcamp;
pub mod jeol;
pub mod nmrml;
pub mod pipe;
pub mod sparky;
pub mod spinsolve;
pub mod varian;
pub mod virtual_fs;

pub use error::{NmrIoError, Result};
pub use virtual_fs::{LoadedDataset, VirtualTree, load_from_bytes, load_from_virtual_tree};

#[cfg(not(target_arch = "wasm32"))]
pub use virtual_fs::load_from_path;
