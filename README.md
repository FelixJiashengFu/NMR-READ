<div align="center">

```text
 _   _ __  __ ____       ____  _____    _    ____
| \ | |  \/  |  _ \     |  _ \| ____|  / \  |  _ \
|  \| | |\/| | |_) |____| |_) |  _|   / _ \ | | | |
| |\  | |  | |  _ <_____|  _ <| |___ / ___ \| |_| |
|_| \_|_|  |_|_| \_\    |_| \_\_____/_/   \_\____/
```

<h3>Pure-Rust Cross-Platform Multi-Vendor NMR Reader, Node-Graph Processor &amp; Denoising Studio</h3>

<p>
  <strong>Zero C / Fortran FFI</strong> &bull;
  <strong>Native Desktop (macOS / Linux / Windows)</strong> &bull;
  <strong>Browser WebAssembly</strong> &bull;
  <strong>MIT OR Apache-2.0</strong>
</p>

</div>

---

## Why `NMR-READ`?

`NMR-READ` is a clean-room **pure-Rust** NMR data workbench inspired by node-graph NMR denoising workflows and multi-vendor spectroscopy parsers across the open-source NMR ecosystem (`nmrglue`, `cheminfo/nmr-parser`, `loom`, `NMRTools.jl`, `SpectroChemPy`, `spinsolveproc`, `nmrML`, and `0ffz/node-nmr`).

By eliminating external dynamic libraries (`FFTW3`, `wavelib`, Fortran `PROPACK`, `OpenBLAS`), **a single Rust codebase compiles natively across Apple Silicon, Linux, Windows, and WebAssembly (`wasm32-unknown-unknown`)** with zero runtime dependencies.

---

## Key Capabilities & Improvements

### 1. Multi-Vendor NMR Format Reader (`nmr-io`)
Works identically from native disk folders, single files, or **in-memory `.zip` archives** (including browser drag-and-drop in WebAssembly):

| Vendor / Format | Files / Extensions | Capabilities |
| :--- | :--- | :--- |
| **Bruker TopSpin** | `acqus`, `acqu2s`, `procs`, `fid`, `ser`, `pdata/1/1r`+`1i` | 1D & 2D raw/processed data, `BYTORDP`/`DTYPA`/`NC_proc`, `GRPDLY` & `DSPFVS`/`DECIM` digital filter group delay |
| **JCAMP-DX** | `.dx`, `.jdx` | AFFN + compressed **ASDF** (`SQZ` / `DIF` / `DUP`) decoding, `NTUPLES`/`PAGE` complex spectra, plus exporter |
| **Varian / Agilent** | `procpar` + `fid` | Big-endian 32-byte file & 28-byte block headers (`i16`, `i32`, `f32`), 1D & 2D |
| **JEOL Delta** | `.jdf` | Big/little-endian (`JEOL.NMR` / `RMN.LOEJ`), parameter block, 1D & 2D submatrix reordering |
| **NMRPipe** | `.fid`, `.ft1`, `.ft2` | Auto-endian 2048-byte `fdata` header, 1D & 2D real/complex unappended arrays, plus 1D exporter |
| **Magritek Spinsolve** | `acqu.par` + `data.1d` / `spectrum.1d` / `nmr_fid.dx` | Benchtop 1D binary & JCAMP-DX import |
| **Sparky UCSF** | `.ucsf` | Big-endian 180B global + 128B axis headers + tiled 2D matrix |
| **nmrML** | `.nmrML` | Open XML schema + Base64 (`complex128` / `complex64` / `float64`) arrays |
| **CSV / TSV / Workflow** | `.csv`, `.tsv`, `.nmrflow.json` | Tabular import/export, peak table CSV export, node-graph workflow save/load |

### 2. Pure-Rust DSP, Denoising & Analysis (`nmr-dsp` & `nmr-core`)
- **FFT / IFFT & 2D FFT**: `rustfft` with standard NMR downfield-left / upfield-right spectral ordering, half-dwell first-point scaling, zero filling, and FID DC offset removal.
- **Bruker Digital Filter Removal**: Exact Fourier linear phase shift + tail trim for `GRPDLY`.
- **Apodization & Smoothing**: Exponential (`EM`), Lorentz-to-Gauss (`GM`), Shifted Sine-Bell (`SSB`), Hanning, Hamming, and Savitzky-Golay FIR polynomial smoothing.
- **Manual & Auto Phase Correction**: 0th/1st-order manual phasing + entropy-minimization auto-phasing.
- **Baseline Correction**: $O(N)$ banded symmetric pentadiagonal **arPLS** (Asymmetrically Reweighted Penalized Least Squares) + iterative polynomial baseline fit.
- **Wavelet Denoising (SWT)**: Symmetric mirror-extended 1D Stationary Wavelet Transform (`Haar`, `Db2`, `Db4`, `Db6`, `Sym4`) with Soft / Hard / Garrote VisuShrink thresholding.
- **Cadzow Hankel Low-Rank SVD**: $O(N \log N)$ FFT-accelerated Hankel matrix-vector operator + Golub-Kahan-Lanczos bidiagonalization + `faer` thin SVD + FFT antidiagonal averaging.
- **Peak Picking, Integration & Quality Metrics**: Parabolic sub-bin interpolation, FWHM estimation, region integration, MAD noise estimation, SNR (dB), RMSE, PSNR, and 1D SSIM.
- **Synthetic NMR Generator**: Built-in `1H` Ethanol, Ethylbenzene, Two-Peak, and 2D COSY FID generators for instant testing and benchmarking.

### 3. Node-Graph Workstation (`nmr-read-gui`) & CLI (`nmr-read-cli`)
- Interactive DAG editor powered by `egui` + `egui-snarl` + `egui_plot`.
- Multi-spectrum 1D overlay with chemistry-standard descending ppm axis, peak markers, and baseline overlays.
- 2D magnitude/real heatmap viewer for 2D NMR matrices.
- Built-in one-click presets & drag-and-drop file/ZIP loading on Desktop and Browser WebAssembly.

---

## Workspace Structure

```text
NMR-READ/
├── crates/
│   ├── nmr-core/       # Signal1D, Signal2D, AxisMeta, PeakList, SSIM/SNR metrics
│   ├── nmr-io/         # Multi-vendor readers (Bruker, Varian, JEOL, Pipe, Spinsolve, Sparky, JCAMP-DX, nmrML, ZIP)
│   ├── nmr-dsp/        # Pure-Rust FFT, Cadzow SVD, SWT Wavelet, arPLS baseline, auto-phase, peak picker, synthetic FIDs
│   ├── nmr-read-cli/   # Headless CLI (inspect, process, convert, benchmark)
│   └── nmr-read-gui/   # Cross-platform Desktop + WebAssembly node-graph UI
├── Trunk.toml          # WebAssembly bundle configuration
└── index.html          # Browser Wasm canvas shell
```

---

## Quick Start

### Run Unit Tests
```bash
cargo test --workspace
```

### Run Desktop GUI
```bash
cargo run --release -p nmr-read-gui
```

### Run CLI Tools & Denoising Benchmark
```bash
# Built-in 1H Ethanol denoising benchmark (Wavelet vs Cadzow)
cargo run --release -p nmr-read-cli -- benchmark --preset ethanol --points 1024 --noise 0.25

# Inspect any Bruker/Varian/Spinsolve folder, .jdf, .ft1/.ft2, .jdx, .ucsf, .nmrML, or .zip
cargo run --release -p nmr-read-cli -- inspect path/to/experiment

# Process & export CSV + peaks
cargo run --release -p nmr-read-cli -- process path/to/bruker_exp --remove-df --fft --auto-phase --baseline-lambda 1000000 --cadzow-rank 16 --out-csv spec.csv --out-peaks peaks.csv

# Convert between formats (CSV, JCAMP-DX .jdx, NMRPipe .ft1)
cargo run --release -p nmr-read-cli -- convert path/to/experiment --to jdx --output spectrum.jdx
```

### Build / Serve WebAssembly (Browser Edition)
```bash
rustup target add wasm32-unknown-unknown
cargo install trunk
trunk serve --release
```

---

## License

Dual-licensed under either of:
- [MIT License](LICENSE-MIT)
- [Apache License, Version 2.0](LICENSE-APACHE)

at your option.
