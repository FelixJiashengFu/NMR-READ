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
  <a href="https://felixjiashengfu.github.io/NMR-READ/"><strong>🌐 Live Web App (WebAssembly)</strong></a>
  &nbsp;&bull;&nbsp;
  <a href="#english">English Guide</a>
  &nbsp;&bull;&nbsp;
  <a href="#简体中文">简体中文指南</a>
  &nbsp;&bull;&nbsp;
  <a href="https://github.com/FelixJiashengFu/NMR-READ/actions">CI &amp; Builds</a>
</p>

<p>
  <strong>Zero C / Fortran FFI</strong> &bull;
  <strong>Desktop (macOS / Linux / Windows)</strong> &bull;
  <strong>Browser WebAssembly</strong> &bull;
  <strong>MIT OR Apache-2.0</strong>
</p>

</div>

---

## English

### Overview

**`NMR-READ`** is a clean-room, pure-Rust NMR data workbench inspired by node-graph processing workflows and open-source spectroscopy readers (`nmrglue`, `cheminfo/nmr-parser`, `loom`, `NMRTools.jl`, `SpectroChemPy`, `spinsolveproc`, `nmrML`, and `0ffz/node-nmr`).

By implementing all transforms, wavelets, SVD denoising, and vendor parsers directly in Rust—without external dynamic libraries (`FFTW3`, `wavelib`, Fortran `PROPACK`, or `OpenBLAS`)—one codebase compiles natively across **macOS (Apple Silicon & Intel)**, **Linux**, **Windows**, and **WebAssembly (`wasm32-unknown-unknown`)**.

---

### How to Use

#### 1. Browser Edition (Zero Install — WebAssembly)

Open directly in your browser: **https://felixjiashengfu.github.io/NMR-READ/**

1. **Try Built-in Presets**: Open **`Presets`** in the top menu bar to switch between:
   - `1H Ethanol: Noisy vs Wavelet vs Cadzow` — compares raw noisy FFT, Sym4 Wavelet denoising, and Cadzow Hankel SVD denoising side by side with SNR / RMSE / SSIM metrics.
   - `1H Ethylbenzene: Apod + Auto-Phase + arPLS Baseline` — full 1D conditioning, entropy auto-phasing, arPLS baseline subtraction, and peak picking.
   - `2D COSY FFT & Heatmap` — synthetic 2D time-domain matrix, 2D FFT, and interactive 2D heatmap.
2. **Load Your Own NMR Data**:
   - **Drag & Drop** a single NMR file (`.jdx`, `.dx`, `.jdf`, `.ft1`, `.ft2`, `.ucsf`, `.nmrML`, `.csv`), multiple files from an experiment folder, or a **`.zip` archive** of a Bruker / Varian / Spinsolve folder directly onto the canvas.
   - Or click **`File → Open File / ZIP...`**.
   - A new **`Loaded 1D Dataset`** or **`Loaded 2D Dataset`** source node will appear on the graph.
3. **Build & Tune the Node Graph**:
   - Click any node in the left **`Node Palette`** to add it to the canvas.
   - Drag wires between matching pin colors (`Signal1D` cyan, `Signal2D` purple, `PeakList` amber).
   - Adjust sliders/dropdowns on any node to re-evaluate downstream nodes in real time; right-click a node to delete it.
4. **Inspect Plots, Peaks & Metrics**:
   - Use the bottom dock tabs: **`1D Spectrum / FID`** (high-to-low ppm axis, zoom/pan with mouse wheel & drag, Re/Im/Magnitude/Baseline overlays), **`2D Heatmap`**, and **`Peaks & SNR Metrics`**.
5. **Export Results & Workflows**:
   - Use **`File → Export Active Spectrum`** to download `.csv`, **JCAMP-DX (`.jdx`)**, **NMRPipe (`.ft1`)**, or `peaks.csv`.
   - Use **`File → Save Workflow (.nmrflow.json)...`** to save your node graph and drag it back onto the app anytime to reload.

---

#### 2. Native Desktop GUI

Requires the [Rust toolchain](https://rustup.rs/). On Linux, install GTK3 dev headers (`sudo apt-get install libgtk-3-dev` on Ubuntu/Debian) for native file dialogs.

```bash
git clone https://github.com/FelixJiashengFu/NMR-READ.git
cd NMR-READ
cargo run --release -p nmr-read-gui
```

> **Tip**: On Desktop, you can also open uncompressed experiment directories directly via **`File → Open Bruker / Varian / Spinsolve Folder...`** or by dragging the folder onto the window.

---

#### 3. Command-Line Interface (`nmr-read-cli`)

```bash
# 1. Inspect any experiment folder or NMR file
cargo run --release -p nmr-read-cli -- inspect path/to/bruker_or_file

# 2. Full 1D processing pipeline (digital filter removal -> FFT -> auto-phase -> arPLS baseline -> Cadzow -> CSV & peaks)
cargo run --release -p nmr-read-cli -- process path/to/bruker_exp \
  --remove-df --fft --auto-phase --baseline-lambda 1000000 --cadzow-rank 16 \
  --out-csv spectrum.csv --out-peaks peaks.csv

# 3. Convert between supported NMR formats (csv, jdx, pipe)
cargo run --release -p nmr-read-cli -- convert path/to/experiment --to jdx --output spectrum.jdx

# 4. Run the built-in 1H Ethanol denoising benchmark (Noisy vs Wavelet vs Cadzow)
cargo run --release -p nmr-read-cli -- benchmark --preset ethanol --points 1024 --noise 0.25
```

---

#### 4. Supported NMR Formats & DSP Features

| Vendor / Format | Files / Extensions | Capabilities |
| :--- | :--- | :--- |
| **Bruker TopSpin** | `acqus`, `acqu2s`, `procs`, `fid`, `ser`, `pdata/1/1r`+`1i` | 1D & 2D raw/processed data, `BYTORDP`/`DTYPA`/`NC_proc`, `GRPDLY` & `DSPFVS`/`DECIM` digital filter removal |
| **JCAMP-DX** | `.dx`, `.jdx` | AFFN + compressed **ASDF** (`SQZ` / `DIF` / `DUP`) decoding, `NTUPLES`/`PAGE` complex spectra, plus exporter |
| **Varian / Agilent** | `procpar` + `fid` | Big-endian 32B file & 28B block headers (`i16`, `i32`, `f32`), 1D & 2D |
| **JEOL Delta** | `.jdf` | Big/little-endian (`JEOL.NMR` / `RMN.LOEJ`), parameter block, 1D & 2D submatrix reordering |
| **NMRPipe** | `.fid`, `.ft1`, `.ft2` | Auto-endian 2048B `fdata` header, 1D & 2D real/complex unappended arrays, plus 1D exporter |
| **Magritek Spinsolve** | `acqu.par` + `data.1d` / `spectrum.1d` / `nmr_fid.dx` | Benchtop 1D binary & JCAMP-DX import |
| **Sparky UCSF** | `.ucsf` | Big-endian 180B global + 128B axis headers + tiled 2D matrix |
| **nmrML** | `.nmrML` | Open XML schema + Base64 (`complex128` / `complex64` / `float64`) arrays |
| **CSV / TSV / JSON** | `.csv`, `.tsv`, `.nmrflow.json` | Tabular import/export, peak table export, serialized node-graph workflows |

- **DSP & Denoising (`nmr-dsp`)**: `rustfft` 1D/2D FFT, half-dwell scaling, DC offset removal, zero filling, Bruker `GRPDLY` phase-ramp removal, EM/GM/Sine-Bell/Hanning/Hamming apodization, Savitzky-Golay smoothing, manual & entropy auto-phasing, $O(N)$ banded pentadiagonal **arPLS** & polynomial baseline subtraction, symmetric-padded 1D **SWT** wavelet denoising (`Haar`/`Db2`/`Db4`/`Db6`/`Sym4`), $O(N \log N)$ FFT-Hankel + Lanczos + `faer` SVD **Cadzow** denoising, parabolic sub-bin peak picking & FWHM, and MAD SNR / RMSE / PSNR / 1D SSIM evaluation.

---

#### 5. Development, Testing & Self-Hosting Wasm

```bash
# Run all workspace unit tests
cargo test --workspace

# Build & serve the WebAssembly app locally with Trunk
rustup target add wasm32-unknown-unknown
cargo install trunk
trunk serve --release
```

---

## 简体中文

### 项目简介

**`NMR-READ`** 是一个完全基于 **纯 Rust** 独立编写的跨平台核磁共振（NMR）数据读取、节点图处理与去噪工作台。项目汲取了多个开源波谱读取与处理生态的设计精华（如 `nmrglue`、`cheminfo/nmr-parser`、`loom`、`NMRTools.jl`、`SpectroChemPy`、`spinsolveproc`、`nmrML` 与 `0ffz/node-nmr`）。

由于所有傅里叶变换、小波、SVD 低秩去噪与厂商二进制解析均使用纯 Rust 实现，**不依赖任何外部 C/Fortran 动态库**（无需 `FFTW3`、`wavelib`、`PROPACK` 或 `OpenBLAS`），一套代码即可直接编译运行在 **macOS（Apple Silicon / Intel）**、**Linux**、**Windows** 以及 **浏览器 WebAssembly** 中。

---

### 使用方法

#### 1. 直接网页打开（零安装 · WebAssembly 在线版）

在浏览器中直接打开：**https://felixjiashengfu.github.io/NMR-READ/**

1. **体验内置预设**：点击顶部菜单栏 **`Presets`**，可一键切换 3 套完整流程：
   - `1H Ethanol: Noisy vs Wavelet vs Cadzow`：对比含噪乙醇 `1H` 谱、Sym4 小波去噪和 Cadzow Hankel SVD 去噪，并实时计算 SNR / RMSE / SSIM 指标。
   - `1H Ethylbenzene: Apod + Auto-Phase + arPLS Baseline`：完整演示切趾、自动相位校正、arPLS 自动基线扣除与自动寻峰。
   - `2D COSY FFT & Heatmap`：合成二维 COSY 矩阵、2D FFT 与交互式 2D 热图可视化。
2. **加载你自己的 NMR 实验数据**：
   - **直接拖拽**：把单个 NMR 文件（`.jdx` / `.dx`、`.jdf`、`.ft1` / `.ft2`、`.ucsf`、`.nmrML`、`.csv`）或**打包成 `.zip` 的 Bruker / Varian / Spinsolve 实验文件夹**直接拖入网页画布；
   - 或点击左上角 **`File → Open File / ZIP...`** 选择文件或压缩包；
   - 加载成功后，画布会自动生成对应的 **`Loaded 1D Dataset`** 或 **`Loaded 2D Dataset`** 数据源节点，你也可以再连入任意后续处理节点。
3. **自由搭建与调节节点管线**：
   - 点击左侧 **`Node Palette`** 中的任意节点即可添加到画布；
   - 拖拽相同颜色的端口进行连线（青色 `Signal1D`、紫色 `Signal2D`、琥珀色 `PeakList`）；
   - 实时拖动节点内的滑块或下拉选项即可自动重算下游结果；右键点击节点可删除。
4. **查看谱图、峰表与信噪比指标**：
   - 底部面板可在 **`1D Spectrum / FID`**（遵循化学位移高场在右、低场在左的 ppm 坐标系，支持滚轮缩放、拖拽平移及 Re / Im / Magnitude / Baseline 叠加）、**`2D Heatmap`** 和 **`Peaks & SNR Metrics`** 之间切换。
5. **导出结果与保存工作流**：
   - 通过 **`File`** 菜单可将当前激活的谱图或峰表下载为 `.csv`、**JCAMP-DX (`.jdx`)**、**NMRPipe (`.ft1`)** 或 `peaks.csv`；
   - 也可以点击 **`Save Workflow (.nmrflow.json)...`** 保存整张节点图配置，下次直接拖回网页即可恢复。

---

#### 2. 本地原生桌面版运行

需预先安装 [Rust 工具链](https://rustup.rs/)（Linux 下若使用系统文件对话框需安装 `libgtk-3-dev`）：

```bash
git clone https://github.com/FelixJiashengFu/NMR-READ.git
cd NMR-READ
cargo run --release -p nmr-read-gui
```

> **提示**：桌面版除了支持单文件与 `.zip`，还支持通过 **`File → Open Bruker / Varian / Spinsolve Folder...`** 或直接把实验目录拖进窗口来读取未压缩文件夹。

---

#### 3. 命令行工具使用（`nmr-read-cli`）

```bash
# 1. 查看任意实验目录或 NMR 文件的格式、维度、核素、谱宽、GRPDLY 与信噪比
cargo run --release -p nmr-read-cli -- inspect path/to/bruker_or_file

# 2. 执行完整 1D 批处理流程（去数字滤波 -> FFT -> 自动相位 -> arPLS 基线 -> Cadzow 去噪 -> 导出 CSV 与峰表）
cargo run --release -p nmr-read-cli -- process path/to/bruker_exp \
  --remove-df --fft --auto-phase --baseline-lambda 1000000 --cadzow-rank 16 \
  --out-csv spectrum.csv --out-peaks peaks.csv

# 3. 格式互转（支持输出 csv / jdx / pipe）
cargo run --release -p nmr-read-cli -- convert path/to/experiment --to jdx --output spectrum.jdx

# 4. 运行内置 1H 乙醇去噪基准测试（对比 Noisy / Wavelet / Cadzow）
cargo run --release -p nmr-read-cli -- benchmark --preset ethanol --points 1024 --noise 0.25
```

---

#### 4. 工作区结构

```text
NMR-READ/
├── crates/
│   ├── nmr-core/       # 1D/2D 信号模型、坐标轴、峰表与 SNR/SSIM 评估指标
│   ├── nmr-io/         # 多厂商读取器（Bruker/Varian/JEOL/Pipe/Spinsolve/Sparky/JCAMP-DX/nmrML/ZIP）
│   ├── nmr-dsp/        # 纯 Rust FFT、Cadzow SVD、SWT 小波、arPLS 基线、自动相位、寻峰与合成谱
│   ├── nmr-read-cli/   # 命令行工具（inspect / process / convert / benchmark）
│   └── nmr-read-gui/   # 桌面端 + 浏览器 WebAssembly 节点图工作台
├── Trunk.toml          # WebAssembly 打包配置
└── index.html          # 浏览器 Wasm 入口页面
```

#### 5. 本地测试与 WebAssembly 本地预览

```bash
# 运行全工作区单元测试
cargo test --workspace

# 本地构建并启动浏览器 Wasm 服务
rustup target add wasm32-unknown-unknown
cargo install trunk
trunk serve --release
```

---

## License

Dual-licensed under either [MIT License](LICENSE-MIT) or [Apache License, Version 2.0](LICENSE-APACHE), at your option.
