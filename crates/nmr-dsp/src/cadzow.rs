use nmr_core::Signal1D;
use num_complex::Complex64;
use rustfft::FftPlanner;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CadzowConfig {
    /// Target matrix rank (number of complex damped exponentials retained).
    pub rank: usize,
    /// Number of Cadzow alternating projection iterations.
    pub iterations: usize,
    /// Hankel matrix pencil rows `L` as fraction of `N` (typically 0.35..0.5).
    pub pencil_fraction: f64,
}

impl Default for CadzowConfig {
    fn default() -> Self {
        Self {
            rank: 16,
            iterations: 3,
            pencil_fraction: 0.4,
        }
    }
}

/// Pure-Rust FFT-accelerated Hankel matrix operator (`L × K` where `K = N - L + 1`, `H[r, c] = x[r + c]`).
///
/// Computes `H * v` and `H^H * u` in `O(N log N)` without allocating the `L × K` matrix.
struct FftHankelOperator {
    rows: usize,
    cols: usize,
    n_fft: usize,
    fft_x: Vec<Complex64>,
    fft_x_conj_rev: Vec<Complex64>,
    fft_fwd: std::sync::Arc<dyn rustfft::Fft<f64>>,
    fft_inv: std::sync::Arc<dyn rustfft::Fft<f64>>,
}

impl FftHankelOperator {
    fn new(signal: &[Complex64], rows: usize) -> Self {
        let n = signal.len();
        let rows = rows.clamp(1, n);
        let cols = n - rows + 1;
        let n_fft = (n + cols.max(rows)).next_power_of_two();

        let mut planner = FftPlanner::<f64>::new();
        let fft_fwd = planner.plan_fft_forward(n_fft);
        let fft_inv = planner.plan_fft_inverse(n_fft);

        // For H * v: (H v)[r] = sum_{c=0}^{cols-1} x[r + c] v[c]
        // Reverse v -> v_rev[j] = v[cols - 1 - j], then (x * v_rev)[r + cols - 1] = sum_c x[r + c] v[c].
        let mut fft_x = vec![Complex64::new(0.0, 0.0); n_fft];
        fft_x[..n].copy_from_slice(signal);
        fft_fwd.process(&mut fft_x);

        // For H^H * u: (H^H u)[c] = sum_{r=0}^{rows-1} conj(x[r + c]) u[r]
        let mut fft_x_conj_rev = vec![Complex64::new(0.0, 0.0); n_fft];
        for r in 0..n {
            fft_x_conj_rev[r] = signal[r].conj();
        }
        fft_fwd.process(&mut fft_x_conj_rev);

        Self {
            rows,
            cols,
            n_fft,
            fft_x,
            fft_x_conj_rev,
            fft_fwd,
            fft_inv,
        }
    }

    /// Compute `y = H * v` (`v.len() == cols`, returns `rows`).
    fn mul_v(&self, v: &[Complex64]) -> Vec<Complex64> {
        let mut buf = vec![Complex64::new(0.0, 0.0); self.n_fft];
        for c in 0..self.cols {
            buf[self.cols - 1 - c] = v[c];
        }
        self.fft_fwd.process(&mut buf);
        for i in 0..self.n_fft {
            buf[i] *= self.fft_x[i];
        }
        self.fft_inv.process(&mut buf);
        let inv = 1.0 / self.n_fft as f64;
        (0..self.rows)
            .map(|r| buf[r + self.cols - 1] * inv)
            .collect()
    }

    /// Compute `z = H^H * u` (`u.len() == rows`, returns `cols`).
    fn mul_hermitian_u(&self, u: &[Complex64]) -> Vec<Complex64> {
        let mut buf = vec![Complex64::new(0.0, 0.0); self.n_fft];
        for r in 0..self.rows {
            buf[self.rows - 1 - r] = u[r];
        }
        self.fft_fwd.process(&mut buf);
        for i in 0..self.n_fft {
            buf[i] *= self.fft_x_conj_rev[i];
        }
        self.fft_inv.process(&mut buf);
        let inv = 1.0 / self.n_fft as f64;
        (0..self.cols)
            .map(|c| buf[c + self.rows - 1] * inv)
            .collect()
    }
}

fn dot(a: &[Complex64], b: &[Complex64]) -> Complex64 {
    a.iter().zip(b.iter()).map(|(x, y)| x.conj() * y).sum()
}

fn norm(v: &[Complex64]) -> f64 {
    v.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt()
}

/// Golub-Kahan-Lanczos bidiagonalization with full reorthogonalization + `faer` small SVD,
/// approximating the top `rank` singular triplets `(sigma_j, u_j, v_j)` of `H`.
fn top_k_svd_lanczos(
    op: &FftHankelOperator,
    rank: usize,
) -> Vec<(f64, Vec<Complex64>, Vec<Complex64>)> {
    let max_dim = op.rows.min(op.cols);
    let r_target = rank.clamp(1, max_dim);
    let m = (r_target * 2 + 12).min(max_dim);

    // Deterministic starting unit vector in C^cols (avoids any RNG dependency on Wasm)
    let mut v0 = vec![Complex64::new(0.0, 0.0); op.cols];
    for (i, slot) in v0.iter_mut().enumerate() {
        let angle = (i as f64) * 1.61803398875;
        *slot = Complex64::from_polar(1.0 / ((i + 1) as f64).sqrt(), angle);
    }
    let n0 = norm(&v0);
    if n0 == 0.0 {
        return Vec::new();
    }
    for x in &mut v0 {
        *x /= n0;
    }

    let mut v_basis: Vec<Vec<Complex64>> = Vec::with_capacity(m);
    let mut u_basis: Vec<Vec<Complex64>> = Vec::with_capacity(m);
    let mut alphas = Vec::with_capacity(m);
    let mut betas = Vec::with_capacity(m);

    v_basis.push(v0);

    for k in 0..m {
        // u_k = H v_k - beta_{k-1} u_{k-1}
        let mut u = op.mul_v(&v_basis[k]);
        if k > 0 {
            let b_prev = betas[k - 1];
            for (ui, &up) in u.iter_mut().zip(u_basis[k - 1].iter()) {
                *ui -= b_prev * up;
            }
        }
        // Full reorthogonalization of u against previous u_0..u_{k-1}
        for u_prev in &u_basis {
            let proj = dot(u_prev, &u);
            for (ui, &up) in u.iter_mut().zip(u_prev.iter()) {
                *ui -= proj * up;
            }
        }
        let alpha = norm(&u);
        if alpha < 1e-13 {
            break;
        }
        for ui in &mut u {
            *ui /= alpha;
        }
        u_basis.push(u);
        alphas.push(alpha);

        if k + 1 >= m {
            break;
        }

        // v_{k+1} = H^H u_k - alpha_k v_k
        let mut v_next = op.mul_hermitian_u(&u_basis[k]);
        for (vi, &vp) in v_next.iter_mut().zip(v_basis[k].iter()) {
            *vi -= alpha * vp;
        }
        // Full reorthogonalization of v_next against v_0..v_k
        for v_prev in &v_basis {
            let proj = dot(v_prev, &v_next);
            for (vi, &vp) in v_next.iter_mut().zip(v_prev.iter()) {
                *vi -= proj * vp;
            }
        }
        let beta = norm(&v_next);
        if beta < 1e-13 {
            break;
        }
        for vi in &mut v_next {
            *vi /= beta;
        }
        betas.push(beta);
        v_basis.push(v_next);
    }

    let u_steps = u_basis.len();
    let v_steps = v_basis.len();
    if u_steps == 0 || v_steps == 0 {
        return Vec::new();
    }

    // Build real upper-bidiagonal matrix B of size (u_steps × v_steps) and compute SVD via faer
    let mut b_mat = faer::Mat::<f64>::zeros(u_steps, v_steps);
    for i in 0..u_steps {
        if i < v_steps {
            b_mat[(i, i)] = alphas[i];
        }
        if i + 1 < v_steps && i < betas.len() {
            b_mat[(i, i + 1)] = betas[i];
        }
    }
    let svd = b_mat.thin_svd().expect("bidiagonal thin SVD");
    let u_small = svd.U();
    let v_small = svd.V();
    let s_diag = svd.S();

    let keep = r_target.min(u_steps).min(v_steps);
    let mut triplets = Vec::with_capacity(keep);
    for j in 0..keep {
        let sigma = s_diag[j];
        let mut u_full = vec![Complex64::new(0.0, 0.0); op.rows];
        for r_idx in 0..u_steps {
            let coeff = u_small[(r_idx, j)];
            for (dst, &src) in u_full.iter_mut().zip(u_basis[r_idx].iter()) {
                *dst += src * coeff;
            }
        }
        let mut v_full = vec![Complex64::new(0.0, 0.0); op.cols];
        for c_idx in 0..v_steps {
            let coeff = v_small[(c_idx, j)];
            for (dst, &src) in v_full.iter_mut().zip(v_basis[c_idx].iter()) {
                *dst += src * coeff;
            }
        }
        triplets.push((sigma, u_full, v_full));
    }
    triplets
}

/// Reconstruct 1D signal via FFT-accelerated antidiagonal averaging of `sum_j sigma_j u_j v_j^H`.
fn reconstruct_antidiagonal_fft(
    n: usize,
    rows: usize,
    cols: usize,
    triplets: &[(f64, Vec<Complex64>, Vec<Complex64>)],
) -> Vec<Complex64> {
    let n_fft = (rows + cols).next_power_of_two();
    let mut planner = FftPlanner::<f64>::new();
    let fft_fwd = planner.plan_fft_forward(n_fft);
    let fft_inv = planner.plan_fft_inverse(n_fft);

    let mut acc_freq = vec![Complex64::new(0.0, 0.0); n_fft];
    let mut u_buf = vec![Complex64::new(0.0, 0.0); n_fft];
    let mut v_buf = vec![Complex64::new(0.0, 0.0); n_fft];

    for (sigma, u, v) in triplets {
        u_buf.fill(Complex64::new(0.0, 0.0));
        v_buf.fill(Complex64::new(0.0, 0.0));
        u_buf[..rows].copy_from_slice(u);
        for c in 0..cols {
            v_buf[c] = v[c].conj();
        }
        fft_fwd.process(&mut u_buf);
        fft_fwd.process(&mut v_buf);
        for i in 0..n_fft {
            acc_freq[i] += (*sigma) * u_buf[i] * v_buf[i];
        }
    }

    fft_inv.process(&mut acc_freq);
    let inv_fft = 1.0 / n_fft as f64;
    let min_rc = rows.min(cols);

    (0..n)
        .map(|k| {
            // Number of elements along antidiagonal r + c = k
            let count = (k + 1).min(min_rc).min(n - k).max(1) as f64;
            acc_freq[k] * (inv_fft / count)
        })
        .collect()
}

/// Apply Cadzow Hankel low-rank iterative denoising to a 1D time-domain FID (or frequency signal).
pub fn cadzow_denoise(sig: &Signal1D, cfg: &CadzowConfig) -> Signal1D {
    let n = sig.len();
    if n < 12 {
        return sig.clone();
    }
    let rows = ((n as f64) * cfg.pencil_fraction.clamp(0.15, 0.5))
        .round()
        .clamp(4.0, (n / 2) as f64) as usize;
    let cols = n - rows + 1;

    let mut current = sig.data.clone();
    for _ in 0..cfg.iterations.max(1) {
        let op = FftHankelOperator::new(&current, rows);
        let triplets = top_k_svd_lanczos(&op, cfg.rank);
        if triplets.is_empty() {
            break;
        }
        current = reconstruct_antidiagonal_fft(n, rows, cols, &triplets);
    }

    let mut out = sig.clone();
    out.data = current;
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use nmr_core::{AcquisitionParams, Domain};

    #[test]
    fn fft_hankel_operator_matches_dense_multiplication() {
        let sig: Vec<Complex64> = (0..16)
            .map(|i| Complex64::new((i as f64 * 0.4).sin(), (i as f64 * 0.25).cos()))
            .collect();
        let rows = 6;
        let cols = sig.len() - rows + 1;
        let op = FftHankelOperator::new(&sig, rows);

        let v: Vec<Complex64> = (0..cols)
            .map(|c| Complex64::new(c as f64 + 1.0, -(c as f64) * 0.3))
            .collect();
        let y_fft = op.mul_v(&v);
        for r in 0..rows {
            let mut y_ref = Complex64::new(0.0, 0.0);
            for c in 0..cols {
                y_ref += sig[r + c] * v[c];
            }
            assert!((y_fft[r] - y_ref).norm() < 1e-10);
        }

        let u: Vec<Complex64> = (0..rows)
            .map(|r| Complex64::new(0.5 * (r as f64), 1.0))
            .collect();
        let z_fft = op.mul_hermitian_u(&u);
        for c in 0..cols {
            let mut z_ref = Complex64::new(0.0, 0.0);
            for r in 0..rows {
                z_ref += sig[r + c].conj() * u[r];
            }
            assert!((z_fft[c] - z_ref).norm() < 1e-10);
        }
    }

    #[test]
    fn cadzow_reconstructs_noiseless_rank2_fid() {
        let n = 128;
        let fid: Vec<Complex64> = (0..n)
            .map(|i| {
                let t = i as f64 * 0.01;
                Complex64::from_polar((-1.5 * t).exp(), 12.0 * t)
                    + Complex64::from_polar(0.7 * (-2.0 * t).exp(), -25.0 * t)
            })
            .collect();
        let sig = Signal1D::new(fid.clone(), Domain::Time, AcquisitionParams::default());
        let denoised = cadzow_denoise(
            &sig,
            &CadzowConfig {
                rank: 2,
                iterations: 2,
                pencil_fraction: 0.4,
            },
        );
        let max_err = fid
            .iter()
            .zip(denoised.data.iter())
            .map(|(a, b)| (*a - *b).norm())
            .fold(0.0_f64, f64::max);
        assert!(max_err < 1e-6, "max_err={max_err}");
    }
}
