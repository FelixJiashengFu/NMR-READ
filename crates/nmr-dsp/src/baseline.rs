use faer::linalg::solvers::SolveLstsq;
use nmr_core::Signal1D;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BaselineMethod {
    /// Asymmetrically Reweighted Penalized Least Squares (Baek et al., Analyst 2015).
    ArPls,
    /// Weighted low-degree polynomial fit rejecting upward peak outliers.
    Polynomial,
}

/// Estimate smooth NMR baseline via **arPLS** in $O(N)$ per iteration using an exact
/// symmetric pentadiagonal solver (no dense $N \times N$ matrices).
///
/// Typical `lambda`: `1e5`..`1e8` (larger = stiffer/smoother baseline).
pub fn arpls_baseline(y: &[f64], lambda: f64, max_iter: usize) -> Vec<f64> {
    let n = y.len();
    if n < 5 {
        return y.to_vec();
    }
    let lam = lambda.max(1.0);
    let mut w = vec![1.0_f64; n];
    let mut z = vec![0.0_f64; n];

    for _ in 0..max_iter.max(1) {
        z = solve_whittaker_d2(&w, y, lam);

        let mut neg_diffs = Vec::new();
        for i in 0..n {
            let d = y[i] - z[i];
            if d < 0.0 {
                neg_diffs.push(d);
            }
        }
        if neg_diffs.is_empty() {
            break;
        }
        let mean_d = neg_diffs.iter().sum::<f64>() / neg_diffs.len() as f64;
        let var_d = neg_diffs
            .iter()
            .map(|&d| (d - mean_d).powi(2))
            .sum::<f64>()
            / neg_diffs.len() as f64;
        let std_d = var_d.sqrt().max(1e-12);

        let mut weight_diff = 0.0_f64;
        let mut weight_norm = 0.0_f64;
        for i in 0..n {
            let d = y[i] - z[i];
            let arg = (2.0 * (d - (-mean_d + 2.0 * std_d)) / std_d).clamp(-60.0, 60.0);
            let w_new = 1.0 / (1.0 + arg.exp());
            weight_diff += (w[i] - w_new).powi(2);
            weight_norm += w[i].powi(2);
            w[i] = w_new.max(1e-8);
        }
        if weight_norm > 0.0 && (weight_diff / weight_norm).sqrt() < 1e-4 {
            break;
        }
    }
    z
}

/// Solve `(diag(w) + lambda * D_2^T D_2) z = diag(w) y` in $O(N)$ using banded Cholesky $L D L^T$
/// for a symmetric positive-definite pentadiagonal matrix.
fn solve_whittaker_d2(w: &[f64], y: &[f64], lam: f64) -> Vec<f64> {
    let n = y.len();
    // D2 rows are [1, -2, 1] for r = 0..n-2.
    // Build main diagonal (d0), first super/sub diagonal (d1), second diagonal (d2).
    let mut a0 = vec![0.0_f64; n];
    let mut a1 = vec![0.0_f64; n]; // a1[i] = A[i, i+1]
    let mut a2 = vec![0.0_f64; n]; // a2[i] = A[i, i+2]
    let mut rhs = vec![0.0_f64; n];

    for i in 0..n {
        a0[i] = w[i];
        rhs[i] = w[i] * y[i];
    }
    for r in 0..(n - 2) {
        a0[r] += lam;
        a0[r + 1] += 4.0 * lam;
        a0[r + 2] += lam;
        a1[r] += -2.0 * lam;
        a1[r + 1] += -2.0 * lam;
        a2[r] += lam;
    }

    // Banded LDL^T factorization (bandwidth 2):
    // L is unit lower triangular with l1[i] = L[i+1, i], l2[i] = L[i+2, i]
    let mut d = vec![0.0_f64; n];
    let mut l1 = vec![0.0_f64; n];
    let mut l2 = vec![0.0_f64; n];

    for i in 0..n {
        let mut diag = a0[i];
        if i >= 1 {
            diag -= l1[i - 1] * l1[i - 1] * d[i - 1];
        }
        if i >= 2 {
            diag -= l2[i - 2] * l2[i - 2] * d[i - 2];
        }
        d[i] = diag.max(1e-18);

        if i + 1 < n {
            let mut v1 = a1[i];
            if i >= 1 {
                v1 -= l2[i - 1] * l1[i - 1] * d[i - 1];
            }
            l1[i] = v1 / d[i];
        }
        if i + 2 < n {
            l2[i] = a2[i] / d[i];
        }
    }

    // Forward solve L u = rhs
    let mut u = vec![0.0_f64; n];
    for i in 0..n {
        let mut s = rhs[i];
        if i >= 1 {
            s -= l1[i - 1] * u[i - 1];
        }
        if i >= 2 {
            s -= l2[i - 2] * u[i - 2];
        }
        u[i] = s;
    }

    // Diagonal scale u <- u / d
    for i in 0..n {
        u[i] /= d[i];
    }

    // Backward solve L^T z = u
    let mut z = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = u[i];
        if i + 1 < n {
            s -= l1[i] * z[i + 1];
        }
        if i + 2 < n {
            s -= l2[i] * z[i + 2];
        }
        z[i] = s;
    }
    z
}

/// Iteratively reweighted polynomial baseline fit (rejecting positive peaks).
pub fn polynomial_baseline(y: &[f64], degree: usize, iters: usize) -> Vec<f64> {
    let n = y.len();
    if n < 4 {
        return y.to_vec();
    }
    let deg = degree.min(8);
    let cols = deg + 1;
    let mut mask = vec![true; n];
    let mut baseline = vec![0.0; n];

    for _ in 0..iters.max(1) {
        let active: Vec<usize> = (0..n).filter(|&i| mask[i]).collect();
        if active.len() <= cols {
            break;
        }
        let mut a = faer::Mat::<f64>::zeros(active.len(), cols);
        let mut b = faer::Mat::<f64>::zeros(active.len(), 1);
        for (row, &idx) in active.iter().enumerate() {
            let x = (idx as f64 / (n - 1).max(1) as f64) * 2.0 - 1.0;
            let mut xp = 1.0;
            for c in 0..cols {
                a[(row, c)] = xp;
                xp *= x;
            }
            b[(row, 0)] = y[idx];
        }
        let sol = a.qr().solve_lstsq(&b);
        for (i, slot) in baseline.iter_mut().enumerate() {
            let x = (i as f64 / (n - 1).max(1) as f64) * 2.0 - 1.0;
            let mut xp = 1.0;
            let mut val = 0.0;
            for c in 0..cols {
                val += sol[(c, 0)] * xp;
                xp *= x;
            }
            *slot = val;
        }

        let mut res: Vec<f64> = (0..n).map(|i| y[i] - baseline[i]).collect();
        res.sort_by(|a, b| a.total_cmp(b));
        let med = res[res.len() / 2];
        let mut abs_dev: Vec<f64> = (0..n).map(|i| (y[i] - baseline[i] - med).abs()).collect();
        abs_dev.sort_by(|a, b| a.total_cmp(b));
        let sigma = (abs_dev[abs_dev.len() / 2] / 0.6745).max(1e-12);

        for i in 0..n {
            mask[i] = (y[i] - baseline[i]) < 2.5 * sigma;
        }
    }
    baseline
}

/// Subtract estimated baseline from the real part of a [`Signal1D`] and attach the baseline curve.
pub fn subtract_baseline(sig: &Signal1D, method: BaselineMethod, param: f64) -> Signal1D {
    if sig.is_empty() {
        return sig.clone();
    }
    let y = sig.real();
    let bl = match method {
        BaselineMethod::ArPls => arpls_baseline(&y, param.max(10.0), 25),
        BaselineMethod::Polynomial => polynomial_baseline(&y, param.round().clamp(1.0, 8.0) as usize, 10),
    };
    let mut out = sig.clone();
    for (c, &b) in out.data.iter_mut().zip(bl.iter()) {
        c.re -= b;
    }
    out.baseline = Some(bl);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arpls_recovers_curved_baseline_under_sharp_peak() {
        let n = 256;
        let true_bl: Vec<f64> = (0..n)
            .map(|i| {
                let x = i as f64 / n as f64;
                0.5 * x * x - 0.3 * x + 1.2
            })
            .collect();
        let mut y = true_bl.clone();
        // Add sharp Lorentzian peak around index 120
        for (i, slot) in y.iter_mut().enumerate() {
            let d = (i as f64 - 120.0) / 3.0;
            *slot += 15.0 / (1.0 + d * d);
        }
        let est = arpls_baseline(&y, 1e5, 30);
        // Baseline error away from peak should be very small
        let err_edge = (est[20] - true_bl[20]).abs();
        let err_under_peak = (est[120] - true_bl[120]).abs();
        assert!(err_edge < 0.05, "edge err={err_edge}");
        assert!(err_under_peak < 0.25, "under-peak err={err_under_peak}");
    }
}
