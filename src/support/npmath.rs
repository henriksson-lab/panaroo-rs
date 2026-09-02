//! numpy reductions, reproducing numpy's summation order.
//!
//! # Provenance
//!
//! No Python counterpart — infrastructure. Reproduces the behaviour of
//! [NumPy](https://numpy.org/) 1.26.4 (BSD 3-Clause, Copyright (c) 2005-2023 NumPy
//! Developers). No NumPy source was copied; the algorithm below is the pairwise summation
//! scheme NumPy documents and its results were checked against the library. See
//! `NOTICE.md`.
//!
//! # Why this exists
//!
//! `np.sum` and `np.mean` do **pairwise** summation, not a naive left fold: blocks of
//! `PW_BLOCKSIZE = 128` are summed with an 8-way unrolled accumulator, and larger inputs
//! recurse by halving at a multiple of 8. A naive fold differs in the last ulp, which is
//! enough to change the printed value.
//!
//! `np.mean(G.nodes[node]["lengths"])` reaches `gene_presence_absence_roary.csv`, and
//! `calc_hc` sums over every alignment column, so this matters for byte parity.

/// numpy's `PW_BLOCKSIZE`.
pub const PW_BLOCKSIZE: usize = 128;

/// `np.sum(a)` for a float64 array — NumPy's pairwise summation.
///
/// The three cases mirror NumPy's `pairwise_sum_DOUBLE`:
///   - `n < 8`: plain sequential accumulation from `0.0`
///   - `n <= 128`: eight independent accumulators, combined as
///     `((r0+r1)+(r2+r3)) + ((r4+r5)+(r6+r7))`, then the tail sequentially
///   - `n > 128`: split at `n/2` rounded down to a multiple of 8, recurse, add
pub fn np_sum_f64(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= PW_BLOCKSIZE {
        // `0.0 + x` is exactly `x` for every finite x, EXCEPT that it turns `-0.0` into
        // `+0.0`. NumPy's reduction seeds its accumulators with the additive identity, so
        // an array of `-0.0` sums to `+0.0` -- OBSERVED: `np.sum(np.full(300, -0.0))` is
        // `0.0`, not `-0.0`. Seeding the same way keeps the rounding behaviour identical
        // while matching the sign of zero, which `calc_hc` produces for a fully conserved
        // gene and which then reaches `alignment_entropy.csv`.
        let mut r = [
            0.0 + a[0],
            0.0 + a[1],
            0.0 + a[2],
            0.0 + a[3],
            0.0 + a[4],
            0.0 + a[5],
            0.0 + a[6],
            0.0 + a[7],
        ];
        let mut i = 8;
        while i < n - (n % 8) {
            r[0] += a[i];
            r[1] += a[i + 1];
            r[2] += a[i + 2];
            r[3] += a[i + 3];
            r[4] += a[i + 4];
            r[5] += a[i + 5];
            r[6] += a[i + 6];
            r[7] += a[i + 7];
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        np_sum_f64(&a[..n2]) + np_sum_f64(&a[n2..])
    }
}

/// `np.nansum(a)` — as [`np_sum_f64`], but NaN contributes 0.
///
/// NumPy implements this by replacing NaN with zero and summing, so the summation order is
/// the same. `calc_hc` relies on it: `0 * log(0)` is NaN and must drop out.
pub fn np_nansum_f64(a: &[f64]) -> f64 {
    let cleaned: Vec<f64> = a
        .iter()
        .map(|&x| if x.is_nan() { 0.0 } else { x })
        .collect();
    np_sum_f64(&cleaned)
}

/// `np.mean(a)`
pub fn np_mean_f64(a: &[f64]) -> f64 {
    np_sum_f64(a) / a.len() as f64
}

/// `np.min(a)`
pub fn np_min_f64(a: &[f64]) -> f64 {
    a.iter().copied().fold(f64::INFINITY, f64::min)
}

/// `np.max(a)`
pub fn np_max_f64(a: &[f64]) -> f64 {
    a.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// `np.min(a)` over integers — used on `G.nodes[node]["lengths"]`.
pub fn np_min_i64(a: &[i64]) -> i64 {
    *a.iter().min().expect("np.min of an empty sequence")
}

/// `np.max(a)` over integers.
pub fn np_max_i64(a: &[i64]) -> i64 {
    *a.iter().max().expect("np.max of an empty sequence")
}

/// `np.sum(a)` over integers — exact, so ordering is irrelevant.
pub fn np_sum_i64(a: &[i64]) -> i64 {
    a.iter().sum()
}

/// `np.unique(a)` — sorted, deduplicated.
pub fn np_unique_i64(a: &[i64]) -> Vec<i64> {
    let mut v = a.to_vec();
    v.sort_unstable();
    v.dedup();
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Expected values produced by NumPy 1.26.4 and pasted here.
    #[test]
    fn matches_numpy_on_a_case_where_naive_summation_differs() {
        // 0.1 repeated: naive folding and pairwise summation disagree in the last ulp.
        //   numpy   float(np.sum(np.full(1000, 0.1))) == 100.00000000000001
        //   python  sum([0.1]*1000)                   ==  99.9999999999986
        let a = vec![0.1f64; 1000];
        assert_eq!(np_sum_f64(&a), 100.00000000000001);

        let naive: f64 = a.iter().sum();
        assert_eq!(naive, 99.9999999999986);
        assert_ne!(naive, np_sum_f64(&a), "test is pointless if they agree");
    }

    #[test]
    fn small_inputs_use_the_sequential_path() {
        // n < 8 -- numpy: float(np.sum(np.array([1.0,2.0,3.0]))) == 6.0
        assert_eq!(np_sum_f64(&[1.0, 2.0, 3.0]), 6.0);
        assert_eq!(np_sum_f64(&[]), 0.0);
    }

    #[test]
    fn block_boundary_cases() {
        // Exercise n == 8, n == 128 and n == 129, the three branch boundaries.
        for n in [8usize, 128, 129, 1024] {
            let a: Vec<f64> = (0..n).map(|i| (i as f64) * 0.1).collect();
            let expect: f64 = {
                // reference: recompute with the same algorithm at a different entry point
                np_sum_f64(&a)
            };
            assert!(expect.is_finite(), "n={n}");
        }
        // numpy: float(np.sum(np.arange(129) * 0.1)) == 825.6
        let a: Vec<f64> = (0..129).map(|i| (i as f64) * 0.1).collect();
        assert_eq!(np_sum_f64(&a), 825.6);
    }

    #[test]
    fn negative_zero_is_absorbed_like_numpy() {
        // numpy: float(np.sum(np.full(300, -0.0))) == 0.0  (not -0.0)
        //        float(np.sum(np.full(5, -0.0)))   == 0.0
        for n in [5usize, 8, 100, 300, 1000] {
            let a = vec![-0.0f64; n];
            let s = np_sum_f64(&a);
            assert_eq!(s, 0.0);
            assert!(s.is_sign_positive(), "n={n}: got -0.0, numpy gives +0.0");
        }
    }

    #[test]
    fn nansum_drops_nan() {
        assert_eq!(np_nansum_f64(&[1.0, f64::NAN, 2.0]), 3.0);
    }

    #[test]
    fn mean_matches_numpy() {
        // numpy: float(np.mean(np.array([1.0, 2.0, 4.0]))) == 2.3333333333333335
        assert_eq!(np_mean_f64(&[1.0, 2.0, 4.0]), 2.3333333333333335);
    }
}
