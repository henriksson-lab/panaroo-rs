//! `scipy.sparse.csr_matrix` and `scipy.sparse.csgraph.connected_components`.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [SciPy](https://scipy.org/) 1.14.1 (`scipy.sparse`, `scipy.sparse.csgraph`),
//! BSD 3-Clause, Copyright (c) 2001-2002 Enthought, Inc., 2003-2024 SciPy Developers. No
//! SciPy source was copied; the rules below were established by running the library. See
//! `NOTICE.md`.
//!
//! # Why the *label numbering* matters
//!
//! `distances_bwtn_centroids` is an `ncentroids x ncentroids` boolean matrix built in
//! `cdhit::pwdist_edlib` and sliced per node in `clean_network::single_linkage`. That
//! function feeds the returned labels to `np.unique(labels)`, which fixes cluster order,
//! which fixes merge order, which fixes the new node IDs written to both GML files. So
//! reproducing a *correct partition* is not enough — the numbering has to match.
//!
//! It does, and simply: SciPy scans vertices `0..n-1`, and each vertex not yet labelled
//! starts a new component taking the next label in sequence. OBSERVED on SciPy 1.14.1:
//!
//! | graph (n=5) | labels |
//! |---|---|
//! | edges 0-2, 1-4 | `[0, 1, 0, 2, 1]` |
//! | edges 4-0, 3-1, 1-2 | `[0, 1, 1, 1, 0]` |
//! | no edges (n=4) | `[0, 1, 2, 3]` |
//!
//! # Constructor details
//!
//! `csr_matrix((data, (row, col)), shape=...)` **sums duplicate coordinates** and sorts
//! column indices within each row. OBSERVED: two entries at `(0,1)` with value 1 give a
//! single stored value of 2.

/// `scipy.sparse.csr_matrix`
#[derive(Debug, Clone, Default)]
pub struct CsrMatrix {
    pub shape: (usize, usize),
    /// Row pointers, `shape.0 + 1` entries.
    pub indptr: Vec<usize>,
    /// Column indices, sorted within each row.
    pub indices: Vec<usize>,
    pub data: Vec<i64>,
}

impl CsrMatrix {
    /// `csr_matrix((data, (row_ind, col_ind)), shape=shape)`
    ///
    /// Duplicate coordinates are summed; column indices are sorted within each row; explicit
    /// zeros produced by summing are kept, as SciPy keeps them until `eliminate_zeros()`.
    pub fn from_coo(
        data: Vec<i64>,
        row_ind: Vec<usize>,
        col_ind: Vec<usize>,
        shape: (usize, usize),
    ) -> Self {
        assert_eq!(data.len(), row_ind.len());
        assert_eq!(data.len(), col_ind.len());

        let mut rows: Vec<Vec<(usize, i64)>> = vec![Vec::new(); shape.0];
        for i in 0..data.len() {
            rows[row_ind[i]].push((col_ind[i], data[i]));
        }

        let mut indptr = Vec::with_capacity(shape.0 + 1);
        let mut indices = Vec::new();
        let mut out_data = Vec::new();
        indptr.push(0);
        for r in rows.iter_mut() {
            r.sort_by_key(|(c, _)| *c);
            let mut j = 0;
            while j < r.len() {
                let c = r[j].0;
                let mut acc = 0i64;
                while j < r.len() && r[j].0 == c {
                    acc += r[j].1;
                    j += 1;
                }
                indices.push(c);
                out_data.push(acc);
            }
            indptr.push(indices.len());
        }

        CsrMatrix {
            shape,
            indptr,
            indices,
            data: out_data,
        }
    }

    /// `m.nonzero()` — `(rows, cols)` in row-major order, columns ascending within a row.
    ///
    /// SciPy omits explicitly-stored zeros here.
    pub fn nonzero(&self) -> (Vec<usize>, Vec<usize>) {
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        for r in 0..self.shape.0 {
            for k in self.indptr[r]..self.indptr[r + 1] {
                if self.data[k] != 0 {
                    rows.push(r);
                    cols.push(self.indices[k]);
                }
            }
        }
        (rows, cols)
    }

    /// `m[index][:, index]` — the double fancy-index in `single_linkage`.
    ///
    /// Row `i` of the result is row `index[i]` of the original, restricted to the columns
    /// in `index` and renumbered to their positions there. `index` may repeat or reorder.
    ///
    /// Kept as one operation because that is what the Python does; PORTING_PLAN.md §9 item
    /// 6 explains why it is slow and why it stays slow until parity is signed off.
    pub fn submatrix(&self, index: &[usize]) -> CsrMatrix {
        // Column original -> list of positions in `index` (a column may appear twice).
        let mut col_pos: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (p, &c) in index.iter().enumerate() {
            col_pos.entry(c).or_default().push(p);
        }

        let n = index.len();
        let mut indptr = Vec::with_capacity(n + 1);
        let mut indices = Vec::new();
        let mut data = Vec::new();
        indptr.push(0);
        for &r in index {
            let mut row: Vec<(usize, i64)> = Vec::new();
            for k in self.indptr[r]..self.indptr[r + 1] {
                if let Some(ps) = col_pos.get(&self.indices[k]) {
                    for &p in ps {
                        row.push((p, self.data[k]));
                    }
                }
            }
            row.sort_by_key(|(c, _)| *c);
            for (c, v) in row {
                indices.push(c);
                data.push(v);
            }
            indptr.push(indices.len());
        }

        CsrMatrix {
            shape: (n, n),
            indptr,
            indices,
            data,
        }
    }

    /// The transpose half of [`Self::undirected_neighbours`], precomputed once.
    ///
    /// `rev[v]` is exactly the sequence the old inner `for r in 0..shape.0` scan produced
    /// for `v`: every row `r` storing a nonzero at column `v`, in ascending `r`, once per
    /// stored entry. Built by the same ascending scan with the same `data[k] != 0` filter,
    /// so it is element-for-element identical -- this hoists a loop-invariant scan out of
    /// the flood fill, it does not change traversal.
    fn transpose_lists(&self) -> Vec<Vec<usize>> {
        let mut rev: Vec<Vec<usize>> = vec![Vec::new(); self.shape.0.max(self.shape.1)];
        for r in 0..self.shape.0 {
            for k in self.indptr[r]..self.indptr[r + 1] {
                if self.data[k] != 0 {
                    rev[self.indices[k]].push(r);
                }
            }
        }
        rev
    }

    /// Neighbours of `v`, treating the matrix as undirected (`directed=False`).
    ///
    /// `rev` comes from [`Self::transpose_lists`] on this same matrix.
    fn undirected_neighbours(&self, v: usize, rev: &[Vec<usize>]) -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        for k in self.indptr[v]..self.indptr[v + 1] {
            if self.data[k] != 0 {
                out.push(self.indices[k]);
            }
        }
        // the transpose direction, precomputed
        out.extend_from_slice(&rev[v]);
        out
    }
}

/// `scipy.sparse.csgraph.connected_components(csgraph=m, directed=directed,
/// return_labels=True)` -> `(n_components, labels)`.
///
/// Labels are assigned by scanning vertices `0..n-1`: the first unlabelled vertex starts
/// component 0, the next unlabelled one component 1, and so on. See the module docs for the
/// observed cases this reproduces.
pub fn connected_components(m: &CsrMatrix, _directed: bool) -> (usize, Vec<i64>) {
    let n = m.shape.0;
    let rev = m.transpose_lists();
    let mut labels = vec![-1i64; n];
    let mut ncomp = 0i64;
    for v in 0..n {
        if labels[v] != -1 {
            continue;
        }
        // flood fill
        let mut stack = vec![v];
        labels[v] = ncomp;
        while let Some(u) = stack.pop() {
            for w in m.undirected_neighbours(u, &rev) {
                if labels[w] == -1 {
                    labels[w] = ncomp;
                    stack.push(w);
                }
            }
        }
        ncomp += 1;
    }
    (ncomp as usize, labels)
}

/// `scipy.sparse.csgraph.shortest_path` — imported by `clean_network` but never called on
/// the phase-1 path (the live call is `nx.shortest_path_length`). Stubbed for completeness.
pub fn shortest_path(_m: &CsrMatrix) -> Vec<Vec<f64>> {
    panic!("noimpl: support::sparse::shortest_path (unused on the phase-1 path)")
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every expectation is a value produced by SciPy 1.14.1.

    fn coo(data: Vec<i64>, r: Vec<usize>, c: Vec<usize>, n: usize) -> CsrMatrix {
        CsrMatrix::from_coo(data, r, c, (n, n))
    }

    #[test]
    fn label_numbering_matches_scipy() {
        // scipy: labels [0, 1, 0, 2, 1], n_components 3
        let m = coo(vec![1, 1], vec![0, 1], vec![2, 4], 5);
        assert_eq!(connected_components(&m, false), (3, vec![0, 1, 0, 2, 1]));

        // scipy: labels [0, 1, 2, 0, 3], n_components 4
        let m = coo(vec![1], vec![3], vec![0], 5);
        assert_eq!(connected_components(&m, false), (4, vec![0, 1, 2, 0, 3]));

        // scipy: labels [0, 1, 1, 1, 0], n_components 2
        let m = coo(vec![1, 1, 1], vec![4, 3, 1], vec![0, 1, 2], 5);
        assert_eq!(connected_components(&m, false), (2, vec![0, 1, 1, 1, 0]));

        // scipy: labels [0, 1, 2, 3], n_components 4
        let m = coo(vec![], vec![], vec![], 4);
        assert_eq!(connected_components(&m, false), (4, vec![0, 1, 2, 3]));
    }

    #[test]
    fn constructor_sums_duplicate_coordinates() {
        // scipy: data [2], indices [1], indptr [0 1 1 1]
        let m = coo(vec![1, 1], vec![0, 0], vec![1, 1], 3);
        assert_eq!(m.data, [2]);
        assert_eq!(m.indices, [1]);
        assert_eq!(m.indptr, [0, 1, 1, 1]);
    }

    #[test]
    fn nonzero_is_row_major() {
        // scipy: rows [0, 1, 2], cols [2, 1, 0]
        let m = coo(vec![1, 1, 1], vec![2, 0, 1], vec![0, 2, 1], 3);
        assert_eq!(m.nonzero(), (vec![0, 1, 2], vec![2, 1, 0]));
    }

    #[test]
    fn submatrix_matches_scipy_double_fancy_index() {
        // scipy: m = csr((( [1,1], ([0,1],[3,2]) )), shape=(4,4)); idx=[3,1,2]
        //        m[idx][:, idx].toarray() == [[0,0,0],[0,0,1],[0,0,0]]
        //        connected_components -> (2, [0, 1, 1])
        let m = coo(vec![1, 1], vec![0, 1], vec![3, 2], 4);
        let sub = m.submatrix(&[3, 1, 2]);
        assert_eq!(sub.shape, (3, 3));
        let mut dense = vec![vec![0i64; 3]; 3];
        for r in 0..3 {
            for k in sub.indptr[r]..sub.indptr[r + 1] {
                dense[r][sub.indices[k]] = sub.data[k];
            }
        }
        assert_eq!(dense, vec![vec![0, 0, 0], vec![0, 0, 1], vec![0, 0, 0]]);
        assert_eq!(connected_components(&sub, false), (2, vec![0, 1, 1]));
    }
}
