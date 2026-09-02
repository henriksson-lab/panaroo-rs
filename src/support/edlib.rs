//! Bindings to the vendored **edlib** C++ library.
//!
//! # Provenance
//!
//! This is the one third-party library the port **links** rather than reimplements. The
//! source is vendored verbatim at `vendor/edlib/` — edlib v1.2.7, MIT, Copyright (c) 2014
//! Martin Šošić — and compiled by `build.rs`. See `vendor/edlib/README.md` and `NOTICE.md`.
//!
//! The bindings below are hand-written (no bindgen, so no libclang at build time). They
//! mirror `vendor/edlib/edlib.h`; if that file is ever updated, re-check the struct
//! layouts here.
//!
//! # Why linked and not reimplemented
//!
//! Panaroo's Python `edlib` package wraps this same C++ code. `cdhit::run_pw`,
//! `find_missing::search_dna` and `cdhit::align_dna_cdhit` feed edit distances straight
//! into threshold comparisons that decide whether two genes merge, so a different aligner
//! — or a different version of this one — changes the pangenome.
//!
//! # Python API being reproduced
//!
//! ```python
//! edlib.align(query, target, mode="HW", task="distance"|"path",
//!             k=..., additionalEqualities=[(a, b), ...])
//! # -> {"editDistance": int, "locations": [(start, end)], "cigar": str|None}
//! ```
//!
//! Two details of the Python wrapper that the raw C API does not give you:
//!   - `locations` is reported even for `task="distance"`, as `(None, end)` — edlib
//!     computes end locations always and start locations only for `LOC`/`PATH`. So
//!     [`AlignResult::locations`] carries an `Option` start. OBSERVED:
//!     `edlib.align("ACGT", "TTACGTTT", mode="HW", task="distance")` gives
//!     `locations == [(None, 5)]`.
//!   - `cigar` is the **extended** format (`=`/`X`/`I`/`D`), not the standard one.
//!     OBSERVED: `task="path"` on the same input gives `'4='`, not `'4M'`. This matters
//!     because `search_dna` does `re.split(r'(\d+)', aln['cigar'])` and indexes the result
//!     positionally — the two formats split into different numbers of tokens.

use std::os::raw::{c_char, c_int, c_uchar};

#[repr(C)]
#[derive(Clone, Copy)]
struct EdlibEqualityPair {
    first: c_char,
    second: c_char,
}

#[repr(C)]
struct EdlibAlignConfig {
    k: c_int,
    mode: c_int,
    task: c_int,
    additional_equalities: *const EdlibEqualityPair,
    additional_equalities_length: c_int,
}

#[repr(C)]
struct EdlibAlignResult {
    status: c_int,
    edit_distance: c_int,
    end_locations: *mut c_int,
    start_locations: *mut c_int,
    num_locations: c_int,
    alignment: *mut c_uchar,
    alignment_length: c_int,
    alphabet_length: c_int,
}

extern "C" {
    fn edlibAlign(
        query: *const c_char,
        query_length: c_int,
        target: *const c_char,
        target_length: c_int,
        config: EdlibAlignConfig,
    ) -> EdlibAlignResult;
    fn edlibFreeAlignResult(result: EdlibAlignResult);
    fn edlibAlignmentToCigar(
        alignment: *const c_uchar,
        alignment_length: c_int,
        cigar_format: c_int,
    ) -> *mut c_char;
}

const EDLIB_MODE_NW: c_int = 0;
const EDLIB_MODE_SHW: c_int = 1;
const EDLIB_MODE_HW: c_int = 2;
const EDLIB_TASK_DISTANCE: c_int = 0;
const EDLIB_TASK_LOC: c_int = 1;
const EDLIB_TASK_PATH: c_int = 2;
#[allow(dead_code)] // the standard format is not what the Python wrapper uses
const EDLIB_CIGAR_STANDARD: c_int = 0;
const EDLIB_CIGAR_EXTENDED: c_int = 1;
const EDLIB_STATUS_OK: c_int = 0;

/// `mode=` argument.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// `"NW"` — global
    Nw,
    /// `"HW"` — infix; the only one Panaroo uses
    Hw,
    /// `"SHW"` — prefix
    Shw,
}

/// `task=` argument.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Task {
    Distance,
    Locations,
    Path,
}

/// The dict returned by `edlib.align`.
#[derive(Debug, Clone, Default)]
pub struct AlignResult {
    /// `-1` when no alignment under `k` was found.
    pub edit_distance: i64,
    /// `(start, end)` pairs; `end` is inclusive, as in edlib. `start` is `None` for
    /// `Task::Distance`, which computes end locations but not start locations — matching
    /// the Python wrapper's `(None, end)`.
    pub locations: Vec<(Option<i64>, i64)>,
    /// Present only for `Task::Path`.
    pub cigar: Option<String>,
}

/// `edlib.align(query, target, mode=..., task=..., k=..., additionalEqualities=...)`
pub fn align(
    query: &str,
    target: &str,
    mode: Mode,
    task: Task,
    k: i64,
    additional_equalities: &[(char, char)],
) -> AlignResult {
    let eqs: Vec<EdlibEqualityPair> = additional_equalities
        .iter()
        .map(|&(a, b)| EdlibEqualityPair {
            first: a as c_char,
            second: b as c_char,
        })
        .collect();

    let cfg = EdlibAlignConfig {
        k: k as c_int,
        mode: match mode {
            Mode::Nw => EDLIB_MODE_NW,
            Mode::Hw => EDLIB_MODE_HW,
            Mode::Shw => EDLIB_MODE_SHW,
        },
        task: match task {
            Task::Distance => EDLIB_TASK_DISTANCE,
            Task::Locations => EDLIB_TASK_LOC,
            Task::Path => EDLIB_TASK_PATH,
        },
        additional_equalities: if eqs.is_empty() {
            std::ptr::null()
        } else {
            eqs.as_ptr()
        },
        additional_equalities_length: eqs.len() as c_int,
    };

    // SAFETY: pointers are valid for the duration of the call; `eqs` outlives it; the
    // result is freed before returning and nothing borrowed from it escapes.
    unsafe {
        let r = edlibAlign(
            query.as_ptr() as *const c_char,
            query.len() as c_int,
            target.as_ptr() as *const c_char,
            target.len() as c_int,
            cfg,
        );

        if r.status != EDLIB_STATUS_OK {
            edlibFreeAlignResult(r);
            panic!("edlib: alignment failed (status != EDLIB_STATUS_OK)");
        }

        let mut out = AlignResult {
            edit_distance: r.edit_distance as i64,
            ..Default::default()
        };

        if r.edit_distance >= 0 && r.num_locations > 0 && !r.end_locations.is_null() {
            let n = r.num_locations as usize;
            let ends = std::slice::from_raw_parts(r.end_locations, n);
            let starts = if r.start_locations.is_null() {
                None
            } else {
                Some(std::slice::from_raw_parts(r.start_locations, n))
            };
            out.locations = (0..n)
                .map(|i| (starts.map(|s| s[i] as i64), ends[i] as i64))
                .collect();
        }

        if task == Task::Path && !r.alignment.is_null() {
            let c = edlibAlignmentToCigar(r.alignment, r.alignment_length, EDLIB_CIGAR_EXTENDED);
            if !c.is_null() {
                out.cigar = Some(std::ffi::CStr::from_ptr(c).to_string_lossy().into_owned());
                libc_free(c as *mut std::ffi::c_void);
            }
        }

        edlibFreeAlignResult(r);
        out
    }
}

extern "C" {
    #[link_name = "free"]
    fn libc_free(p: *mut std::ffi::c_void);
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values from the Python `edlib` package, which wraps this same library.

    #[test]
    fn infix_distance_matches_python_edlib() {
        // python: edlib.align("ACGT", "TTACGTTT", mode="HW", task="distance", k=-1)
        //   -> {'editDistance': 0, ...}
        let r = align("ACGT", "TTACGTTT", Mode::Hw, Task::Distance, -1, &[]);
        assert_eq!(r.edit_distance, 0);
        // the Python wrapper reports (None, 5) here
        assert_eq!(r.locations, [(None, 5)]);
    }

    #[test]
    fn k_bound_reports_minus_one_when_exceeded() {
        // python: edlib.align("AAAA", "TTTT", mode="HW", task="distance", k=1)
        //   -> {'editDistance': -1, ...}
        let r = align("AAAA", "TTTT", Mode::Hw, Task::Distance, 1, &[]);
        assert_eq!(r.edit_distance, -1);
    }

    #[test]
    fn additional_equalities_are_honoured() {
        // Panaroo declares N equal to every base so ambiguity does not cost distance.
        // python: edlib.align("ACGT", "ACNT", mode="HW", task="distance", k=-1,
        //                     additionalEqualities=[('A','N'),('C','N'),('G','N'),('T','N')])
        //   -> editDistance 0   (without the equalities it is 1)
        let eqs = [('A', 'N'), ('C', 'N'), ('G', 'N'), ('T', 'N')];
        assert_eq!(
            align("ACGT", "ACNT", Mode::Hw, Task::Distance, -1, &eqs).edit_distance,
            0
        );
        assert_eq!(
            align("ACGT", "ACNT", Mode::Hw, Task::Distance, -1, &[]).edit_distance,
            1
        );
    }

    #[test]
    fn path_task_yields_locations_and_a_standard_cigar() {
        // python: edlib.align("ACGT", "TTACGTTT", mode="HW", task="path")
        //   -> {'editDistance': 0, 'locations': [(2, 5)], 'cigar': '4='}
        let r = align("ACGT", "TTACGTTT", Mode::Hw, Task::Path, -1, &[]);
        assert_eq!(r.edit_distance, 0);
        assert_eq!(r.locations, [(Some(2), 5)]);
        assert_eq!(r.cigar.as_deref(), Some("4="));
    }
}
