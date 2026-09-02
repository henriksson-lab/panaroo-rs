//! Translation of `panaroo/panaroo/set_default_args.py`.

use crate::Args;

/// `set_default_args.py::set_default_args`
///
/// Fills in the options left `None` by `get_options`, based on `--mode`. Mutates in place
/// and returns, as the Python does.
///
/// Watch the types: `min_trailing_support`, `min_edge_support_sv` and
/// `edge_support_threshold` get `max(2, math.ceil(0.05 * n))` — an `int` — in strict and
/// moderate mode, but `edge_support_threshold` gets the **float** `0.0` in sensitive mode.
/// That difference is visible if the value is ever printed.
pub fn set_default_args(args: &mut Args) {
    use crate::support::pyfmt::py_ceil;

    let n_samples = args.input_files.len() as f64;

    // The three modes share `id`, `family_threshold` and `len_dif_percent`; they differ in
    // the trailing-end and edge-support defaults.
    if args.id.is_none() {
        args.id = Some(0.98);
    }
    if args.family_threshold.is_none() {
        args.family_threshold = Some(0.7);
    }
    if args.len_dif_percent.is_none() {
        args.len_dif_percent = Some(0.98);
    }

    match args.mode.as_str() {
        "strict" => {
            if args.min_trailing_support.is_none() {
                args.min_trailing_support = Some(2.max(py_ceil(0.05 * n_samples)));
            }
            if args.trailing_recursive.is_none() {
                args.trailing_recursive = Some(99999999);
            }
            if args.min_edge_support_sv.is_none() {
                args.min_edge_support_sv = Some(2.max(py_ceil(0.01 * n_samples)));
            }
            if args.remove_by_consensus.is_none() {
                args.remove_by_consensus = Some(true);
            }
            if args.edge_support_threshold.is_none() {
                args.edge_support_threshold = Some(2.max(py_ceil(0.01 * n_samples)) as f64);
            }
        }
        "moderate" => {
            if args.min_trailing_support.is_none() {
                args.min_trailing_support = Some(2.max(py_ceil(0.01 * n_samples)));
            }
            if args.trailing_recursive.is_none() {
                args.trailing_recursive = Some(99999999);
            }
            if args.min_edge_support_sv.is_none() {
                args.min_edge_support_sv = Some(2.max(py_ceil(0.01 * n_samples)));
            }
            if args.remove_by_consensus.is_none() {
                args.remove_by_consensus = Some(false);
            }
            if args.edge_support_threshold.is_none() {
                args.edge_support_threshold = Some(2.max(py_ceil(0.01 * n_samples)) as f64);
            }
        }
        _ => {
            if args.min_trailing_support.is_none() {
                args.min_trailing_support = Some(2);
            }
            if args.trailing_recursive.is_none() {
                args.trailing_recursive = Some(0);
            }
            if args.min_edge_support_sv.is_none() {
                args.min_edge_support_sv = Some(2);
            }
            if args.remove_by_consensus.is_none() {
                args.remove_by_consensus = Some(false);
            }
            if args.edge_support_threshold.is_none() {
                // Note: the float 0.0 here, where the other two modes set an int. The type
                // difference is visible if the value is ever printed.
                args.edge_support_threshold = Some(0.0);
            }
        }
    }
}
