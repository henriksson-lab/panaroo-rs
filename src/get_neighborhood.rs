//! Translation of `panaroo/panaroo/get_neighborhood.py`.
//!
//! Only `conv_list` is reachable from `__main__::main` (through the `generate_output` star
//! import chain). The rest of the module belongs to `panaroo-gene-neighbourhood`, phase 6.

/// `get_neighborhood.py::conv_list`
///
/// Byte-identical to [`crate::isvalid::conv_list`] upstream. Both are kept per
/// PORTING_PLAN.md §5 — rule 2 says if the Python duplicated it, the Rust duplicates it.
pub fn conv_list(maybe_list: Vec<String>) -> Vec<String> {
    maybe_list
}

/// The scalar arm of `get_neighborhood.py::conv_list`; see
/// [`crate::isvalid::conv_list_scalar`].
pub fn conv_list_scalar(maybe_list: String) -> Vec<String> {
    vec![maybe_list]
}
