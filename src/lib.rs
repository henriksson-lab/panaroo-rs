//! Faithful Rust translation of Panaroo.
//!
//! Layout mirrors the Python package one file per module, one function per function.
//! See `PORTING_PLAN.md`. Progress is tracked in `port_order.csv`.
//!
//! Everything under [`support`] has **no** Python counterpart — it replaces third-party
//! libraries (networkx, gffutils, Biopython, scipy, numpy, joblib, intbitset, edlib).
//! Every function outside `support` must trace back to a named Python function.
//!
//! # Lints deliberately relaxed
//!
//! A faithful transcription is not idiomatic Rust, and several clippy lints push away from
//! the Python being translated. These are relaxed crate-wide, with reasons, rather than
//! sprinkled as local `allow`s:
//!
//! - `needless_range_loop` / `manual_memcpy` — the Python indexes by position and the
//!   indices are load-bearing (`hits_trans_dict[member][i]` is zipped against `hits` by
//!   position, `dna[seq_index]` against `protein[seq_index]`). Iterator adaptors would
//!   obscure that.
//! - `map_entry` — `if k not in d: d[k] = ...` is transcribed literally, because several
//!   sites depend on *when* the key is first inserted (`PyDict` is insertion-ordered and
//!   that order reaches the GML).
//! - `collapsible_if` — nested `if`s mirror the Python's nesting, which matters when
//!   reading the two side by side.
//! - `non_snake_case` on locals — a handful keep the Python's name (`exon_count_by_RNA`)
//!   so the transcription is greppable against the original.
//! - `too_many_arguments` — upstream signatures have up to 16 parameters; splitting them
//!   would break the one-function-per-function rule.

#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::map_entry)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::too_many_arguments)]
#![allow(non_snake_case)]

// --- the argparse Namespace, shared by __main__ and set_default_args -------------------
mod args;
pub use args::Args;

// --- 1:1 with panaroo/panaroo/*.py -----------------------------------------------------
pub mod biocode_convert;
pub mod cdhit;
pub mod clean_network;
pub mod find_missing;
pub mod generate_alignments;
pub mod generate_network;
pub mod generate_output;
pub mod get_neighborhood;
pub mod isvalid;
pub mod merge_nodes;
pub mod prokka;
pub mod set_default_args;

// --- infrastructure (no Python counterpart) --------------------------------------------
pub mod support;

/// `panaroo/__init__.py::__version__`
pub const VERSION: &str = "1.8.0";
