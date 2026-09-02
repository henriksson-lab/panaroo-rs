//! The `argparse.Namespace` produced by `__main__.py::get_options`.
//!
//! Not a translation of a Python *function* — it is the shape of the object that
//! `get_options` returns and `set_default_args` mutates. It lives in the library rather
//! than in `main.rs` because both modules need it.
//!
//! Fields whose `Option` is `None` after `get_options` are filled in by
//! [`crate::set_default_args::set_default_args`] based on `mode`.

#[derive(Debug, Clone)]
pub struct Args {
    // --- Input/output ---
    /// `-i/--input` (`dest="input_files"`, `nargs='+'`, required)
    pub input_files: Vec<String>,
    /// `-o/--out_dir` (required). `main` normalises this to a trailing separator.
    pub output_dir: String,
    /// `--clean-mode` (`dest="mode"`, required): `strict` | `moderate` | `sensitive`
    pub mode: String,
    /// `--remove-invalid-genes` (`dest="filter_invalid"`, default `False`)
    pub filter_invalid: bool,

    // --- Matching ---
    /// `-c/--threshold` (`dest="id"`); default filled by mode
    pub id: Option<f64>,
    /// `-f/--family_threshold`; default filled by mode
    pub family_threshold: Option<f64>,
    /// `--len_dif_percent`; default filled by mode
    pub len_dif_percent: Option<f64>,
    /// `--family_len_dif_percent`
    ///
    /// `type=float` with **no** `default=`, so it is `None` unless the user supplies it,
    /// and `set_default_args` never fills it in. The help text's "default=0.0" is wrong.
    /// The value reaches cd-hit's `-s` flag, which therefore really does receive the
    /// literal string `None` on a default run.
    pub family_len_dif_percent: Option<f64>,
    /// `--merge_paralogs` (default `False`)
    pub merge_paralogs: bool,

    // --- Refind ---
    /// `--search_radius` (default `5000`)
    pub search_radius: i64,
    /// `--refind_prop_match` (default `0.2`)
    pub refind_prop_match: f64,
    /// `--refind-mode` (default `"default"`): `default` | `strict` | `off`
    pub refind_mode: String,

    // --- Graph correction; defaults filled by mode ---
    /// `--min_trailing_support`
    pub min_trailing_support: Option<i64>,
    /// `--trailing_recursive`
    pub trailing_recursive: Option<i64>,
    /// `--edge_support_threshold`
    ///
    /// Note the type is inconsistent upstream: an `int` from `max(2, ceil(...))` in strict
    /// and moderate mode, but the float `0.0` in sensitive mode.
    pub edge_support_threshold: Option<f64>,
    /// `--length_outlier_support_proportion`
    ///
    /// // UPSTREAM BUG: the help text says `default=0.01`; the actual default is `0.1`.
    pub length_outlier_support_proportion: f64,
    /// `--remove_by_consensus` (`type=ast.literal_eval`, choices `[True, False]`);
    /// default filled by mode
    pub remove_by_consensus: Option<bool>,
    /// `--cycle_threshold_min` (default `5`)
    pub cycle_threshold_min: i64,
    /// `--min_edge_support_sv`; default filled by mode
    pub min_edge_support_sv: Option<i64>,
    /// `--all_seq_in_graph` (default `False`)
    pub all_seq_in_graph: bool,
    /// `--no_clean_edges` (`dest="clean_edges"`, `store_false`, default `True`)
    pub clean_edges: bool,

    // --- Alignment ---
    /// `-a/--alignment` (`dest="aln"`, default `None`): `core` | `pan`
    pub aln: Option<String>,
    /// `--aligner` (`dest="alr"`, default `"mafft"`)
    pub alr: String,
    /// `--codons` (default `False`)
    pub codons: bool,
    /// `--strict_codons` (default `False`)
    pub strict_codons: bool,
    /// `--core_threshold` (`dest="core"`, default `0.95`)
    pub core: f64,
    /// `--core_subset` (`dest="subset"`, default `None` = all)
    pub subset: Option<i64>,
    /// `--core_entropy_filter` (`dest="hc_threshold"`, default `None`)
    pub hc_threshold: Option<f64>,

    // --- Misc ---
    /// `-t/--threads` (`dest="n_cpu"`, default `1`)
    ///
    /// Signed, because a negative value is meaningful: `joblib.Parallel(n_jobs=-1)` means
    /// "all cores", `-2` means "all but one". The raw value also goes into the cd-hit
    /// command string verbatim.
    pub n_cpu: i64,
    /// `--codon-table` (`dest="table"`, default `11`)
    pub table: i64,
    /// `--quiet` (`dest="verbose"`, `store_false`, default `True`)
    pub verbose: bool,
}
