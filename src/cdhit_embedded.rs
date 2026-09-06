//! Run cd-hit in-process, via the [`cdhit-rs`] translation, instead of spawning the C++ binary.
//!
//! Enabled by the `cdhit-embedded` Cargo feature, which is off by default.
//!
//! [`cdhit-rs`]: https://github.com/henriksson-lab/cdhit-rs
//!
//! # What this stage does and does not change
//!
//! This is the **first** step of the integration described in `cdhit_backend.rs`, and it is
//! deliberately the conservative one: it removes the *process*, not the *files*. The same
//! input FASTA is written, the same outputs are produced by the same code, at the same
//! paths. What disappears is `fork`/`exec` and the shell.
//!
//! Doing it this way makes the change testable in isolation. If the parity suite still
//! reports 13/13 byte-identical with this enabled, that is evidence about the translation
//! itself — `cdhit-rs` clusters the way the C++ does — uncontaminated by any change to how
//! data reaches it. Going in-memory first would have conflated the two, and a failure would
//! not say which half was wrong.
//!
//! The in-memory path (no temporary FASTA, representatives returned as indices) is the next
//! step, and `cdhit_backend.rs` records its requirements.
//!
//! # Why argv, rather than building `Options` by hand
//!
//! `Options::set_options` is cd-hit's own argument parser, translated. Handing it the same
//! argv the subprocess path would have handed the shell means the two backends cannot drift
//! as flags change, and it keeps the command strings in [`crate::cdhit`] as the single
//! specification. Setting `Options` fields directly would create a second, silently
//! divergent definition of what each Panaroo call means.
//!
//! # Splitting the command string on whitespace
//!
//! [`argv_from_command`] splits on whitespace, which would be wrong if any path could
//! contain a space. It cannot matter here: the call sites interpolate paths unquoted
//! (`-i {input_file}`), so a path with a space already breaks the subprocess path in exactly
//! the same way. The split therefore produces precisely the tokens `/bin/sh` would have
//! produced — no more, no less.

use cdhit_rs::cdhit_common::*;
use std::sync::Mutex;

/// cd-hit is not re-entrant: `set_options` mutates shared alphabet/score-matrix state and
/// the translation keeps the C++'s global-ish structure. Panaroo only ever calls it from
/// one thread, but that is a property of the caller, not of this function, so serialise.
static CDHIT_LOCK: Mutex<()> = Mutex::new(());

struct QuietOutputGuard {
    previous: bool,
}

impl QuietOutputGuard {
    fn new(quiet: bool) -> Self {
        let previous = quiet_output();
        set_quiet_output(quiet);
        Self { previous }
    }
}

impl Drop for QuietOutputGuard {
    fn drop(&mut self) {
        set_quiet_output(self.previous);
    }
}

/// Turn a shell command string into argv, dropping the output redirection.
///
/// `> /dev/null` is a shell construct; in-process there is no shell to interpret it, and
/// suppression is handled by the caller instead. See the module note on whitespace.
pub fn argv_from_command(cmd: &str) -> Vec<String> {
    cmd.split_whitespace()
        .take_while(|t| *t != ">")
        .map(|t| t.to_string())
        .collect()
}

/// `MAX_UAA` for the nucleotide alphabet.
///
/// Defined in `cdhit-rs/src/bin/cd-hit-est.rs:13` as a private `const` of that binary, so it
/// is not importable from the library. Duplicated here with its source noted; if `cdhit-rs`
/// ever exports it, use theirs.
const MAX_UAA_EST: i32 = 4;

fn num_procs() -> i32 {
    std::thread::available_parallelism()
        .map(|n| n.get() as i32)
        .unwrap_or(1)
}

/// `cd-hit` — the protein front-end, `src/bin/cd-hit.rs` translated to a callable.
///
/// Mirrors that `main` step for step. The only omissions are its `println!` banners and the
/// `print_usage` exit paths, which are process-level behaviour rather than clustering.
pub fn run_cd_hit_main(argv: &[String], quiet: bool) {
    let _guard = CDHIT_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let _quiet_guard = QuietOutputGuard::new(quiet);

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let mut mat = ScoreMatrix::new();
    let mut alphabet = Alphabet::new();

    if !options.set_options(argv, false, false, num_procs(), &mut mat, &mut alphabet) {
        panic!("cd-hit: could not parse options from {argv:?}");
    }
    options.validate();

    let db_in = options.input.clone();
    let db_out = options.output.clone();

    let naa_tab = Naa::new(MAX_UAA);
    options.naan = naa_tab.naan_array[options.naa as usize];
    seq_db.naan = naa_tab.naan_array[options.naa as usize];

    seq_db.read(&db_in, &options);
    seq_db.sort_divide(&mut options, true, &alphabet);
    seq_db.do_clustering(&options, &mat, &naa_tab, &[]);

    seq_db.write_clusters(&db_in, &db_out, &options);
    seq_db.write_extra_1d(&options);
}

/// `cd-hit-est` — the nucleotide front-end, `src/bin/cd-hit-est.rs` translated to a callable.
///
/// Differs from [`run_cd_hit_main`] exactly as the two front-ends differ: nucleotide
/// defaults set *before* parsing, `set_options(.., est = true)`, `MAX_UAA_EST`, and the
/// complementary word index when `-r` is non-zero. Panaroo never uses paired-end mode, but
/// the branch is kept so this stays a faithful mirror of the front-end.
pub fn run_cd_hit_est_main(argv: &[String], quiet: bool) {
    let _guard = CDHIT_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let _quiet_guard = QuietOutputGuard::new(quiet);

    let mut options = Options::default();
    let mut seq_db = SequenceDB::default();
    let mut mat = ScoreMatrix::new();
    let mut alphabet = Alphabet::new();

    // Set before parsing, so a command-line flag still overrides them.
    options.cluster_thd = 0.95;
    options.naa = 10;
    options.naa_top_limit = 12;
    alphabet.setaa_to_na();
    mat.set_to_na();

    if !options.set_options(argv, false, true, num_procs(), &mut mat, &mut alphabet) {
        panic!("cd-hit-est: could not parse options from {argv:?}");
    }
    options.validate();

    let db_in = options.input.clone();
    let db_in_pe = options.input_pe.clone();
    let db_out = options.output.clone();
    let db_out_pe = options.output_pe.clone();

    let naa_tab = Naa::new(MAX_UAA_EST);
    options.naan = naa_tab.naan_array[options.naa as usize];
    seq_db.naan = naa_tab.naan_array[options.naa as usize];

    let mut comp_aan_idx: Vec<i32> = Vec::new();
    if options.option_r != 0 {
        comp_aan_idx.resize(seq_db.naan as usize, 0);
        make_comp_short_word_index(options.naa, &naa_tab.naan_array, &mut comp_aan_idx);
    }

    if options.pe_mode != 0 {
        seq_db.read_pe(&db_in, &db_in_pe, &options);
    } else {
        seq_db.read(&db_in, &options);
    }
    seq_db.sort_divide(&mut options, true, &alphabet);
    seq_db.do_clustering(&options, &mat, &naa_tab, &comp_aan_idx);

    if options.pe_mode != 0 {
        seq_db.write_clusters_pe(&db_in, &db_in_pe, &db_out, &db_out_pe, &options);
    } else {
        seq_db.write_clusters(&db_in, &db_out, &options);
    }
    seq_db.write_extra_1d(&options);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn argv_drops_the_shell_redirection() {
        // exactly the shape `run_cdhit` builds
        let cmd = "cd-hit -T 4 -i IN -o OUT -c 0.98 -s None -M 0 -d 999 -g 1 -n 2 > /dev/null";
        let argv = argv_from_command(cmd);
        assert_eq!(argv[0], "cd-hit");
        assert_eq!(argv.last().unwrap(), "2");
        assert!(!argv.iter().any(|t| t.contains("dev/null")));
        // `-s None` survives as a token: the flag really does receive the string "None".
        let i = argv.iter().position(|t| t == "-s").unwrap();
        assert_eq!(argv[i + 1], "None");
    }

    #[test]
    fn argv_is_unchanged_when_there_is_no_redirection() {
        let argv = argv_from_command("cd-hit-est -T 1 -i A -o B");
        assert_eq!(argv, ["cd-hit-est", "-T", "1", "-i", "A", "-o", "B"]);
    }
}
