//! An in-memory alternative path into cd-hit.
//!
//! # Status: specification, not yet wired in
//!
//! Nothing in the pipeline calls this yet. It exists so that [`cdhit-rs`] has a concrete,
//! compile-checkable target to build against while it is still being translated, and so
//! that the shape of the seam is settled before either side commits to it. The live path
//! remains [`crate::cdhit::run_cdhit`] / [`run_cdhit_est`], which shell out.
//!
//! [`cdhit-rs`]: https://github.com/henriksson-lab/cdhit-rs
//!
//! # Why an in-memory path is worth having
//!
//! Today every cd-hit round round-trips through the filesystem. In `iterative_cdhit` that
//! is entirely gratuitous: it writes a temporary FASTA, cd-hit reads it and writes a
//! representative FASTA plus a `.clstr`, the next round reads that back, and only the final
//! `.clstr` is parsed. **None of those files is a Panaroo deliverable** — they are pure
//! interprocess plumbing, repeated once per threshold (11 rounds on a default run).
//!
//! The goal is to write nothing to disk on this path. Spilling to a file is acceptable
//! where memory genuinely forces it, but it should be a deliberate fallback, not the
//! default mechanism.
//!
//! # Who actually needs the files
//!
//! Established by reading the call sites rather than assumed:
//!
//! - **`iterative_cdhit` (11 of the 12 invocations per run): nothing needs any file.** The
//!   temp input, the per-round representatives FASTA and the intermediate `.clstr` are pure
//!   interprocess plumbing. Round N+1 reads round N's representatives off disk only because
//!   the two rounds are separate processes; in-process it can take them from the previous
//!   result.
//! - **The initial clustering (1 invocation) writes two deliverables.**
//!   `combined_protein_cdhit_out.txt` and its `.clstr` are among the 13 files the parity
//!   suite compares byte-for-byte, so they must still appear on disk with identical bytes.
//!
//! Note what is *not* true: `__main__.py:378` passes only `cd_hit_out + ".clstr"` to
//! `generate_network`. **Panaroo never reads the representatives FASTA it just produced.**
//! It is an unread by-product that happens to sit in the output directory. Its bytes matter
//! for the parity claim and for nothing else.
//!
//! # Writing the deliverables is Panaroo's job, not the backend's
//!
//! Hence the API below takes no output paths and writes nothing. `cdhit-core`'s
//! `write_clusters` produces the representatives FASTA by seeking back into the *input* file
//! and copying each representative's record verbatim (`des_begin` + `tot_length` bytes) —
//! which is why it needs the source to still exist and be seekable. An embedded backend
//! should not do that at all: it returns which sequences are representatives, and the caller
//! reconstructs the file if it wants one.
//!
//! That is cheap here because the caller already has everything needed. The 60-column
//! wrapping in `combined_protein_cdhit_out.txt` is not cd-hit's: it was copied verbatim out
//! of `combined_protein_CDS.fasta`, which **Panaroo itself wrote**. Reproducing it means
//! calling the same writer again, not reverse-engineering a foreign format.
//!
//! # Integration-level speedups, ranked — and an honest ceiling
//!
//! Measured inputs on the `ci` dataset: 5,362 centroids reach `iterative_cdhit`; the DNA
//! rounds carry ~5.1 MB of sequence (4 rounds), the protein rounds ~1.7 MB (7 rounds), and
//! a `.clstr` for ~5,300 clusters is ~675 KB. Eleven rounds per run.
//!
//! Per round the same sequence data is currently: serialised to FASTA text (by cd-hit's
//! verbatim record copy), read back, re-parsed, and re-encoded to residue indices — and the
//! cluster assignment, which cd-hit already holds in `Sequence::cluster_id`, is formatted to
//! `.clstr` text and immediately re-parsed by `parse_clstr_strings`.
//!
//! 1. **Keep one `SequenceDB` alive across the rounds.** Round N+1's input is exactly round
//!    N's representatives, so writing them out and reading them back is pure round-tripping.
//!    Filtering the existing db to `rep_seqs` removes, for 10 of 11 rounds: the verbatim copy
//!    (which re-reads the input a second time), the file read, the FASTA parse, and the
//!    residue encode. This is the largest of these items and the one that needs the most care
//!    — `write_clusters` emits representatives in ascending original index, so the filtered db
//!    must preserve that order before `sort_divide` re-sorts by length.
//! 2. **Return cluster assignments as data.** `Sequence::cluster_id` is already the answer;
//!    formatting ~675 KB of `.clstr` and parsing it back (allocating a `String` per member,
//!    ~5,300 x 11 = ~58k allocations) is avoidable in full.
//! 3. **Skip the round-1 input file.** `iterative_cdhit` builds the whole FASTA in a `String`
//!    and writes it; cd-hit then parses it straight back.
//! 4. **`clusters = clust_dict.values().cloned().collect()`** in `iterative_cdhit` deep-clones
//!    every cluster every round — another ~58k `String` clones. Purely Panaroo-side, and
//!    independent of any cd-hit change.
//!
//! **Ceiling: all four together are on the order of a few hundred milliseconds of a ~50 s
//! run — well under 1%.** Everything above is O(total residues): parsing, encoding and text
//! formatting are linear in ~25 MB, while `do_clustering` is O(pairs x alignment) and is what
//! actually dominates. Integration cannot change that ratio; only cd-hit's internal work can.
//!
//! So these are worth doing for the reasons the rest of the port's optimisations are — free,
//! provable, and linear in genome count — and **not** because they will move the benchmark.
//! Anyone hoping for a step change should look at cd-hit's internals, not at this seam.
//!
//! One thing that is NOT an integration win, despite looking like one: lowering cd-hit's `-T`
//! independently of Panaroo's thread count. cd-hit's CPU cost is superlinear in threads
//! (9.67 CPU-s at `-T 1` vs 19.96 at `-T 20`), so capping it would cut total CPU — but it
//! changes the flags, and therefore diverges from the Python. It would have to be an opt-in
//! flag, not a default.
//!
//! # Fidelity requirement
//!
//! The command strings built in [`crate::cdhit`] are the specification, not this struct.
//! Whatever an embedded backend does must be what cd-hit would have done given exactly those
//! flags — including the two upstream oddities the port preserves deliberately: `-s None`
//! (`--family_len_dif_percent` has no argparse default) and `aS=AS` passing an **int** where
//! a float is expected, which changes the rendered flag. Hence [`PyNum`] here rather than
//! `f64`: the distinction is observable in the command line and must survive into any
//! backend that claims to reproduce it.
//!
//! A backend is therefore best implemented by handing cd-hit's own argument parser the same
//! argv the subprocess path would build, rather than by setting fields by hand. That way the
//! two backends cannot drift apart as flags change.
//!
//! # Threading
//!
//! Upstream cd-hit is **nondeterministic when run multithreaded**, and `cdhit-rs`, being a
//! faithful translation, inherits that. So `threads > 1` output is not reproducible and
//! cannot be validated against anything: under `-T > 1` the C++ has no single correct
//! answer, and a translation that produced one would be diverging, not improving. `-T 1` is
//! the only configuration in which an embedded backend can be shown correct. See the README
//! section "A second source of nondeterminism".

use crate::support::pyfmt::PyNum;

/// Which cd-hit executable the request corresponds to.
///
/// Panaroo only ever uses these two: `cd-hit` for protein and `cd-hit-est` for nucleotide.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqKind {
    /// `cd-hit`
    Protein,
    /// `cd-hit-est`
    Nucleotide,
}

/// One cd-hit invocation, expressed as data rather than as a command line.
///
/// Field names follow the Python parameter names in `panaroo/cdhit.py`, not cd-hit's flag
/// letters, so a reader can line this up against the call sites. The flag each maps to is
/// given in the doc comment.
#[derive(Debug, Clone)]
pub struct CdhitRequest<'a> {
    /// Input sequences as `(name, sequence)`, **in input order**.
    ///
    /// Order is load-bearing: cd-hit sorts by decreasing length but its tie-breaking falls
    /// back on input order, so shuffling this changes which sequence becomes a cluster
    /// representative. Reproducing the file the subprocess path writes means reproducing
    /// this order exactly.
    ///
    /// Borrowed rather than owned **only because the subprocess backend must re-serialise
    /// these to a file anyway**. An embedded backend wants ownership: cd-hit's residue
    /// encoding is a length-preserving byte-to-byte map, so given an owned `String` it can
    /// `into_bytes()` and encode **in place**, with no allocation and no copy. If the
    /// embedded path becomes the default, change this to `Vec<(String, String)>` — the
    /// borrow is the one thing here that forces a copy that need not exist.
    pub sequences: &'a [(String, String)],
    /// Selects `cd-hit` vs `cd-hit-est`.
    pub kind: SeqKind,
    /// `-c`
    pub id: f64,
    /// `-T`. See the threading note above: only `1` is verifiable.
    pub n_cpu: i64,
    /// `-s`. `PyNum::None` renders as the literal `None`, which is what a default run sends.
    pub s: PyNum,
    /// `-aL`
    pub a_l: PyNum,
    /// `-AL`
    pub big_a_l: i64,
    /// `-aS`. May legitimately hold an `Int` — see the fidelity note above.
    pub a_s: PyNum,
    /// `-AS`
    pub big_a_s: i64,
    /// `accurate=True` sets `-g 1 -n 2`, which is why the low-identity rounds are expensive.
    pub accurate: bool,
    /// `use_local=True` sets `-G 0`.
    pub use_local: bool,
    /// `-n`, when the caller pins it rather than letting cd-hit choose from `id`.
    pub word_length: Option<i64>,
    /// `-l`. Protein only.
    pub min_length: Option<i64>,
    /// `-r`. Nucleotide only.
    pub strand: Option<i64>,
    /// `-p`. Nucleotide only.
    pub print_aln: bool,
    /// `-mask NX`. Nucleotide only.
    pub mask: bool,
}

/// The result of one invocation.
#[derive(Debug, Clone, Default)]
pub struct CdhitResult {
    /// Representatives as **indices into [`CdhitRequest::sequences`]**, in the order cd-hit
    /// would have written them — which is ascending original input index, since
    /// `write_clusters` sorts by `Sequence::index` before copying records out.
    ///
    /// Indices rather than copied strings: the caller already owns the sequences, so this
    /// costs `size_of::<usize>()` per representative instead of two `String` clones, and it
    /// is everything the caller needs both to feed the next round and to reconstruct the
    /// representatives FASTA byte-identically.
    pub representatives: Vec<usize>,
    /// Clusters in `.clstr` order, each listing its members in `.clstr` order.
    ///
    /// This is what `parse_cdhit_clusters` currently reconstructs by parsing text.
    pub clusters: Vec<Vec<String>>,
    /// The verbatim `.clstr` text, when the caller must persist a byte-identical artefact.
    ///
    /// `None` when the caller does not need it. The initial clustering in `main.rs` does —
    /// `combined_protein_cdhit_out.txt.clstr` is a compared output file — and asking the
    /// backend for the exact bytes is safer than re-serialising [`Self::clusters`] and
    /// hoping the formatting matches.
    pub clstr_text: Option<String>,
}

/// How a [`CdhitRequest`] gets executed.
///
/// Two implementations are intended: the current subprocess path, and an embedded one
/// calling `cdhit-core` directly. Keeping both behind one trait is what makes them
/// differentially testable — the same request through both must give the same result, which
/// is a far sharper test than comparing whole-pipeline output.
pub trait CdhitBackend {
    /// Run one clustering. Writes nothing to disk.
    ///
    /// Spilling to a temporary file is acceptable where memory genuinely forces it, but it
    /// must be a deliberate fallback. Note that cd-hit's own spill mechanism cannot be
    /// relied on for this: `Sequence::swap` / `-B store_disk` are **dead code in 4.8.1**
    /// (nothing ever assigns `swap`), and Panaroo passes `-M 0`, so upstream holds the
    /// entire database in RAM unconditionally.
    fn run(&self, req: &CdhitRequest<'_>) -> CdhitResult;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The default protein round from `iterative_cdhit`, as data.
    ///
    /// Pins the two upstream oddities so a backend implementer sees them in a test rather
    /// than having to find them in a doc comment: `s` is `None`, and `a_s` is an **Int**.
    #[test]
    fn a_default_protein_round_carries_the_upstream_oddities() {
        let seqs = vec![("c1".to_string(), "MKV".to_string())];
        let req = CdhitRequest {
            sequences: &seqs,
            kind: SeqKind::Protein,
            id: 0.99,
            n_cpu: 1,
            s: PyNum::None,
            a_l: PyNum::Float(0.0),
            big_a_l: 99999999,
            a_s: PyNum::Int(99999999),
            big_a_s: 99999999,
            accurate: true,
            use_local: false,
            word_length: None,
            min_length: None,
            strand: None,
            print_aln: false,
            mask: false,
        };
        // `-s None`, not `-s 0.0`: the flag really does receive the string "None".
        assert_eq!(req.s.to_py_string(), "None");
        // `-aS 99999999`, not `-aS 99999999.0`: the int/float distinction is observable.
        assert_eq!(req.a_s.to_py_string(), "99999999");
        assert_eq!(req.kind, SeqKind::Protein);
    }
}
