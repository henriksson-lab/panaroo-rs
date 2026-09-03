//! Run MAFFT in-process, via the [`rust-MAFFT`] translation, instead of spawning the
//! `mafft` shell wrapper.
//!
//! Enabled by the `mafft-embedded` Cargo feature, which is off by default.
//!
//! [`rust-MAFFT`]: https://github.com/luksgrin/rust-MAFFT
//!
//! # Why this is worth doing
//!
//! `--alignment core` on the `ci` dataset issues **~5,100 mafft invocations**, one per gene
//! cluster. Upstream `mafft` is a ~3,000-line shell script that then execs its own binaries,
//! so each of those is a shell startup plus a process spawn for an alignment of a handful of
//! short sequences. That fixed overhead, repeated 5,100 times, is most of what the alignment
//! stage costs once `output_sequence` stopped re-parsing the combined FASTA.
//!
//! # Why argv, and not the `MafftEngine` API
//!
//! `mafft` crate exposes `MafftEngine::new(AlignmentMode::FftNs2).align(..)`, which takes the
//! alignment mode as an argument. But Panaroo passes `--auto`, and `--auto` *selects* that
//! mode from the sequence count and length; `--adjustdirection` likewise runs k-mer strand
//! detection before alignment. Calling the engine directly would mean re-deriving both here,
//! which is exactly how a caller silently diverges from C MAFFT. `run_from` takes the same
//! argv the shell would have received, so the flag semantics stay single-sourced in the
//! translated CLI layer.
//!
//! # Output
//!
//! `mafft` writes the alignment to stdout and Panaroo captures those bytes verbatim into
//! `<gene>.aln.fas`. `run_from` takes a `&mut dyn Write`, so the bytes are collected into a
//! buffer rather than going through a pipe.
//!
//! # Case
//!
//! C MAFFT folds case at *read* time — nucleotide input to lowercase, protein to uppercase
//! (`core/io.c:1462-1467`), driven by the sequence type, which `--nuc`/`--amino` force. The
//! fork reproduces that; before it did, DNA output came back uppercased and every one of the
//! ~5,100 comparisons would have failed on case alone.

/// Turn the assembled shell command into argv.
///
/// Panaroo builds `mafft --auto --adjustdirection --thread 1 --nuc {file}` (DNA) or
/// `mafft --auto --amino {file}` (protein) by interpolation, with no quoting — so a path
/// containing a space already breaks the subprocess path in the same way, and splitting on
/// whitespace yields exactly the tokens `/bin/sh` would have produced.
pub fn argv_from_command(cmd: &str) -> Vec<String> {
    cmd.split_whitespace().map(|t| t.to_string()).collect()
}

/// Run one MAFFT command in-process, returning what it would have written to stdout.
///
/// Panics with the message MAFFT would have printed, mirroring the subprocess path's
/// failure mode rather than silently producing an empty alignment.
pub fn run_mafft(cmd: &str) -> Vec<u8> {
    let argv = argv_from_command(cmd);
    let mut out: Vec<u8> = Vec::new();
    // Silence MAFFT's progress lines (`Alignment: N columns`, strategy banners). The
    // subprocess path captured and discarded them via `popen_communicate`; in-process they
    // would otherwise reach the user's terminal ~5,100 times from 20 threads. Errors are
    // unaffected -- they come back through `MafftError`, never through the sink.
    let quiet = mafft_rs::progress::SilentProgress;
    match mafft_rs::run_from_with_progress(&argv, &mut out, &quiet) {
        Ok(()) => out,
        Err(e) => panic!("RuntimeError: mafft failed ({}): {}", e.code(), e.message()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn argv_matches_the_tokens_the_shell_would_produce() {
        let argv = argv_from_command("mafft --auto --adjustdirection --thread 1 --nuc /tmp/g.fa");
        assert_eq!(
            argv,
            [
                "mafft",
                "--auto",
                "--adjustdirection",
                "--thread",
                "1",
                "--nuc",
                "/tmp/g.fa"
            ]
        );
        // the protein form Panaroo uses for codon realignment
        assert_eq!(
            argv_from_command("mafft --auto --amino /tmp/p.fa"),
            ["mafft", "--auto", "--amino", "/tmp/p.fa"]
        );
    }
}
