//! `Bio.Seq` — the handful of sequence operations Panaroo uses.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [Biopython](https://biopython.org/) 1.84 (`Bio.Seq`), Biopython License Agreement /
//! BSD 3-Clause. No Biopython source was copied; the rules below were established by
//! running the library. The genetic code data lives in [`super::codon_table`], generated
//! from Biopython. See `NOTICE.md`.
//!
//! # Two translators, not interchangeable
//!
//! Panaroo has its own numpy lookup-table translator, duplicated as `prokka::translate` and
//! `generate_alignments::translate`. Those are Panaroo functions and live in those modules.
//! [`translate`] here is **Biopython's**, used by `find_missing::translate_to_match` and the
//! codon-alignment path. They disagree: Panaroo's honours a start codon by substituting `M`
//! for the first residue and uses the table selected by `--codon-table`; Biopython's does
//! neither and always uses NCBI table 1. Keep the call sites straight.

/// `Bio.Seq.reverse_complement(s)`
///
/// IUPAC-aware and case-preserving; any character with no complement (`-`, `*`, digits) is
/// passed through unchanged. OBSERVED: `"ACGTRYKM"` -> `"KMRYACGT"`, `"acgt"` -> `"acgt"`,
/// `"ACGT-N"` -> `"N-ACGT"`.
pub fn reverse_complement(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| complement_base(b) as char)
        .collect()
}

/// The complement of one base, preserving case. IUPAC ambiguity codes included.
fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'U' => b'A',
        b'u' => b'a',
        b'M' => b'K',
        b'K' => b'M',
        b'm' => b'k',
        b'k' => b'm',
        b'R' => b'Y',
        b'Y' => b'R',
        b'r' => b'y',
        b'y' => b'r',
        b'W' => b'W',
        b'S' => b'S',
        b'w' => b'w',
        b's' => b's',
        b'B' => b'V',
        b'V' => b'B',
        b'b' => b'v',
        b'v' => b'b',
        b'D' => b'H',
        b'H' => b'D',
        b'd' => b'h',
        b'h' => b'd',
        b'N' => b'N',
        b'n' => b'n',
        other => other,
    }
}

/// `Bio.Seq.translate(s)` — NCBI table 1, `*` for stop, `X` for an ambiguous codon.
///
/// A trailing partial codon is dropped (Biopython emits a `BiopythonWarning` and continues).
/// OBSERVED: `"ATGAAAT"` -> `"MK"`, `"ATGNNNAAA"` -> `"MXK"`, `"atgaaa"` -> `"MK"`.
///
/// Unlike `prokka::translate`, this does **not** substitute `M` for a start codon.
pub fn translate(s: &str) -> String {
    let table = super::codon_table::generic_by_id(1).expect("NCBI table 1");
    let b = s.as_bytes();
    let mut out = String::with_capacity(b.len() / 3);
    for codon in b.chunks_exact(3) {
        let up: String = codon
            .iter()
            .map(|&c| c.to_ascii_uppercase() as char)
            .collect();
        out.push(translate_codon(table, &up));
    }
    out
}

/// One codon, with Biopython's IUPAC ambiguity handling.
///
/// **This is the subtle part.** Biopython does not simply give `X` for any codon containing
/// an ambiguity code: it expands the ambiguity and, if every possibility yields the *same*
/// residue, returns that residue. OBSERVED against Biopython 1.84:
///
/// ```text
/// CGN -> R     (CGA/CGC/CGG/CGT are all Arg)
/// GCN -> A     (all Ala)
/// CGR -> R
/// AAN -> X     (Lys or Asn)
/// TTN -> X
/// NNN -> X
/// ```
///
/// Getting this wrong is not cosmetic: `translate_to_match` writes the result into
/// `combined_protein_CDS.fasta` and `gene_data.csv` for every refound gene.
fn translate_codon(table: &super::codon_table::CodonTable, codon: &str) -> char {
    let lookup = |c: &str| -> Option<char> {
        if table.stop_codons.contains(&c) {
            Some('*')
        } else {
            table
                .forward_table
                .iter()
                .find(|(k, _)| *k == c)
                .map(|&(_, aa)| aa)
        }
    };

    if let Some(aa) = lookup(codon) {
        return aa;
    }

    // expand IUPAC ambiguity codes
    let sets: Vec<&str> = codon.chars().map(iupac_bases).collect();
    if sets.iter().any(|s| s.is_empty()) {
        return 'X';
    }
    let mut resolved: Option<char> = None;
    for a in sets[0].chars() {
        for b in sets[1].chars() {
            for c in sets[2].chars() {
                let expanded: String = [a, b, c].iter().collect();
                match lookup(&expanded) {
                    Some(aa) => match resolved {
                        None => resolved = Some(aa),
                        Some(prev) if prev == aa => {}
                        Some(_) => return 'X',
                    },
                    None => return 'X',
                }
            }
        }
    }
    resolved.unwrap_or('X')
}

/// The bases an IUPAC nucleotide code stands for. Empty for anything unrecognised.
fn iupac_bases(c: char) -> &'static str {
    match c.to_ascii_uppercase() {
        'A' => "A",
        'C' => "C",
        'G' => "G",
        'T' | 'U' => "T",
        'R' => "AG",
        'Y' => "CT",
        'S' => "CG",
        'W' => "AT",
        'K' => "GT",
        'M' => "AC",
        'B' => "CGT",
        'D' => "AGT",
        'H' => "ACT",
        'V' => "ACG",
        'N' => "ACGT",
        _ => "",
    }
}

/// The genetic code tables live in [`super::codon_table`], not here: they are Biopython
/// *data* rather than a function, and third-party data belongs in its own attributed file.
pub use super::codon_table::{generic_by_id, CodonTable};

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values from running Biopython 1.84.

    #[test]
    fn reverse_complement_matches_biopython() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("acgt"), "acgt");
        assert_eq!(reverse_complement("ACGTN"), "NACGT");
        assert_eq!(reverse_complement("ACGTRYKM"), "KMRYACGT");
        assert_eq!(reverse_complement("ACGT-N"), "N-ACGT");
    }

    #[test]
    fn translate_matches_biopython() {
        assert_eq!(translate("ATGAAATAA"), "MK*");
        assert_eq!(translate("ATGAAAT"), "MK"); // trailing partial codon dropped
        assert_eq!(translate("ATGNNNAAA"), "MXK");
        assert_eq!(translate("ATGTGA"), "M*");
        assert_eq!(translate("atgaaa"), "MK");
    }

    #[test]
    fn translate_resolves_iupac_ambiguity_like_biopython() {
        // Biopython expands ambiguity and returns the residue when all options agree.
        assert_eq!(translate("CGN"), "R");
        assert_eq!(translate("GCN"), "A");
        assert_eq!(translate("CGR"), "R");
        // ... and X when they do not
        assert_eq!(translate("AAN"), "X");
        assert_eq!(translate("TTN"), "X");
        assert_eq!(translate("NNN"), "X");
        assert_eq!(translate("ATN"), "X");
    }

    #[test]
    fn translate_does_not_do_panaroos_start_codon_substitution() {
        // prokka::translate turns a leading TTG into M; Bio.Seq.translate does not.
        assert_eq!(translate("TTGAAA"), "LK");
    }
}
