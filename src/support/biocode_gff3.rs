//! biocode's GFF3 writer.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable output of
//! [biocode](https://github.com/jorvis/biocode) 0.12.1 (MIT), specifically
//! `biocode.things.Gene.print_as(source='GenBank', format='gff3')` and the
//! `biocode.utils.wrapped_fasta` helper. No biocode source was copied; the format below was
//! established by running the real converter over a RefSeq GenBank file and reading the
//! output. See `NOTICE.md`.
//!
//! biocode is a **runtime dependency of Panaroo**, not just a vendored file:
//! `panaroo/biocode_convert.py` imports `biocode.annotation`, `biocode.things` and
//! `biocode.utils` at module scope, and `prokka.py` imports that module, so
//! `python -m panaroo` fails at import time without it.
//!
//! # Format reproduced (OBSERVED on GCF_000759575.2)
//!
//! One block per gene, tab-separated, with `source` = `GenBank`, coordinates 1-based
//! inclusive (`fmin + 1 .. fmax`), and `.` for score and phase except on `CDS`:
//!
//! ```text
//! seqid  GenBank  gene         s  e  .  ±  .  ID={lt};locus_tag={lt}
//! seqid  GenBank  mRNA         s  e  .  ±  .  ID={lt}.mRNA.{n};Parent={lt}
//! seqid  GenBank  CDS          s  e  .  ±  p  ID={rna}.CDS.{k};Parent={rna}
//! seqid  GenBank  exon         s  e  .  ±  .  ID={rna}.exon.{k};Parent={rna}
//! seqid  GenBank  polypeptide  s  e  .  ±  .  ID={lt}.polypeptide.{n};Parent={rna}[;gene_symbol=..][;product_name=..]
//! ```
//!
//! and for non-coding RNA genes:
//!
//! ```text
//! seqid  GenBank  rRNA  s  e  .  ±  .  ID={lt}.rRNA.{n};Parent={lt};product_name=..
//! seqid  GenBank  tRNA  s  e  .  ±  .  ID={lt}.tRNA.{n};Parent={lt};anticodon=..
//! ```
//!
//! An absent `gene_symbol` or `product_name` is omitted entirely rather than written empty.

use std::fmt::Write as _;

/// One GFF3 row.
pub struct Row {
    pub seqid: String,
    pub feature_type: String,
    /// 0-based inclusive start; written as `fmin + 1`.
    pub fmin: i64,
    /// Exclusive end; written as-is.
    pub fmax: i64,
    pub strand: char,
    /// Column 8. `None` renders as `.`.
    pub phase: Option<u8>,
    /// Column 9, already in the order biocode emits.
    pub attributes: Vec<(String, String)>,
}

impl Row {
    pub fn render(&self) -> String {
        let mut attrs = String::new();
        for (i, (k, v)) in self.attributes.iter().enumerate() {
            if i > 0 {
                attrs.push(';');
            }
            let _ = write!(attrs, "{k}={v}");
        }
        format!(
            "{}\tGenBank\t{}\t{}\t{}\t.\t{}\t{}\t{}",
            self.seqid,
            self.feature_type,
            self.fmin + 1,
            self.fmax,
            self.strand,
            match self.phase {
                Some(p) => p.to_string(),
                None => ".".to_string(),
            },
            attrs
        )
    }
}

/// `biocode.utils.wrapped_fasta(residues)` — 60 characters per line, no trailing newline
/// (the caller adds one).
pub fn wrapped_fasta(residues: &str) -> String {
    let b = residues.as_bytes();
    let mut out = String::with_capacity(b.len() + b.len() / 60 + 1);
    for (i, chunk) in b.chunks(60).enumerate() {
        if i > 0 {
            out.push('\n');
        }
        out.push_str(std::str::from_utf8(chunk).expect("ASCII residues"));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_renders_like_biocode() {
        // OBSERVED first line of the reference conversion of GCF_000759575.2:
        //   NZ_CP073712\tGenBank\tgene\t1\t1323\t.\t+\t.\tID=LT41_RS00005;locus_tag=LT41_RS00005
        let r = Row {
            seqid: "NZ_CP073712".into(),
            feature_type: "gene".into(),
            fmin: 0,
            fmax: 1323,
            strand: '+',
            phase: None,
            attributes: vec![
                ("ID".into(), "LT41_RS00005".into()),
                ("locus_tag".into(), "LT41_RS00005".into()),
            ],
        };
        assert_eq!(
            r.render(),
            "NZ_CP073712\tGenBank\tgene\t1\t1323\t.\t+\t.\tID=LT41_RS00005;locus_tag=LT41_RS00005"
        );
    }

    #[test]
    fn cds_carries_a_phase() {
        let r = Row {
            seqid: "c".into(),
            feature_type: "CDS".into(),
            fmin: 0,
            fmax: 9,
            strand: '-',
            phase: Some(0),
            attributes: vec![("ID".into(), "x".into())],
        };
        assert_eq!(r.render(), "c\tGenBank\tCDS\t1\t9\t.\t-\t0\tID=x");
    }

    #[test]
    fn fasta_wraps_at_sixty() {
        assert_eq!(
            wrapped_fasta(&"A".repeat(130))
                .lines()
                .map(|l| l.len())
                .collect::<Vec<_>>(),
            [60, 60, 10]
        );
        assert_eq!(wrapped_fasta("ACGT"), "ACGT");
    }
}
