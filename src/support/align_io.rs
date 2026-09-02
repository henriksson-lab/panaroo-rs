//! `Bio.AlignIO` and `Bio.Align.MultipleSeqAlignment`.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [Biopython](https://biopython.org/) 1.84 (`Bio.AlignIO`, `Bio.Align`), Biopython
//! License Agreement / BSD 3-Clause. No Biopython source was copied. See `NOTICE.md`.
//!
//! Used by the alignment stage: reading aligner output back in, reordering a protein
//! alignment to match DNA record order, and writing the concatenated core alignment.
//! `AlignIO` is `SeqIO` plus one rule: **every record must be the same length**, and
//! `AlignIO.read` raises otherwise. That check is what distinguishes
//! `generate_alignments::is_valid_alignment` from `is_valid_fasta`.

use super::seqio::SeqRecord;

/// `Bio.Align.MultipleSeqAlignment`
#[derive(Debug, Clone, Default)]
pub struct MultipleSeqAlignment {
    pub records: Vec<SeqRecord>,
}

impl MultipleSeqAlignment {
    pub fn new(records: Vec<SeqRecord>) -> Self {
        MultipleSeqAlignment { records }
    }

    /// `len(aln)` — number of records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// `aln.get_alignment_length()`
    pub fn get_alignment_length(&self) -> usize {
        self.records.first().map(|r| r.seq.len()).unwrap_or(0)
    }
}

/// `AlignIO.read(path, "fasta")` — panics on ragged records, as Biopython raises.
pub fn read_fasta_alignment(path: &str) -> MultipleSeqAlignment {
    let records = super::seqio::parse_fasta_file(path);
    if let Some(first) = records.first() {
        if records.iter().any(|r| r.seq.len() != first.seq.len()) {
            panic!("ValueError: Sequences must all be the same length ({path})");
        }
    }
    MultipleSeqAlignment { records }
}

/// `AlignIO.write(aln, path, "fasta")`
pub fn write_fasta_alignment(aln: &MultipleSeqAlignment, path: &str) {
    super::seqio::write_fasta_file(&aln.records, path);
}
