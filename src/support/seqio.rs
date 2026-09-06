//! `Bio.SeqIO` and `Bio.SeqRecord` — FASTA reading and writing.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [Biopython](https://biopython.org/) 1.84 (`Bio.SeqIO`, `Bio.SeqRecord`), Biopython
//! License Agreement / BSD 3-Clause. No Biopython source was copied; every rule below was
//! established by running the library and recording what it produced. See `NOTICE.md`.
//!
//! # Behaviour reproduced (all OBSERVED against Biopython 1.84)
//!
//! **Parsing** (`SeqIO.parse(handle, "fasta")`):
//!   - `record.id` is the header up to the first whitespace
//!   - `record.description` is the **whole** header line after `>`, id included — so for
//!     `>abc def ghi` the id is `abc` and the description is `abc def ghi`
//!   - `record.name` equals the id
//!   - a header with no description gives `description == id`, not `""`
//!   - sequence lines are concatenated with no separator
//!
//! **Writing** (`SeqIO.write(records, handle, "fasta")`):
//!   - sequence lines wrap at **60** characters
//!   - the header is chosen by this rule, which is easy to get wrong:
//!     - `description` empty  -> `>{id}`
//!     - `description` starts with a token equal to `id` -> `>{description}`
//!     - otherwise -> `>{id} {description}`
//!
//!     So `id="i3", description="i3 extra words"` writes `>i3 extra words` (not
//!     `>i3 i3 extra words`), and `id="i4", description="other text"` writes
//!     `>i4 other text`.
//!
//! Every FASTA Panaroo writes sets `description == id`, so in practice the header is always
//! `>{id}` — but the general rule is implemented because `output_sequence` builds ids of the
//! form `{isolate};{seq_id}` and a future change to descriptions would otherwise silently
//! alter every output file.

use std::io::Write;

/// `Bio.SeqRecord.SeqRecord`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SeqRecord {
    pub id: String,
    pub name: String,
    pub description: String,
    pub seq: String,
    /// `record.annotations` — Panaroo only ever stores `"scaffold"`.
    pub annotations: Vec<(String, String)>,
}

impl SeqRecord {
    /// `SeqRecord(seq, id=..., description=...)` — `name` defaults to the id, as Biopython
    /// does when `name` is not given.
    pub fn new(seq: String, id: String, description: String) -> Self {
        SeqRecord {
            name: id.clone(),
            id,
            description,
            seq,
            annotations: Vec::new(),
        }
    }

    /// `record.annotations["scaffold"]`
    pub fn scaffold(&self) -> &str {
        self.annotations
            .iter()
            .find(|(k, _)| k == "scaffold")
            .map(|(_, v)| v.as_str())
            .expect("KeyError: 'scaffold'")
    }

    pub fn set_annotation(&mut self, key: &str, value: &str) {
        self.annotations.push((key.to_string(), value.to_string()));
    }
}

/// `SeqIO.parse(handle, "fasta")`
pub fn parse_fasta(text: &str) -> Vec<SeqRecord> {
    let mut out: Vec<SeqRecord> = Vec::new();
    let mut cur: Option<SeqRecord> = None;
    for line in text.lines() {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if let Some(header) = line.strip_prefix('>') {
            if let Some(r) = cur.take() {
                out.push(r);
            }
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            cur = Some(SeqRecord {
                name: id.clone(),
                id,
                description: header.to_string(),
                seq: String::new(),
                annotations: Vec::new(),
            });
        } else if let Some(r) = cur.as_mut() {
            r.seq.push_str(line.trim());
        }
        // Text before the first '>' is ignored, as Biopython does.
    }
    if let Some(r) = cur {
        out.push(r);
    }
    out
}

/// `SeqIO.parse(handle, "fasta-pearson")` — the tolerant variant `prokka.py` falls back to.
///
/// Differs from the strict parser only in accepting `*` and `-` inside sequences, which our
/// parser already does (it does no alphabet validation), so this delegates. Kept as a named
/// entry point so the call site reads like the Python.
pub fn parse_fasta_pearson(text: &str) -> Vec<SeqRecord> {
    parse_fasta(text)
}

/// `SeqIO.parse(path, "fasta")`
pub fn parse_fasta_file(path: &str) -> Vec<SeqRecord> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("could not read {path}: {e}"));
    parse_fasta(&text)
}

/// The header line Biopython writes for a record — see the module docs.
fn fasta_title(r: &SeqRecord) -> String {
    if r.description.is_empty() {
        r.id.clone()
    } else if r.description.split_whitespace().next() == Some(r.id.as_str()) {
        r.description.clone()
    } else {
        format!("{} {}", r.id, r.description)
    }
}

/// Biopython's FASTA line width.
pub const FASTA_WRAP: usize = 60;

/// `SeqIO.write(records, handle, "fasta")`
pub fn write_fasta<W: Write>(records: &[SeqRecord], out: &mut W) {
    for r in records {
        writeln!(out, ">{}", fasta_title(r)).expect("write");
        write_fasta_sequence(&r.seq, out);
    }
}

/// `SeqIO.write(records, path, "fasta")`
pub fn write_fasta_file(records: &[SeqRecord], path: &str) {
    let f = std::fs::File::create(path).unwrap_or_else(|e| panic!("could not create {path}: {e}"));
    let mut w = std::io::BufWriter::new(f);
    write_fasta(records, &mut w);
}

fn write_fasta_sequence<W: Write>(seq: &str, out: &mut W) {
    let b = seq.as_bytes();
    if b.is_empty() {
        writeln!(out).expect("write");
        return;
    }
    for chunk in b.chunks(FASTA_WRAP) {
        out.write_all(chunk).expect("write");
        out.write_all(b"\n").expect("write");
    }
}

/// Write one FASTA record whose description is empty, i.e. the header is exactly `>{id}`.
pub fn write_fasta_record<W: Write>(id: &str, seq: &str, out: &mut W) {
    writeln!(out, ">{id}").expect("write");
    write_fasta_sequence(seq, out);
}

/// `SeqIO.parse(open(path), "genbank")` — implemented in [`super::genbank`], which is its
/// own file because it is a separate third-party format reader.
pub use super::genbank::parse_genbank_file;

/// A GenBank record, reduced to what `biocode_convert` reads.
#[derive(Debug, Clone, Default)]
pub struct GenBankRecord {
    pub id: String,
    pub name: String,
    pub description: String,
    pub seq: String,
    pub features: Vec<GenBankFeature>,
    pub annotations: Vec<(String, String)>,
}

#[derive(Debug, Clone, Default)]
pub struct GenBankFeature {
    pub feature_type: String,
    /// `feature.location.start` — 0-based.
    pub start: i64,
    /// `feature.location.end` — exclusive.
    pub end: i64,
    /// `+1` or `-1`; `0` for an unstranded location, which `convert_gbk_gff3` rejects.
    pub strand: i32,
    /// `feature.location.parts` as `(start, end)` pairs, 0-based half-open.
    pub parts: Vec<(i64, i64)>,
    pub qualifiers: Vec<(String, Vec<String>)>,
}

impl GenBankFeature {
    /// `feature.qualifiers[key][0]`, or `None` when the key is absent.
    pub fn qual(&self, key: &str) -> Option<&str> {
        self.qualifiers
            .iter()
            .find(|(k, _)| k == key)
            .and_then(|(_, v)| v.first())
            .map(|s| s.as_str())
    }

    /// `feature.qualifiers[key]`
    pub fn quals(&self, key: &str) -> &[String] {
        self.qualifiers
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_slice())
            .unwrap_or(&[])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values from running Biopython 1.84.

    #[test]
    fn parse_splits_id_from_description_the_biopython_way() {
        // python: id='abc' name='abc' desc='abc def ghi' seq='ACGTAC'
        //         id='xyz' name='xyz' desc='xyz'         seq='TTTT'
        let recs = parse_fasta(">abc def ghi\nACGT\nAC\n>xyz\nTTTT\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "abc");
        assert_eq!(recs[0].name, "abc");
        assert_eq!(recs[0].description, "abc def ghi");
        assert_eq!(recs[0].seq, "ACGTAC");
        // a header with no description yields description == id, not ""
        assert_eq!(recs[1].description, "xyz");
    }

    #[test]
    fn write_matches_biopython_header_rule_and_wrapping() {
        // python: SeqIO.write([...], out, "fasta") produced exactly this
        let recs = vec![
            SeqRecord::new("A".repeat(130), "i1".into(), "i1".into()),
            SeqRecord::new("ACGT".into(), "i2".into(), "".into()),
            SeqRecord::new("ACGT".into(), "i3".into(), "i3 extra words".into()),
            SeqRecord::new("ACGT".into(), "i4".into(), "other text".into()),
        ];
        let mut out = Vec::new();
        write_fasta(&recs, &mut out);
        let expect = format!(
            ">i1\n{}\n{}\n{}\n>i2\nACGT\n>i3 extra words\nACGT\n>i4 other text\nACGT\n",
            "A".repeat(60),
            "A".repeat(60),
            "A".repeat(10)
        );
        assert_eq!(String::from_utf8(out).unwrap(), expect);
    }

    #[test]
    fn round_trip_preserves_sequence() {
        let recs = vec![SeqRecord::new(
            "ACGTACGT".repeat(20),
            "x".into(),
            "x".into(),
        )];
        let mut out = Vec::new();
        write_fasta(&recs, &mut out);
        let back = parse_fasta(&String::from_utf8(out).unwrap());
        assert_eq!(back[0].seq, recs[0].seq);
        assert_eq!(back[0].id, "x");
    }
}
