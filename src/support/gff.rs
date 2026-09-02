//! GFF3 reader — a behavioural stand-in for `gffutils`.
//!
//! # Provenance and licence
//!
//! This module has **no Python counterpart in Panaroo**. It replaces the third-party
//! library Panaroo calls:
//!
//! > **gffutils** 0.13 — <https://github.com/daler/gffutils>
//! > The MIT License (MIT), Copyright (c) 2013 Ryan Dale
//!
//! **No gffutils source code was copied.** The observable behaviour reproduced here was
//! derived empirically, by running gffutils 0.13 against probe inputs and recording what it
//! returned; each such observation is cited inline below as `OBSERVED`. Panaroo itself
//! (<https://github.com/gtonkinhill/panaroo>, MIT, Copyright (c) 2019 Gerry Tonkin-Hill) is
//! the consumer whose usage defines the required surface.
//!
//! Format reference: GFF3 specification, Sequence Ontology,
//! <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>.
//!
//! See `NOTICE.md` at the repository root.
//!
//! # Why a reader rather than a transcription
//!
//! Upstream calls
//! `gff.create_db(text, dbfn=":memory:", force=True, keep_order=True, from_string=True,
//! merge_strategy="create_unique")`, which builds an in-memory SQLite database per genome
//! and then does a handful of point lookups against it. There is no Rust gffutils, so
//! "translate it as-is" is not on the menu — the only choice is *which* faithful
//! reimplementation to write, and a direct parser is both simpler and faster than emulating
//! a SQLite round trip. This is a translation decision, not an early optimisation; see
//! `PORTING_PLAN.md` §9 item 10.
//!
//! # Required surface
//!
//! Panaroo touches only this much of a gffutils `FeatureDB` / `Feature`:
//!
//! | usage | site |
//! |---|---|
//! | `db.all_features(featuretype=())` | `prokka.py::get_gene_sequences` |
//! | `db[feature_id]` | `find_missing.py::search_gff` |
//! | `.id`, `.seqid`, `.featuretype`, `.start`, `.stop`, `.strand`, `.frame` | both |
//! | `.attributes["gene" \| "name" \| "product"]`, raising `KeyError` when absent | `get_gene_sequences` |
//! | `gene[0]` — positional indexing, column 1 (seqid) | `search_gff` |
//! | `gene.start`, `gene.end` — `end` is an alias of `stop` | `search_gff` |
//!
//! # Behaviour reproduced (all `OBSERVED` against gffutils 0.13)
//!
//! 1. **File order is preserved** under `keep_order=True`: `all_features` yields features in
//!    the order they appear in the text. This is load-bearing — `get_gene_sequences` numbers
//!    genes `{file}_{scaffold_index}_{gene_index}`, and `scaffold_index` comes from the
//!    order contigs are first seen while walking features. Those IDs reach every output
//!    file. See `PORTING_PLAN.md` §6.6.
//! 2. **`id` is the `ID` attribute** when present.
//! 3. **`merge_strategy="create_unique"`** appends `_1`, `_2`, … to a repeated `ID`, in
//!    order of appearance. The `attributes["ID"]` value keeps the *original*, so
//!    `entry.id != entry.attributes["ID"][0]` for a duplicate. Lookup by the uniquified id
//!    works.
//! 4. **A feature with no `ID`** gets `{featuretype}_{n}`, where `n` counts from 1 **per
//!    featuretype** (`CDS_1`, `tRNA_1`, `CDS_2`).
//! 5. **Columns 1–8 are *not* URL-decoded** — `seqid` `c%2C1` stays `c%2C1`.
//! 6. **Column 9 is split before it is decoded.** The `;` and `,` separators are found in the
//!    raw text, and only then is each value percent-decoded. So `product=foo%2Cbar` is one
//!    value `foo,bar`, while `product=foo,bar` is two values `["foo", "bar"]`. Getting this
//!    order backwards changes the annotation strings Panaroo writes.
//! 7. **Attribute key order is preserved**, and repeated keys accumulate into one list.
//! 8. **Surrounding double quotes are stripped** from a value.
//! 9. **A trailing `;`** produces no empty attribute.
//! 10. `start` and `stop` are 1-based inclusive integers; `score`, `strand` and `frame` are
//!     the raw column strings (`frame` keeps `"."`).
//!
//! Percent-escapes are not hypothetical: the four `ci` test genomes contain 8565 `%2C`
//! sequences, all inside `product=` values that Panaroo joins into the `annotation` and
//! `description` node attributes.
//!
//! # Known deviations
//!
//! gffutils *sniffs* an attribute dialect (separator, quoting, trailing semicolon) from the
//! first lines of a file and then applies it uniformly. This module instead implements the
//! GFF3 dialect directly and handles quoting and trailing separators per field. For the
//! inputs Panaroo accepts — Prokka GFF3, and the output of Panaroo's own
//! `convert_refseq_to_prokka_gff.py` — the two agree; for a file using `key value`
//! GTF-style attributes they would not. [`GffDb::create_db`] returns an error on a line it
//! cannot parse as GFF3 rather than guessing.

use super::pydict::PyDict;
use std::fmt;

/// One GFF3 feature. A gffutils `Feature`, reduced to the fields Panaroo reads.
#[derive(Debug, Clone, Default)]
pub struct GffEntry {
    /// `entry.id` — the `ID` attribute, uniquified per behaviour 3, or `{featuretype}_{n}`
    /// per behaviour 4. **Not** necessarily equal to `attributes["ID"][0]`.
    pub id: String,
    /// Column 1. `entry.seqid`, also `gene[0]` in `find_missing::search_gff`.
    pub seqid: String,
    /// Column 2. `entry.source`
    pub source: String,
    /// Column 3. `entry.featuretype`
    pub featuretype: String,
    /// Column 4, 1-based inclusive. `entry.start`
    pub start: i64,
    /// Column 5, 1-based inclusive. `entry.stop`, aliased as `entry.end`.
    pub stop: i64,
    /// Column 6, raw. `entry.score`
    pub score: String,
    /// Column 7. `entry.strand`
    pub strand: String,
    /// Column 8, raw — keeps `"."`. `entry.frame`
    pub frame: String,
    /// Column 9, parsed. Values are lists, key order preserved.
    pub attributes: PyDict<String, Vec<String>>,
}

impl GffEntry {
    /// `entry.end` — gffutils exposes `stop` under both names.
    pub fn end(&self) -> i64 {
        self.stop
    }

    /// `gene[i]` — positional access to the nine columns. Panaroo only ever uses `gene[0]`.
    pub fn column(&self, i: usize) -> String {
        match i {
            0 => self.seqid.clone(),
            1 => self.source.clone(),
            2 => self.featuretype.clone(),
            3 => self.start.to_string(),
            4 => self.stop.to_string(),
            5 => self.score.clone(),
            6 => self.strand.clone(),
            7 => self.frame.clone(),
            8 => format_attributes(&self.attributes),
            _ => panic!("gff: column index {i} out of range"),
        }
    }

    /// `entry.attributes[key]`, which raises `KeyError` when absent.
    ///
    /// Returns `None` for the missing case so call sites can transcribe Python's
    /// `try/except KeyError` as a match — see `prokka::attr_first_or_empty`.
    pub fn attribute(&self, key: &str) -> Option<&Vec<String>> {
        self.attributes.get(&key.to_string())
    }
}

fn format_attributes(attrs: &PyDict<String, Vec<String>>) -> String {
    let mut parts = Vec::new();
    for (k, vs) in attrs.items() {
        parts.push(format!("{}={}", k, vs.join(",")));
    }
    parts.join(";")
}

/// Failure to parse a GFF3 line. Panaroo surfaces these as `RuntimeError`.
#[derive(Debug, Clone)]
pub struct GffError {
    pub line_number: usize,
    pub message: String,
    pub line: String,
}

impl fmt::Display for GffError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "GFF3 parse error on line {}: {} -- {:?}",
            self.line_number, self.message, self.line
        )
    }
}

impl std::error::Error for GffError {}

/// The result of `gff.create_db(...)`.
#[derive(Debug, Clone, Default)]
pub struct GffDb {
    entries: Vec<GffEntry>,
    /// `entry.id` -> index into `entries`.
    by_id: PyDict<String, usize>,
}

impl GffDb {
    /// `gff.create_db(text, dbfn=":memory:", force=True, keep_order=True, from_string=True,
    /// merge_strategy="create_unique")`
    ///
    /// `text` is the annotation block only: Panaroo splits off `##FASTA` and drops
    /// `##sequence-region` lines (`prokka::clean_gff_string`) before calling this.
    pub fn create_db(text: &str) -> Result<GffDb, GffError> {
        let mut entries: Vec<GffEntry> = Vec::new();
        let mut by_id: PyDict<String, usize> = PyDict::new();
        // Behaviour 4: the auto-id counter runs per featuretype.
        let mut auto_counter: PyDict<String, usize> = PyDict::new();

        for (i, raw) in text.lines().enumerate() {
            let line_number = i + 1;
            let line = raw.strip_suffix('\r').unwrap_or(raw);
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                return Err(GffError {
                    line_number,
                    message: format!("expected 9 tab-separated columns, found {}", cols.len()),
                    line: line.to_string(),
                });
            }

            let parse_coord = |s: &str, what: &str| -> Result<i64, GffError> {
                s.trim().parse::<i64>().map_err(|_| GffError {
                    line_number,
                    message: format!("{what} is not an integer: {s:?}"),
                    line: line.to_string(),
                })
            };

            // Behaviour 5: columns 1-8 are kept verbatim, not percent-decoded.
            let mut entry = GffEntry {
                id: String::new(),
                seqid: cols[0].to_string(),
                source: cols[1].to_string(),
                featuretype: cols[2].to_string(),
                start: parse_coord(cols[3], "start")?,
                stop: parse_coord(cols[4], "stop")?,
                score: cols[5].to_string(),
                strand: cols[6].to_string(),
                frame: cols[7].to_string(),
                // Column 9 may itself contain tabs only if the file is malformed; rejoin
                // any extra splits so we do not silently truncate an attribute.
                attributes: parse_attributes(&cols[8..].join("\t")),
            };

            // Behaviour 2 and 4.
            let base_id = match entry.attributes.get(&"ID".to_string()) {
                Some(vs) if !vs.is_empty() => vs[0].clone(),
                _ => {
                    let n = auto_counter.get(&entry.featuretype).copied().unwrap_or(0) + 1;
                    auto_counter.insert(entry.featuretype.clone(), n);
                    format!("{}_{}", entry.featuretype, n)
                }
            };

            // Behaviour 3: merge_strategy="create_unique". `attributes["ID"]` is left alone.
            let mut unique = base_id.clone();
            if by_id.contains_key(&unique) {
                let mut suffix = 1usize;
                loop {
                    let candidate = format!("{base_id}_{suffix}");
                    if !by_id.contains_key(&candidate) {
                        unique = candidate;
                        break;
                    }
                    suffix += 1;
                }
            }
            entry.id = unique.clone();

            by_id.insert(unique, entries.len());
            entries.push(entry);
        }

        Ok(GffDb { entries, by_id })
    }

    /// `db.all_features(featuretype=())` — file order (behaviour 1).
    pub fn all_features(&self) -> &[GffEntry] {
        &self.entries
    }

    /// `db[feature_id]`. gffutils raises `FeatureNotFoundError`; Panaroo never guards, so a
    /// miss is a bug and panicking matches the Python's observable behaviour.
    pub fn get(&self, id: &str) -> &GffEntry {
        match self.by_id.get(&id.to_string()) {
            Some(&i) => &self.entries[i],
            None => panic!("gff: FeatureNotFoundError: {id}"),
        }
    }

    /// Non-panicking lookup, for callers that want to test presence.
    pub fn try_get(&self, id: &str) -> Option<&GffEntry> {
        self.by_id.get(&id.to_string()).map(|&i| &self.entries[i])
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Parse column 9.
///
/// Behaviour 6 is the subtle one: split on the raw `;` and `,`, *then* percent-decode each
/// value. Decoding first would turn `product=foo%2Cbar` into two values instead of one.
/// Behaviours 7, 8 and 9 also live here.
fn parse_attributes(field: &str) -> PyDict<String, Vec<String>> {
    let mut attrs: PyDict<String, Vec<String>> = PyDict::new();

    for part in field.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue; // behaviour 9: trailing ';' yields nothing
        }
        let (key, raw_value) = match part.split_once('=') {
            Some((k, v)) => (k.trim(), v),
            // A bare token with no '=' is not GFF3; gffutils keeps it as a valueless key.
            None => (part, ""),
        };

        let values: Vec<String> = if raw_value.is_empty() {
            Vec::new()
        } else {
            raw_value
                .split(',') // behaviour 6: split raw ...
                .map(|v| percent_decode(strip_quotes(v))) // ... then decode (behaviour 8 first)
                .collect()
        };

        // Behaviour 7: repeated keys accumulate; first occurrence fixes the key's position.
        match attrs.get_mut(&key.to_string()) {
            Some(existing) => existing.extend(values),
            None => {
                attrs.insert(key.to_string(), values);
            }
        }
    }

    attrs
}

/// Behaviour 8: a value wrapped in double quotes has them removed.
fn strip_quotes(v: &str) -> &str {
    let v = v.trim();
    if v.len() >= 2 && v.starts_with('"') && v.ends_with('"') {
        &v[1..v.len() - 1]
    } else {
        v
    }
}

/// Percent-decode a GFF3 attribute value.
///
/// Bytes are decoded then reassembled as UTF-8, so a multi-byte character written as a run
/// of escapes (`%C3%A9`) round-trips. An escape that is not two hex digits is left verbatim,
/// which is what a permissive reader must do for files that use a bare `%`.
fn percent_decode(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    let bytes = s.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            let hi = (bytes[i + 1] as char).to_digit(16);
            let lo = (bytes[i + 2] as char).to_digit(16);
            if let (Some(h), Some(l)) = (hi, lo) {
                out.push((h * 16 + l) as u8);
                i += 3;
                continue;
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every expectation here is an `OBSERVED` value: it was produced by running
    /// gffutils 0.13 on the same input, not by reading its source.
    fn db(text: &str) -> GffDb {
        GffDb::create_db(text).expect("parse")
    }

    #[test]
    fn file_order_is_preserved() {
        let t = "##gff-version 3\n\
                 c1\tp\tCDS\t40\t49\t.\t+\t0\tID=c\n\
                 c1\tp\tCDS\t1\t9\t.\t+\t0\tID=a\n\
                 c1\tp\tCDS\t20\t29\t.\t+\t0\tID=b\n";
        let d = db(t);
        let ids: Vec<&str> = d.all_features().iter().map(|e| e.id.as_str()).collect();
        assert_eq!(ids, ["c", "a", "b"]);
    }

    #[test]
    fn duplicate_ids_get_create_unique_suffixes() {
        let t = "##gff-version 3\n\
                 c1\tp\tCDS\t1\t9\t.\t+\t0\tID=dup\n\
                 c1\tp\tCDS\t20\t29\t.\t+\t0\tID=dup\n\
                 c1\tp\tCDS\t40\t49\t.\t+\t0\tID=dup\n";
        let d = db(t);
        let ids: Vec<&str> = d.all_features().iter().map(|e| e.id.as_str()).collect();
        assert_eq!(ids, ["dup", "dup_1", "dup_2"]);
        // The ID attribute keeps the original on every one of them.
        for e in d.all_features() {
            assert_eq!(e.attribute("ID").unwrap(), &vec!["dup".to_string()]);
        }
        assert_eq!(d.get("dup").start, 1);
        assert_eq!(d.get("dup_1").start, 20);
        assert_eq!(d.get("dup_2").start, 40);
    }

    #[test]
    fn missing_id_counts_per_featuretype() {
        let t = "##gff-version 3\n\
                 c1\tp\tCDS\t1\t9\t.\t+\t0\tgene=a\n\
                 c1\tp\ttRNA\t20\t29\t.\t+\t0\tgene=b\n\
                 c1\tp\tCDS\t40\t49\t.\t+\t0\tgene=c\n";
        let d = db(t);
        let ids: Vec<&str> = d.all_features().iter().map(|e| e.id.as_str()).collect();
        assert_eq!(ids, ["CDS_1", "tRNA_1", "CDS_2"]);
    }

    #[test]
    fn columns_one_to_eight_are_not_decoded() {
        let t = "##gff-version 3\nc%2C1\tsrc%3Bx\tCDS\t1\t9\t.\t+\t0\tID=a\n";
        let d = db(t);
        let e = &d.all_features()[0];
        assert_eq!(e.seqid, "c%2C1");
        assert_eq!(e.source, "src%3Bx");
        assert_eq!(e.column(0), "c%2C1");
    }

    #[test]
    fn split_happens_before_decode() {
        // The distinction that matters: %2C is one value, a literal comma is two.
        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=x;note=foo%2Cbar\n";
        assert_eq!(
            db(t).all_features()[0].attribute("note").unwrap(),
            &["foo,bar"]
        );

        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=x;note=foo,bar\n";
        assert_eq!(
            db(t).all_features()[0].attribute("note").unwrap(),
            &["foo", "bar"]
        );

        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=x;note=a%3Bb\n";
        assert_eq!(db(t).all_features()[0].attribute("note").unwrap(), &["a;b"]);
    }

    #[test]
    fn attribute_key_order_is_preserved() {
        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tzz=1;ID=a;aa=2;mm=3\n";
        let d = db(t);
        let keys: Vec<&String> = d.all_features()[0].attributes.keys().collect();
        assert_eq!(keys, ["zz", "ID", "aa", "mm"]);
    }

    #[test]
    fn quotes_stripped_and_trailing_semicolon_ignored() {
        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=a;note=\"hi there\";\n";
        let d = db(t);
        let e = &d.all_features()[0];
        assert_eq!(e.attribute("note").unwrap(), &["hi there"]);
        let keys: Vec<&String> = e.attributes.keys().collect();
        assert_eq!(keys, ["ID", "note"]);
    }

    #[test]
    fn repeated_keys_accumulate() {
        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=a;product=x;product=y\n";
        assert_eq!(
            db(t).all_features()[0].attribute("product").unwrap(),
            &["x", "y"]
        );
    }

    #[test]
    fn coordinates_and_raw_columns() {
        let t = "##gff-version 3\nc1\tp\tCDS\t40\t49\t.\t-\t.\tID=a\n";
        let d = db(t);
        let e = &d.all_features()[0];
        assert_eq!((e.start, e.stop, e.end()), (40, 49, 49));
        assert_eq!(
            (e.score.as_str(), e.strand.as_str(), e.frame.as_str()),
            (".", "-", ".")
        );
    }

    #[test]
    fn missing_attribute_is_none_not_panic() {
        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=a\n";
        assert!(db(t).all_features()[0].attribute("gene").is_none());
    }

    #[test]
    fn short_line_is_an_error_not_a_guess() {
        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\n";
        let err = GffDb::create_db(t).unwrap_err();
        assert_eq!(err.line_number, 2);
        assert!(err.message.contains("9 tab-separated columns"));
    }

    #[test]
    fn multibyte_percent_escapes_round_trip() {
        let t = "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=a;note=caf%C3%A9\n";
        assert_eq!(
            db(t).all_features()[0].attribute("note").unwrap(),
            &["café"]
        );
    }
}
