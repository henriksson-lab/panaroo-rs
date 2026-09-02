//! GenBank flat-file reader.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [Biopython](https://biopython.org/) 1.84's `SeqIO.parse(handle, "genbank")`,
//! Biopython License Agreement / BSD 3-Clause. No Biopython source was copied; the fields
//! and location semantics below were established by running the library and by comparing
//! the resulting GFF3 against the reference. See `NOTICE.md`.
//!
//! # Scope
//!
//! Only what `biocode_convert::convert_gbk_gff3` reads:
//!
//!   - `record.name` — the first token of the `LOCUS` line, used as the GFF3 seqid
//!   - `record.seq` — the `ORIGIN` block, used for the appended `##FASTA`
//!   - `record.features` — each with `type`, a location, and qualifiers
//!   - per feature: `location.start` / `location.end` (0-based half-open, as Biopython
//!     reports them), `location.strand` (`+1` / `-1`), and `location.parts` for a
//!     `join(...)`
//!
//! Everything else in a GenBank record (REFERENCE, COMMENT, the header annotations) is
//! skipped, because nothing downstream reads it.
//!
//! # Location syntax handled
//!
//! `123..456`, `complement(...)`, `join(...)`, `order(...)`, and the `<`/`>` partial
//! markers, which Biopython strips when reporting the numeric bounds. A location with no
//! resolvable strand is reported as `strand: 0`, which `convert_gbk_gff3` turns into the
//! upstream `Exception("ERROR: unstranded feature encountered")`.

use super::seqio::{GenBankFeature, GenBankRecord};

/// `SeqIO.parse(open(path), "genbank")`
pub fn parse_genbank_file(path: &str) -> Vec<GenBankRecord> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("could not read {path}: {e}"));
    parse_genbank(&text)
}

/// `SeqIO.parse(handle, "genbank")`
pub fn parse_genbank(text: &str) -> Vec<GenBankRecord> {
    let mut records = Vec::new();
    let mut cur: Option<GenBankRecord> = None;
    let mut section = Section::Header;

    // A feature is accumulated across its continuation lines before being parsed.
    let mut feat_type = String::new();
    let mut feat_loc = String::new();
    let mut feat_quals: Vec<String> = Vec::new();

    let flush_feature = |cur: &mut Option<GenBankRecord>,
                         feat_type: &mut String,
                         feat_loc: &mut String,
                         feat_quals: &mut Vec<String>| {
        if feat_type.is_empty() {
            return;
        }
        if let Some(rec) = cur.as_mut() {
            let (start, end, strand, parts) = parse_location(feat_loc);
            rec.features.push(GenBankFeature {
                feature_type: std::mem::take(feat_type),
                start,
                end,
                strand,
                parts,
                qualifiers: parse_qualifiers(feat_quals),
            });
        }
        feat_type.clear();
        feat_loc.clear();
        feat_quals.clear();
    };

    for raw in text.lines() {
        let line = raw.strip_suffix('\r').unwrap_or(raw);

        if line.starts_with("LOCUS") {
            flush_feature(&mut cur, &mut feat_type, &mut feat_loc, &mut feat_quals);
            if let Some(r) = cur.take() {
                records.push(r);
            }
            let name = line.split_whitespace().nth(1).unwrap_or("").to_string();
            cur = Some(GenBankRecord {
                id: name.clone(),
                name,
                ..Default::default()
            });
            section = Section::Header;
            continue;
        }
        if line.starts_with("//") {
            flush_feature(&mut cur, &mut feat_type, &mut feat_loc, &mut feat_quals);
            if let Some(r) = cur.take() {
                records.push(r);
            }
            section = Section::Header;
            continue;
        }
        if line.starts_with("FEATURES") {
            section = Section::Features;
            continue;
        }
        if line.starts_with("ORIGIN") {
            flush_feature(&mut cur, &mut feat_type, &mut feat_loc, &mut feat_quals);
            section = Section::Origin;
            continue;
        }
        // any other column-0 keyword ends the features block
        if !line.starts_with(' ') && !line.is_empty() {
            if matches!(section, Section::Features) {
                flush_feature(&mut cur, &mut feat_type, &mut feat_loc, &mut feat_quals);
            }
            section = Section::Header;
            continue;
        }

        match section {
            Section::Features => {
                // A feature key sits at column 5 (0-based), its location at column 21.
                let is_new_feature = line.len() > 5
                    && line.as_bytes()[..5].iter().all(|&c| c == b' ')
                    && line.as_bytes()[5] != b' ';
                if is_new_feature {
                    flush_feature(&mut cur, &mut feat_type, &mut feat_loc, &mut feat_quals);
                    let rest = &line[5..];
                    let mut it = rest.splitn(2, char::is_whitespace);
                    feat_type = it.next().unwrap_or("").to_string();
                    feat_loc = it.next().unwrap_or("").trim().to_string();
                } else {
                    let t = line.trim();
                    if t.starts_with('/') {
                        feat_quals.push(t.to_string());
                    } else if let Some(last) = feat_quals.last_mut() {
                        // Continuation of a qualifier value. Biopython joins with a single
                        // space **unconditionally** -- including when the previous line ends
                        // in a hyphen, which is tempting to special-case but wrong.
                        // OBSERVED: `product_name=2-amino-4-hydroxy-6- hydroxymethyl...`
                        // keeps the space in the reference output.
                        last.push(' ');
                        last.push_str(t);
                    } else {
                        // continuation of the location
                        feat_loc.push_str(t);
                    }
                }
            }
            Section::Origin => {
                if let Some(rec) = cur.as_mut() {
                    for tok in line.split_whitespace() {
                        if tok.chars().all(|c| c.is_ascii_digit()) {
                            continue; // the column offset
                        }
                        rec.seq.push_str(&tok.to_ascii_uppercase());
                    }
                }
            }
            Section::Header => {}
        }
    }

    flush_feature(&mut cur, &mut feat_type, &mut feat_loc, &mut feat_quals);
    if let Some(r) = cur {
        records.push(r);
    }
    records
}

enum Section {
    Header,
    Features,
    Origin,
}

/// Parse a GenBank location into `(start, end, strand, parts)`.
///
/// `start`/`end` are 0-based half-open, matching `feature.location.start` / `.end` in
/// Biopython. `parts` holds each `join(...)` segment in the same convention; a simple
/// location yields one part.
fn parse_location(loc: &str) -> (i64, i64, i32, Vec<(i64, i64)>) {
    let loc = loc.trim();
    let (body, strand) = if let Some(inner) = strip_call(loc, "complement") {
        (inner.to_string(), -1)
    } else {
        (loc.to_string(), 1)
    };

    let inner = strip_call(&body, "join")
        .or_else(|| strip_call(&body, "order"))
        .map(|s| s.to_string())
        .unwrap_or(body);

    let mut parts: Vec<(i64, i64)> = Vec::new();
    for seg in split_top_level(&inner) {
        // a nested complement inside a join keeps the outer strand for our purposes --
        // Panaroo's input never mixes strands within one feature
        let seg = strip_call(seg.trim(), "complement")
            .unwrap_or(seg.trim())
            .to_string();
        if let Some((a, b)) = seg.split_once("..") {
            let a = a.trim_start_matches(['<', '>']).trim();
            let b = b.trim_start_matches(['<', '>']).trim();
            if let (Ok(a), Ok(b)) = (a.parse::<i64>(), b.parse::<i64>()) {
                parts.push((a - 1, b)); // 1-based inclusive -> 0-based half-open
            }
        } else {
            let s = seg.trim_start_matches(['<', '>']).trim();
            if let Ok(p) = s.parse::<i64>() {
                parts.push((p - 1, p));
            }
        }
    }

    if parts.is_empty() {
        return (0, 0, 0, parts);
    }
    let start = parts.iter().map(|p| p.0).min().unwrap();
    let end = parts.iter().map(|p| p.1).max().unwrap();

    // Biopython reports parts of a complement location in reverse (biological) order.
    if strand == -1 {
        parts.reverse();
    }
    (start, end, strand, parts)
}

/// `name(...)` -> the inner text, if the whole string is that call.
fn strip_call<'a>(s: &'a str, name: &str) -> Option<&'a str> {
    let prefix = format!("{name}(");
    if s.starts_with(&prefix) && s.ends_with(')') {
        Some(&s[prefix.len()..s.len() - 1])
    } else {
        None
    }
}

/// Split on commas that are not inside parentheses.
fn split_top_level(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut start = 0usize;
    for (i, c) in s.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            ',' if depth == 0 => {
                out.push(&s[start..i]);
                start = i + c.len_utf8();
            }
            _ => {}
        }
    }
    out.push(&s[start..]);
    out
}

/// Parse `/key="value"` and `/flag` lines into Biopython's `qualifiers` shape.
fn parse_qualifiers(lines: &[String]) -> Vec<(String, Vec<String>)> {
    let mut out: Vec<(String, Vec<String>)> = Vec::new();
    for l in lines {
        let l = l.trim_start_matches('/');
        let (key, value) = match l.split_once('=') {
            Some((k, v)) => {
                let v = v.trim();
                let v = v.strip_prefix('"').unwrap_or(v);
                let v = v.strip_suffix('"').unwrap_or(v);
                (k.to_string(), v.to_string())
            }
            None => (l.to_string(), String::new()),
        };
        match out.iter_mut().find(|(k, _)| *k == key) {
            Some((_, vs)) => vs.push(value),
            None => out.push((key, vec![value])),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_locations_like_biopython() {
        // Biopython reports 0-based half-open bounds and +1/-1 strand.
        assert_eq!(parse_location("1..1323"), (0, 1323, 1, vec![(0, 1323)]));
        assert_eq!(parse_location("complement(5..9)"), (4, 9, -1, vec![(4, 9)]));
        assert_eq!(parse_location("<1..100"), (0, 100, 1, vec![(0, 100)]));
        assert_eq!(parse_location("100..>200"), (99, 200, 1, vec![(99, 200)]));
        let (s, e, st, parts) = parse_location("join(1..10,20..30)");
        assert_eq!((s, e, st), (0, 30, 1));
        assert_eq!(parts, vec![(0, 10), (19, 30)]);
        // a complement join reports its parts in biological order
        let (_, _, st, parts) = parse_location("complement(join(1..10,20..30))");
        assert_eq!(st, -1);
        assert_eq!(parts, vec![(19, 30), (0, 10)]);
    }

    #[test]
    fn parses_a_minimal_record() {
        // Built line by line on purpose: a `\`-continuation in a Rust string literal
        // strips the leading whitespace of the next line, and GenBank is column-sensitive.
        let gb = [
            "LOCUS       NZ_TEST     30 bp    DNA     circular CON 28-MAY-2025",
            "DEFINITION  Test.",
            "FEATURES             Location/Qualifiers",
            "     source          1..30",
            "                     /organism=\"Testus\"",
            "     gene            1..9",
            "                     /locus_tag=\"T_0001\"",
            "     CDS             1..9",
            "                     /locus_tag=\"T_0001\"",
            "                     /gene=\"abc\"",
            "                     /product=\"a test",
            "                     protein\"",
            "     gene            complement(11..19)",
            "                     /locus_tag=\"T_0002\"",
            "ORIGIN",
            "        1 atgaaattta tgaaatttat gaaatttatg",
            "//",
        ]
        .join("\n");
        let gb = gb.as_str();
        let recs = parse_genbank(gb);
        assert_eq!(recs.len(), 1);
        let r = &recs[0];
        assert_eq!(r.name, "NZ_TEST");
        assert_eq!(r.seq, "ATGAAATTTATGAAATTTATGAAATTTATG");
        let types: Vec<&str> = r.features.iter().map(|f| f.feature_type.as_str()).collect();
        assert_eq!(types, ["source", "gene", "CDS", "gene"]);
        let cds = &r.features[2];
        assert_eq!((cds.start, cds.end, cds.strand), (0, 9, 1));
        let q = |k: &str| -> Option<&Vec<String>> {
            cds.qualifiers.iter().find(|(a, _)| a == k).map(|(_, v)| v)
        };
        assert_eq!(q("locus_tag").unwrap()[0], "T_0001");
        assert_eq!(q("gene").unwrap()[0], "abc");
        // a qualifier value wrapped across lines is rejoined with a single space
        assert_eq!(q("product").unwrap()[0], "a test protein");
    }
}
