//! Dump a GFF3 file's parse in a canonical form, for differential testing against gffutils.
//!
//! Pairs with `tests/parity/gff/dump_gffutils.py`, which emits the identical format from
//! the real library. Any difference between the two is a bug in `support::gff`.
//!
//!     cargo run --release --example dump_gff -- FILE.gff > rs.tsv
//!
//! Input handling mirrors `prokka::get_gene_sequences`: strip commas, split off `##FASTA`,
//! drop `##sequence-region` lines.
use panaroo::support::gff::GffDb;

fn main() {
    let path = std::env::args().nth(1).expect("usage: dump_gff FILE.gff");
    let text = std::fs::read_to_string(&path).expect("read");
    let text = text.replace(',', "");
    let ann = text.split("##FASTA").next().unwrap().to_string();
    let ann: String = ann
        .lines()
        .filter(|l| !l.contains("##sequence-region"))
        .collect::<Vec<_>>()
        .join("\n");

    let db = GffDb::create_db(&ann).expect("parse");
    for e in db.all_features() {
        let attrs: Vec<String> = e
            .attributes
            .items()
            .map(|(k, vs)| format!("{}={}", k, vs.join("\u{1}")))
            .collect();
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.id,
            e.seqid,
            e.source,
            e.featuretype,
            e.start,
            e.stop,
            e.score,
            e.strand,
            e.frame,
            attrs.join("\u{2}")
        );
    }
}
