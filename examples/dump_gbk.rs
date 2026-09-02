//! Convert a GenBank file to GFF3 via `biocode_convert::convert_gbk_gff3`.
//!
//!     cargo run --release --example dump_gbk -- IN.gbff OUT.gff
//! Compare against the reference with `tests/parity/gbk/check.sh`.
fn main() {
    let a: Vec<String> = std::env::args().collect();
    panaroo::biocode_convert::convert_gbk_gff3(&a[1], &a[2], true);
}
