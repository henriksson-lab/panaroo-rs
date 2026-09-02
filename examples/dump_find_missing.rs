//! Emit the same values as `tests/parity/find_missing/cases.py`, for differential testing.
//!
//!     cargo run --release --example dump_find_missing -- CASES.tsv
//! where CASES.tsv is the Python script's output; the input columns are replayed and the
//! result columns recomputed.
use panaroo::find_missing::{search_dna, translate_to_match};

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: dump_find_missing CASES.tsv");
    for line in std::fs::read_to_string(&path).expect("read").lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "search_dna" => {
                let (name, q, db, pm, pid) = (f[1], f[2], f[3], f[4], f[5]);
                let (seq, loc) = search_dna(db, q, pm.parse().unwrap(), pid.parse().unwrap(), true);
                let locs: Vec<String> = loc.iter().map(|x| x.to_string()).collect();
                println!(
                    "search_dna\t{name}\t{q}\t{db}\t{pm}\t{pid}\t{seq}\t{}",
                    locs.join(",")
                );
            }
            "translate_to_match" => {
                let (name, hit, prot) = (f[1], f[2], f[3]);
                println!(
                    "translate_to_match\t{name}\t{hit}\t{prot}\t{}",
                    translate_to_match(hit, prot)
                );
            }
            other => panic!("unknown case kind: {other}"),
        }
    }
}
