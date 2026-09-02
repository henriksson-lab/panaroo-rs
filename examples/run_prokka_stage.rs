//! Run only `prokka::process_prokka_input`, for the Phase 3 parity checkpoint.
//!
//! Pairs with `tests/parity/prokka/run_reference.py`, which runs the reference Python's
//! `process_prokka_input` on the same input. The two output directories must be
//! byte-identical.
//!
//!     cargo run --release --example run_prokka_stage -- INPUT_LIST OUT_DIR [N_CPU]
use panaroo::prokka::process_prokka_input;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let input_list = a
        .get(1)
        .expect("usage: run_prokka_stage INPUT_LIST OUT_DIR [N_CPU]");
    let out_dir = a
        .get(2)
        .expect("usage: run_prokka_stage INPUT_LIST OUT_DIR [N_CPU]");
    let n_cpu: i64 = a.get(3).map(|s| s.parse().unwrap()).unwrap_or(1);

    let files: Vec<String> = std::fs::read_to_string(input_list)
        .expect("read input list")
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();

    std::fs::create_dir_all(out_dir).expect("mkdir out");
    let out_dir = if out_dir.ends_with('/') {
        out_dir.clone()
    } else {
        format!("{out_dir}/")
    };

    // filter_invalid defaults to False; --clean-mode does not affect this stage.
    process_prokka_input(&files, &out_dir, false, true, n_cpu, 11);
}
