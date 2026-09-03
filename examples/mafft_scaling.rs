//! Does the embedded MAFFT path scale across threads?
//!
//! Panaroo aligns ~5,100 genes by running `parallel_map` over N worker threads, each
//! invoking MAFFT with `--thread 1`. The parallelism is ours; MAFFT itself is serial. This
//! measures whether that actually scales in-process.
//!
//!   cargo run --release --features mafft-embedded --example mafft_scaling -- <gene.fa> <reps>
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let path = a.get(1).expect("usage: mafft_scaling <gene.fa> [reps]");
    let reps: usize = a.get(2).map(|s| s.parse().unwrap()).unwrap_or(64);
    let cmd = format!("mafft --auto --adjustdirection --thread 1 --nuc {path}");

    for workers in [1usize, 2, 4, 8, 20] {
        let jobs: Vec<usize> = (0..reps).collect();
        let t = std::time::Instant::now();
        let _ = panaroo::support::parallel::parallel_map(workers as i64, jobs, |_| {
            panaroo::mafft_embedded::run_mafft(&cmd).len()
        });
        let el = t.elapsed().as_secs_f64();
        println!(
            "  workers={workers:<3} {reps} alignments in {el:6.2} s  ({:6.1} aln/s, speedup {:.2}x)",
            reps as f64 / el,
            0.0
        );
    }
}
