# Benchmark harness

Compares `panaroo-rs` against the **patched Python reference** (the same tree
`tests/parity/reference/apply.sh` builds, Tier D + Tier B applied) on **wall-clock time**
and **peak memory**, at 1 and 10 threads.

This is a benchmark, not a parity test — but it refuses to report a number without first
checking the two sides produced the same output, because a speed number for a run that
computed something different is not a speed number.

## Setup

Identical to the parity harness — it reuses the same env, the same reference tree and the
same data:

```sh
tests/parity/env/create.sh          # pinned conda env: cd-hit, mafft, biopython, ...
conda activate panaroo-parity
tests/parity/data/fetch.sh ci       # 4 real M. tuberculosis genomes (skipped if present)
tests/parity/reference/apply.sh     # the patched Python reference
```

## Running

```sh
conda activate panaroo-parity
tests/bench/run.sh                       # dataset ci, 3 repeats, -t 1 and -t 10
tests/bench/run.sh ci -n 5 -T 1,4,10     # 5 repeats, three thread counts
tests/bench/run.sh smoke -n 1            # fast sanity run
tests/bench/run.sh ci -n 3 -- -a core    # with the core-alignment stage (invokes mafft)
tests/bench/run.sh ci -P                 # python only (Rust side unavailable/mid-refactor)
```

Options: `-n` repeats, `-T` comma-separated thread counts, `-m` clean-mode, `-r`
reference tree, `-o` results dir, `-k` keep per-run output dirs, `-V` skip the
output-equivalence check, `-P` python only. Anything after `--` is passed verbatim to
**both** implementations and recorded in the results file.

## Output

Written under `tests/bench/build/run-<timestamp>/` (`tests/bench/build/latest` symlinks
the most recent), which is gitignored scratch:

| file | what |
|---|---|
| `results.tsv` | one row per run — machine readable |
| `env.txt` | CPU, cores, RAM, kernel, load, repo rev, rustc, reference stamp, resolved cd-hit/mafft/python paths + versions, exact flags |
| `logs/*.time` | raw `/usr/bin/time -v` output per run |
| `logs/*.sample.json` | process-tree memory sampler output per run |
| `logs/*.log` | stdout+stderr of each run |
| `logs/verify-t<N>.txt` | `tests/parity/canonicalise.py` comparison of the two sides |

`results.tsv` columns: `run_id dataset n_genomes impl threads rep phase exit_status
wall_s time_maxrss_kb tree_peak_rss_kb tree_peak_pss_kb tree_peak_nproc clean_mode
extra_args log cmdline`.

Re-render the summary from a saved TSV at any time:

```sh
tests/bench/summarise.py tests/bench/build/latest/results.tsv
```

## What the harness controls for

- **Same work, not just the same flags.** One dataset, one `--clean-mode`, one `-t`, and
  any extra args go to both sides verbatim; all of it is recorded in `results.tsv` and
  `env.txt`.
- **Same external binaries.** Both implementations shell out to `cd-hit` (and `mafft`
  under `-a`) and both resolve them through `PATH`. Both are launched from the harness
  shell, so they *provably* see the same `PATH`; the resolved paths and versions are
  asserted at startup and written to `env.txt`. Note that with no `-a/--alignment`,
  **mafft is never invoked by either side.**
- **Same reference code.** The Python side is `$ref` = `tests/parity/build/reference`,
  stamped with the upstream revision and the applied patches; the harness refuses to run
  against a tree with no `PARITY_REFERENCE` stamp.
- **Release build, no stale binary.** The Rust binary is behind a Cargo feature that is
  off by default, and without it cargo silently builds *nothing*. The harness deletes
  `target/release/panaroo` and rebuilds with `cargo build --release --features cli`; if
  the binary is missing afterwards it says so loudly rather than timing a leftover.
- **Output equivalence.** After the warm-up runs (which are discarded for timing anyway),
  the two output trees are compared with `tests/parity/canonicalise.py` at each thread
  count. The summary states the result. If they differ, that is reported instead of a
  clean speedup claim.
- **Caches.** One discarded warm-up run per configuration precedes the timed repeats, so
  neither side is charged for cold page cache or a cold import cache.
- **Drift.** The timed repeats are interleaved (`py, rust, py, rust, …`), not batched per
  side, so machine drift hits both equally.
- **Spread.** Every metric is min / median / max across repeats, never just a mean.
- **Python startup.** `python -c "import panaroo"` is timed separately (best of 3) and
  reported. It stays *inside* the measured runs — it is a real cost of the Python
  implementation — but it is broken out so the reader can see how much of the gap at
  small inputs is fixed overhead rather than compute.
- **Machine state.** `env.txt` and the summary record CPU model, logical core count, and
  the load average. **Nothing enforces that the box is idle** — check the load line: if
  the 1-minute load meaningfully exceeds the thread count you asked for, something else
  was running and the numbers are not trustworthy.

Each run gets a fresh, empty output directory (`build/run-*/out/<impl>-t<N>-r<rep>`),
because Panaroo writes into its output dir and would otherwise collide across runs. They
are deleted after measuring unless `-k` is given.

## Reading the memory columns — the important caveat

**`time -v maxRSS` is not the memory the run used.** GNU `/usr/bin/time` reports
`getrusage(RUSAGE_CHILDREN)`'s "Maximum resident set size", which is the peak of the
**single largest process** in the tree. It is never a sum. Both sides spawn cd-hit; the
Python side additionally forks a `multiprocessing` pool. So at `-t 10` this column can
badly understate what the machine actually had to hold. It is kept as the baseline column
because it is exact and trivially reproducible by anyone with GNU time.

The harness therefore also samples `/proc` every 100 ms over the whole process tree
(root + all live descendants) and reports:

- **`tree peak PSS`** — peak of the summed *proportional* set size (`smaps_rollup`'s
  `Pss`). PSS divides each shared page among the processes mapping it, so a forked worker
  pool is not double-counted. **This is the number to quote** for "how much memory did
  this run need".
- **`tree peak RSS`** — same sampling, summing `VmRSS`. Over-counts pages shared between
  forked children; treat it as an upper bound only.

Both sampled columns are **lower bounds**: a spike shorter than the 100 ms sampling
interval is missed. `procs` is the largest number of live processes seen in one sample.

One asymmetry is real rather than an artefact: Rust parallelises with threads inside one
process, Python with forked processes. At `-t 10` the Rust `maxRSS` is essentially the
whole run's footprint, while the Python `maxRSS` is one worker's share. Compare the
**PSS** column across implementations; compare `maxRSS` only within one implementation.
