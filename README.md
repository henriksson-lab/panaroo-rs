# panaroo-rs

A faithful *mostly*-Rust translation of [Panaroo](https://github.com/gtonkinhill/panaroo), the
prokaryotic pangenome pipeline.

Note that this crate is not fully Rust yet. Translation of upstream dependencies is needed first. The aim of this translation is to improve speed over the original Panaroo code

Also note that the Panaroo code output is not reproduced due to translation challenges and a suspected upstream bug

**not yet tested enough**

* 2026-09-02: Initial translation


## Installing

The command-line tool is behind the **`cli` feature, which is off by default**, because
the crate's primary product is the library — the translated pipeline functions, which you
can drive directly. A library consumer should not have to build `clap`.

To get the `panaroo` executable, ask for the feature explicitly:

```sh
cargo install panaroo-rs --features cli
```

or from a checkout:

```sh
cargo build --release --features cli   # vendored edlib is compiled by build.rs; no cmake, no libclang
```

To use it as a library instead, the default feature set is what you want:

```toml
[dependencies]
panaroo-rs = "0.1"        # library only, no clap
```

### In-process cd-hit and MAFFT (optional)

By default both external tools are invoked exactly as upstream Panaroo invokes them: as
subprocesses, from `PATH`. Two further features run them **in-process** instead, via Rust
translations of each:

| feature | replaces | source |
|---|---|---|
| `cdhit-embedded` | `cd-hit`, `cd-hit-est` | [henriksson-lab/cdhit-rs](https://github.com/henriksson-lab/cdhit-rs) (GPL-2.0-or-later) |
| `mafft-embedded` | `mafft` | [mahogny/rust-MAFFT](https://github.com/mahogny/rust-MAFFT), a fork of luksgrin/rust-MAFFT (MIT AND BSD-3-Clause) |

```sh
# source build only -- see below
cargo install --git https://github.com/henriksson-lab/panaroo-rs --features cli,cdhit-embedded,mafft-embedded
```

**This crate is not published to crates.io, and neither are the embedded backends** — by
decision, not oversight. Each is depended on wherever its latest code lives: cd-hit from
GitHub, MAFFT from a local checkout of the fork until its latest commit is pushed. So
`mafft-embedded` currently builds only on a machine with that checkout, which is also why
CI enables `cli,cdhit-embedded` and not `mafft-embedded`.

Both are verified byte-identical against the real tools on the parity suite — cd-hit on all
13 output files at `-t 1`, MAFFT on all 5,096 per-gene alignments of this dataset. **That
MAFFT figure does not yet generalise:** the parity work found a pre-existing gap-penalty
scaling bug in rust-MAFFT's DNA pairwise phase (L-INS-i, the mode `--auto` selects for
small clusters) that this dataset's conserved, indel-poor clusters happen not to trigger,
while realistic clusters with indels do (~50%). A fix is in progress in the fork; until it
lands and is re-verified, treat `mafft-embedded` as verified on this data, not in general.
Both integrations
drive the translated tool with the **same argv** the subprocess path would have built, so
flag semantics have exactly one definition. Note that `cdhit-embedded` pulls a GPL-2.0
dependency into the build.

Requires Rust **1.85** or newer (a floor set by `clap` and `indexmap`, not by this code).

The binary is named `panaroo`, matching upstream's CLI so it is a drop-in replacement.
**That means `cargo install` can shadow a Python Panaroo already on your `PATH`** — check
`which -a panaroo` if you have both.

`cd-hit` must be on `PATH`, and an aligner (`mafft` by default) if you use `--alignment`.
These are invoked as subprocesses; none of them is bundled.

## Verifying

```sh
tests/parity/env/create.sh && conda activate panaroo-parity
tests/parity/data/fetch.sh ci
tests/parity/reference/apply.sh
tests/parity/e2e.sh ci
```

See `tests/parity/README.md` for the rest of the suite.

## Performance

Same input, same output — all 13 output files byte-identical — on 4 *M. tuberculosis*
genomes (`--clean-mode strict`, 20 threads):

| | Python Panaroo | panaroo-rs | panaroo-rs `+cdhit-embedded` |
|---|---|---|---|
| wall clock | 64.0 s | 42.1 s (**1.5×**) | 44.0 s (1.5×) |
| CPU time (user + sys) | 612 s | 546 s (1.1×) | 530 s (1.2×) |
| peak memory (tree PSS) | 1691 MB | 462 MB (**3.7×**) | **256 MB** (**6.6×**) |

Three caveats, because a benchmark table without them is worse than none:

- **Single run, on a machine that was not idle** (load average 9.1 on 40 logical CPUs).
  Wall-clock time is the number most contaminated by that; treat 1.5× as a rough
  magnitude, not a measurement. Re-run with `tests/bench/run.sh ci -n 5 -T 20` on a quiet
  box for figures with a measured spread.
- **Almost all of this run is cd-hit, not our code.** `perf` puts **99.2%** of CPU inside
  the `cd-hit`/`cd-hit-est` subprocesses and **0.85%** in everything this crate wrote. So the
  wall-clock ratio is mostly a statement about process overhead and scheduling, and no
  amount of optimisation here can move it. The third column runs a Rust cd-hit in-process
  instead: same wall time, but **peak memory drops 462 → 256 MB** because it no longer forks
  a 20-thread C process holding its own word tables.
- **The CPU-time ratio (1.1×) is largely an artefact and should not be read as "the two
  implementations do about the same amount of work".** Both sides invoke the same external
  `cd-hit` binary 12 times per run with identical flags, and that binary is 60–75% of the
  CPU time. Worse, cd-hit's own CPU cost is *superlinear in its thread count* — measured at
  9.67 CPU-s with `-T 1` against 19.96 CPU-s with `-T 20`, a 2.07× penalty — so at `-t 20`
  most of both CPU columns is shared work that no port can change, inflated by cd-hit's
  threading overhead. The pipeline's own code is a minority of the total on this dataset.
  A meaningful CPU-time comparison needs `-t 1`, or needs cd-hit's time subtracted out.
- **The memory figure is tree PSS, not `time -v` max-RSS.** Max-RSS reports the largest
  single process, which for both implementations is the same `cd-hit` child (389 vs 388 MB
  — an artefact, not a result). The real difference is structural: Python parallelises by
  forking a `multiprocessing` pool (26 processes here), Rust with threads in one process
  (4 processes). Summed proportional set size across the whole process tree is what
  captures that.

Hardware: Xeon Gold 6138, 1 socket, 20 physical / 40 logical cores. `-t 20` is one thread
per physical core. Full harness, raw numbers and methodology: `tests/bench/`.

**What this configuration does *not* measure.** It runs without `-a/--alignment`, so mafft
is never invoked by either side and the whole alignment stage is skipped. That stage is
where the largest optimisation in the port lives — `output_sequence` used to re-parse the
whole combined FASTA once per gene cluster, ~99.8 GB per run, now parsed once — worth
**1.24× wall and 2.23× peak RSS** on `--alignment core`, and completely invisible here.
Conversely, that path is dominated by ~5,100 mafft invocations, which are not ours to
speed up. Neither configuration alone is representative.

## Reproducibility and how much parity to expect

Read this before relying on the output.

### Stock Panaroo is not reproducible run to run

Upstream iterates Python `set`s at points where the order reaches the output or the control
flow, and reads two directories in `os.listdir` order. Concretely, running stock Panaroo
twice on identical input gives different `pre_filt_graph.gml`, `final_graph.gml`,
`gene_presence_absence.csv` and `gene_presence_absence_roary.csv` unless `PYTHONHASHSEED` is
pinned — because `hash(str)` is randomised per process. The two `os.listdir` sites are worse
still: directory order is not stable across machines or filesystems, so the same input on
two computers can give a different column order in `core_gene_alignment.aln`.

We measured this: stock upstream under `PYTHONHASHSEED` 0, 1 and 2 produces four differing
files out of thirteen. **The pangenome itself does not move** — `gene_presence_absence.Rtab`
is byte-stable — so this is a rendering problem, not a science problem. But it does mean two
runs are not diffable, which makes any downstream regression check unreliable.

We believe this is an upstream bug rather than a deliberate choice, because Panaroo already
applies `sorted()` at two of these sites (`find_missing.py:52,77`) — the pattern is
established, just not applied consistently. It is written up as B6 in
`ORIGINAL_CODE_BUG.md`, in a form usable as an upstream issue report.

### A second source of nondeterminism: cd-hit itself, when threaded

The above concerns Panaroo's own Python. There is a second, independent source, and it sits
below both implementations: **upstream cd-hit produces nondeterministic output when run
multithreaded.** Panaroo passes `-t/--threads` straight through as cd-hit's `-T`, so any run
with `-t > 1` inherits it.

This matters more than it might appear, because cd-hit's clustering is the *first* stage:
its `.clstr` output becomes the nodes of the pangenome graph, so a difference there
propagates into every downstream file. It also means:

- **Panaroo's output is not reproducible at `-t > 1` regardless of implementation**, and no
  amount of determinism work inside Panaroo — including this port's Tier D patches — can fix
  it, because the nondeterminism is in a separate program.
- **A parity failure at `-t > 1` may be spurious.** Both sides invoke the same cd-hit binary
  independently, so the two runs can legitimately disagree without either implementation
  being wrong. Reproduce any failure at `-t 1` before treating it as a real defect.

What we have actually observed here: the parity suite runs at `-t 8` by default and has been
byte-identical across all 13 output files on the `ci` dataset consistently. That does not
disprove the nondeterminism — it suggests this dataset does not reliably trigger it — and it
should not be read as a guarantee. Treat `-t 1` as the reproducible configuration.

### What this port is compared against

Because a moving target cannot be matched byte-for-byte, parity here is measured against a
**patched reference**, built by `tests/parity/reference/apply.sh` from a pinned upstream
commit (`e96b497`) plus two patch sets, both enabled:

| tier | what | why |
|---|---|---|
| **D — canonicalisation** | 22 one-line edits across 7 files, replacing arbitrary `set` iteration and `os.listdir` order with a canonical order | without it there is no stable reference to compare against |
| **B — bug fixes** | the upstream defects in `ORIGINAL_CODE_BUG.md` | the stated goal was parity with canonicalised *and* bug-fixed behaviour |

`tests/parity/determinism.sh` verifies the Tier D reference gives byte-identical output
across `PYTHONHASHSEED` values.

### So: how much parity can you expect?

| against | expectation |
|---|---|
| the patched reference | **byte-identical**, all 13 files, verified on four real *M. tuberculosis* genomes and on both `--alignment` modes |
| stock Panaroo, same `PYTHONHASHSEED` | the same pangenome, differing by the Tier B bug fixes — **1–2 gene clusters out of ~5100** on our test data |
| stock Panaroo, different `PYTHONHASHSEED` | the same pangenome; four files differ in element ordering only |

If exact agreement with published stock-Panaroo numbers matters more than correctness,
disable the Tier B patches in `tests/parity/reference/enabled.txt` and re-verify — the port
implements the fixed behaviour, so that comparison will show the difference rather than
hide it.

### Known non-parity

- `--core_subset` is not implemented. Reproducing which genes survive needs a CPython
  Mersenne Twister clone for `random.shuffle`. It panics with that message rather than
  silently diverging.
- `alignment_resume_state.json` carries a wall-clock `started_at` and can never match; the
  comparator comes with that one field dropped.

## Citation

If you use this in published work, cite **Panaroo**, not this port:

> Tonkin-Hill, G., MacAlasdair, N., Ruis, C. et al. Producing polished prokaryotic
> pangenomes with the Panaroo pipeline. Genome Biol 21, 180 (2020).
> https://doi.org/10.1186/s13059-020-02090-4

Panaroo's post-processing scripts embed other algorithms that upstream asks you to cite
separately; those entry points are **not** translated here, so this port does not touch them.

## License

**MIT**, matching Panaroo. This is not a free choice: `panaroo-rs` is a derivative work of
Panaroo (MIT, Copyright (c) 2019 Gerry Tonkin-Hill), so it must carry a compatible licence
and preserve upstream's copyright notice. That notice is in `NOTICE.md`.

The `LICENSE` file carries this project's MIT grant **and** reproduces upstream Panaroo's
notice verbatim below it, as MIT requires of a derivative work.

Third-party licences, in full in `NOTICE.md`:

| | licence | how it is used |
|---|---|---|
| Panaroo | MIT | the program being translated |
| edlib v1.2.7 | MIT | **vendored and linked** (`vendor/edlib/`); its notice must ship with any binary |
| gffutils, biocode | MIT | behaviour reimplemented, no code copied |
| NetworkX, SciPy, NumPy, joblib, Biopython | BSD-3 / Biopython License | behaviour reimplemented |
| CPython | PSF-2.0 | `dict` and number-formatting behaviour reimplemented |
| intbitset | **LGPL-3.0-or-later** | behaviour reimplemented; not linked, not redistributed — see the note below |
| cd-hit, MAFFT, MUSCLE, PRANK, Clustal Omega, FAMSA | GPL-2.0 / GPL-3.0 / BSD-3 | invoked as **subprocesses**, never linked, never redistributed |

Two points worth a lawyer's eye before release:

- **intbitset is the one copyleft dependency.** `src/support/intbitset.rs` implements an
  ordinary bitset over non-negative integers from the operations Panaroo calls, without
  consulting intbitset's source. Nothing here links against or redistributes it.
- **The external aligners are GPL.** Running a GPL program as a subprocess is not linking,
  so it does not impose the GPL on this code — but do not bundle those binaries into a
  distribution without checking that separately.

