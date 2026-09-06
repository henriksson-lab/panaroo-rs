# panaroo-rs

A faithful *mostly*-Rust translation of [Panaroo](https://github.com/gtonkinhill/panaroo), the
prokaryotic pangenome pipeline.

The recommended build is mostly Rust: Panaroo itself plus in-process Rust translations of
cd-hit and MAFFT. Some less common aligners remain external tools when requested. The aim
of this translation is near-parity with higher speed.

Exact byte-parity is expected only in the single-threaded configuration (`-t 1`). Above
that, upstream cd-hit itself can produce nondeterministic clustering output.

**not yet tested enough**

* 2026-09-02: Initial translation


## Installing

The command-line tool is behind the **`cli` feature, which is off by default**, because
the crate's primary product is the library — the translated pipeline functions, which you
can drive directly. A library consumer should not have to build `clap`.

To get the recommended `panaroo` executable with in-process cd-hit and MAFFT, ask for
those features explicitly:

```sh
cargo install --git https://github.com/henriksson-lab/panaroo-rs --features cli,cdhit-embedded,mafft-embedded
```

or from a checkout:

```sh
cargo build --release --features cli,cdhit-embedded,mafft-embedded
```

To use it as a library instead, the default feature set is what you want:

```toml
[dependencies]
panaroo-rs = "0.1"        # library only, no clap
```

### In-process cd-hit and MAFFT

For the current recommended build, cd-hit and MAFFT run **in-process** via Rust
translations:

| feature | replaces | source |
|---|---|---|
| `cdhit-embedded` | `cd-hit`, `cd-hit-est` | [henriksson-lab/cdhit-rs](https://github.com/henriksson-lab/cdhit-rs) (GPL-2.0-or-later) |
| `mafft-embedded` | `mafft` | [mahogny/rust-MAFFT](https://github.com/mahogny/rust-MAFFT), a fork of luksgrin/rust-MAFFT (MIT AND BSD-3-Clause) |

**This crate is not published to crates.io, and neither are the embedded backends** — by
decision, not oversight. The embedded backends are pinned to Git revisions that were
checked with this tree, so source installs can build the in-process configuration without
local path dependencies.

Both integrations drive the translated tool with the **same argv** the original external
tool invocation would have used, so flag semantics have exactly one definition. On the
current parity corpus, the fully embedded single-thread path is byte-identical except for
`alignment_resume_state.json`'s wall-clock timestamp field. Note that `cdhit-embedded`
pulls a GPL-2.0 dependency into the build.

Requires Rust **1.85** or newer (a floor set by `clap` and `indexmap`, not by this code).

The binary is named `panaroo`, matching upstream's CLI so it is a drop-in replacement.
**That means `cargo install` can shadow a Python Panaroo already on your `PATH`** — check
`which -a panaroo` if you have both.

## Verifying

```sh
tests/parity/env/create.sh && conda activate panaroo-parity
tests/parity/data/fetch.sh ci
tests/parity/reference/apply.sh
tests/parity/e2e.sh ci
```

See `tests/parity/README.md` for the rest of the suite.

## Performance

Current embedded-backend measurements use the 4-genome *M. tuberculosis* parity input,
`--clean-mode strict`, `-t 1`, and byte-identical outputs unless noted.

| workload | original tool/path | Rust path | result |
|---|---:|---:|---:|
| direct protein `cd-hit` | 7.02 s | 5.01 s | 1.40x faster |
| direct nucleotide `cd-hit-est` | 23.74 s | 17.01 s | 1.40x faster |
| 451 per-gene MAFFT alignments | ~311 s | 60.06-60.18 s | ~5.2x faster |
| full Panaroo, `-a core`, `+cdhit-embedded,+mafft-embedded` | external reference | 78.93 s | timestamp-only diff |

The direct cd-hit comparisons use the exact intermediate FASTA files and flags generated
by Panaroo. Both the clustered FASTA and `.clstr` outputs were byte-identical, so the
current Rust cd-hit path is not slower than bundled CD-HIT on this workload despite doing
the hot diagonal tests with SIMD.

The MAFFT figure is a corpus harness over Panaroo's unaligned gene clusters. It is a good
measurement of the embedded alignment path used by this dataset, not a claim that
`rust-MAFFT` is a complete replacement for every MAFFT mode and input shape.

The older 20-thread no-alignment benchmark is still useful as a throughput smoke test,
but not as the parity baseline:

| 20-thread no-alignment run | Python Panaroo | panaroo-rs | panaroo-rs `+cdhit-embedded` |
|---|---:|---:|---:|
| wall clock | 64.0 s | 42.1 s (1.52x faster) | 44.0 s (1.45x faster) |
| CPU time (user + sys) | 612 s | 546 s (1.12x faster) | 530 s (1.15x faster) |
| peak memory (tree PSS) | 1691 MB | 462 MB (3.7x lower) | 256 MB (6.6x lower) |

That run used `--clean-mode strict`, `-t 20`, and no `-a/--alignment`, so MAFFT was never
invoked. All 13 no-alignment outputs were byte-identical in that run, but exact parity
should still be judged at `-t 1`, because threaded cd-hit can be nondeterministic.

Hardware for these runs: Xeon Gold 6138, 1 socket, 20 physical / 40 logical cores. Full
harness notes, raw numbers and rejected optimisation experiments live under
`tests/bench/`.

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
| cd-hit, MAFFT, MUSCLE, PRANK, Clustal Omega, FAMSA | GPL-2.0 / GPL-3.0 / BSD-3 | cd-hit and MAFFT are linked via Rust translations in the recommended build; other aligners remain external when requested |

Two points worth a lawyer's eye before release:

- **intbitset is the one copyleft dependency.** `src/support/intbitset.rs` implements an
  ordinary bitset over non-negative integers from the operations Panaroo calls, without
  consulting intbitset's source. Nothing here links against or redistributes it.
- **The aligner licensing needs review before release.** The recommended build links the
  GPL-2.0-or-later `cdhit-rs` translation. Less common aligners remain external tools when
  requested, and bundling any original third-party binaries needs separate review.
