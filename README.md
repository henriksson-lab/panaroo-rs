# panaroo-rs

A faithful *mostly*-Rust translation of [Panaroo](https://github.com/gtonkinhill/panaroo), the
prokaryotic pangenome pipeline.

Note that this crate is not fully Rust yet. Translation of upstream dependencies is needed first. The aim of this translation is to improve speed over the original Panaroo code

Also note that the Panaroo code output is not reproduced due to translation challenges and a suspected upstream bug

**not yet tested enough**

* 2026-09-02: Initial translation



## Status

All 92 functions on the `panaroo` entry point are translated, and the pipeline produces
**byte-identical output** to its reference on real data — 13/13 files across four
*M. tuberculosis* genomes, including both `--alignment` modes.

"Its reference" is doing real work in that sentence: it is a *patched* Panaroo, not stock.
[Reproducibility and how much parity to expect](#reproducibility-and-how-much-parity-to-expect)
explains why, and what that means for you.

**No optimisation has been done yet.** The port deliberately reproduces upstream's hot
spots so that parity could be established first; they are catalogued in `PORTING_PLAN.md`
§9 and are the obvious next work. Current runtime is roughly 20% below the Python.

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

## Documentation

| file | what |
|---|---|
| `PORTING_PLAN.md` | the translation plan, parity strategy, and deferred optimisations |
| `ORIGINAL_CODE_BUG.md` | bugs found in upstream Panaroo, with measured effects |
| `NOTICE.md` | third-party attribution — what is reimplemented, vendored, or wrapped |
| `port_order.csv` | the 92-function bottom-up checklist |
