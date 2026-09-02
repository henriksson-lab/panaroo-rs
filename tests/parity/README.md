# Parity harness

Everything needed to check `panaroo-rs` against the Python reference. See
`PORTING_PLAN.md` §6 and §10 for what "parity" means here and why it needs this much
machinery.

## One-time setup

```sh
tests/parity/env/create.sh          # pinned conda env: cd-hit, mafft, biopython, ...
conda activate panaroo-parity
tests/parity/data/fetch.sh ci       # real M. tuberculosis genomes from upstream's CI
tests/parity/data/fetch.sh smoke 3  # fast derived subset for iteration
tests/parity/reference/apply.sh     # build the patched Python reference
```

## Running

```sh
tests/parity/e2e.sh ci                    # full pipeline, both sides, compare
tests/parity/e2e.sh tiny -- -a core       # with the core alignment stage
tests/parity/determinism.sh smoke         # Tier D acceptance: stable across hash seeds

tests/parity/prokka/check.sh ci           # input stage in isolation
tests/parity/gff/check.sh                 # support::gff vs gffutils
tests/parity/graph/check.sh               # support::graph vs networkx
tests/parity/find_missing/check.sh        # search_dna / translate_to_match
tests/parity/gbk/check.sh                 # GenBank -> GFF3 vs biocode
```

`tests/parity/run.sh` is the older combined runner; `e2e.sh` supersedes it.

## Layout

| path | what |
|---|---|
| `env/` | pinned conda environment; tool versions are part of the parity contract |
| `data/fetch.sh` | downloads/derives the datasets; nothing here is committed |
| `reference/apply.sh` | builds the patched Python reference from the pinned upstream commit |
| `reference/patches/` | one `.patch` per change, each carrying its own written rationale |
| `reference/enabled.txt` | which patches apply; Tier D on, Tier B off |
| `canonicalise.py` | Tier-1 normalisation + directory comparison |
| `e2e.sh` | full-pipeline comparison |
| `determinism.sh` | Tier D acceptance test |
| `<module>/check.sh` | per-module differential tests against the library each replaces |
| `run.sh` | older combined runner, superseded by `e2e.sh` |
| `build/` | generated, gitignored |

## Reading the comparison

`canonicalise.py A B` reports per file:

- `identical (tier 0)` — byte-for-byte. The goal for every file.
- `identical after canonicalisation (tier 1)` — same pangenome, different rendering of a
  Python `set`'s iteration order. Acceptable *only* for `seqIDs`/`geneIDs` in the GMLs and
  `;`-token order in the presence/absence CSVs.
- `DIFFERS` — a real difference. Investigate.

A file moving from tier 0 to tier 1 is a regression worth explaining, not a pass.

## Gotchas

- **`input.txt` must contain absolute paths.** The reference runs from its own tree, so
  relative paths resolve against the wrong directory and Panaroo fails with the unhelpful
  `RuntimeError: Error reading prokka input!`. `fetch.sh` writes absolute paths; `run.sh`
  checks.
- **`biocode` is required** despite `panaroo/biocode_convert.py` looking vendored — it does
  `from biocode import ...` and `prokka.py` imports it at module scope, so `python -m
  panaroo` fails at import time without it. It is pip-only.
- `PYTHONHASHSEED=0` is set for every reference run. With the Tier D patch the reference is
  reproducible without it too, and `determinism.sh` checks exactly that.
- `alignment_resume_state.json` carries a wall-clock `started_at`, so it can only ever match
  at tier 1. `canonicalise.py` drops that one field and compares the rest.
