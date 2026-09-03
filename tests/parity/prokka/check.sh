#!/usr/bin/env bash
# Phase 3 checkpoint: the input stage must be byte-identical to the reference.
#
#   tests/parity/prokka/check.sh [DATASET]     (default: ci)
#
# Covers prokka::{get_gene_sequences, translate_sequences, output_files,
# process_prokka_input} plus support::{gff, seqio, seq, codon_table, parallel}.
# Nothing in this stage touches a Python `set`, so tier 0 is the right bar.
#
# Requires: conda activate panaroo-parity
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
# shellcheck source=tests/lib/safe_rm.sh
source "$repo/tests/lib/safe_rm.sh"
dataset="${1:-ci}"
threads="${THREADS:-8}"
input="$repo/tests/parity/build/data/$dataset/input.txt"
[[ -f "$input" ]] || { echo "no input at $input -- run tests/parity/data/fetch.sh $dataset" >&2; exit 1; }

out="$repo/tests/parity/build/out/prokka/$dataset"
safe_rm_rf "$repo/tests" "$out"; mkdir -p "$out/py" "$out/rs"

cargo build --release --example run_prokka_stage --manifest-path "$repo/Cargo.toml" >/dev/null
PYTHONHASHSEED=0 python3 "$here/run_reference.py" "$input" "$out/py" "$threads" >/dev/null 2>&1
"$repo/target/release/examples/run_prokka_stage" "$input" "$out/rs" "$threads" >/dev/null 2>&1

bad=0
for f in gene_data.csv combined_DNA_CDS.fasta combined_protein_CDS.fasta; do
  if cmp -s "$out/py/$f" "$out/rs/$f"; then
    printf "  %-28s IDENTICAL (%s bytes)\n" "$f" "$(stat -c%s "$out/py/$f")"
  else
    printf "  %-28s DIFFERS\n" "$f"
    diff <(head -c 4000 "$out/py/$f") <(head -c 4000 "$out/rs/$f") | head -8 | sed 's/^/      /'
    bad=1
  fi
done
echo
[[ $bad -eq 0 ]] && echo "PASS -- input stage is byte-identical on '$dataset'" || { echo "FAIL"; exit 1; }
