#!/usr/bin/env bash
# Acceptance test for the Tier D patch: the reference must produce byte-identical output
# regardless of PYTHONHASHSEED.
#
#   tests/parity/determinism.sh [DATASET]     (default: smoke)
#
# This is what caught the site the first Tier D pass missed -- custom_stringizer's `set`
# branch in isvalid.py, which only shows up in pre_filt_graph.gml. Re-run it after any
# change to the reference patches.
#
# Requires: conda activate panaroo-parity
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
# shellcheck source=tests/lib/safe_rm.sh
source "$repo/tests/lib/safe_rm.sh"
build="$repo/tests/parity/build"
dataset="${1:-smoke}"
input="$build/data/$dataset/input.txt"
threads="${THREADS:-8}"

[[ -f "$input" ]] || { echo "no input at $input -- run tests/parity/data/fetch.sh $dataset" >&2; exit 1; }
command -v cd-hit >/dev/null || { echo "cd-hit not on PATH -- conda activate panaroo-parity" >&2; exit 1; }

out="$build/out/determinism/$dataset"
safe_rm_rf "$repo/tests" "$out"; mkdir -p "$out"
for s in 0 1 2; do
  echo "== seed $s =="
  ( cd "$build/reference" && PYTHONHASHSEED=$s python -m panaroo \
      -i "$input" -o "$out/seed$s" --clean-mode strict -t "$threads" >/dev/null 2>&1 )
done

echo
bad=0
for f in $(ls "$out/seed0"); do
  if cmp -s "$out/seed0/$f" "$out/seed1/$f" && cmp -s "$out/seed0/$f" "$out/seed2/$f"; then
    printf "  %-36s identical\n" "$f"
  else
    printf "  %-36s DIFFERS\n" "$f"; bad=1
  fi
done
echo
if [[ $bad -eq 0 ]]; then echo "PASS -- reference is reproducible across hash seeds"; else
  echo "FAIL -- a set-iteration site is still unpatched; diff the offending file to find which key moves"; exit 1
fi
