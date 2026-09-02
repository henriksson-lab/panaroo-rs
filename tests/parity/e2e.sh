#!/usr/bin/env bash
# End-to-end parity: the full pipeline must be byte-identical to the reference.
#
#   tests/parity/e2e.sh [DATASET] [-- EXTRA_ARGS...]      (default dataset: ci)
#
# Requires: conda activate panaroo-parity
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
build="$repo/tests/parity/build"
dataset="${1:-ci}"; shift || true
[[ "${1:-}" == "--" ]] && shift
threads="${THREADS:-8}"
mode="${MODE:-strict}"

input="$build/data/$dataset/input.txt"
[[ -f "$input" ]] || { echo "no input at $input -- run tests/parity/data/fetch.sh $dataset" >&2; exit 1; }
command -v cd-hit >/dev/null || { echo "cd-hit not on PATH -- conda activate panaroo-parity" >&2; exit 1; }

out="$build/out/e2e/$dataset"
rm -rf "$out"; mkdir -p "$out/py" "$out/rs"
cargo build --release --features cli --manifest-path "$repo/Cargo.toml" >/dev/null 2>&1

echo "== reference =="
( cd "$build/reference" && PYTHONHASHSEED=0 python -m panaroo \
    -i "$input" -o "$out/py" --clean-mode "$mode" -t "$threads" "$@" >/dev/null )

echo "== panaroo-rs =="
mapfile -t infiles < "$input"
"$repo/target/release/panaroo" -i "${infiles[@]}" -o "$out/rs" --clean-mode "$mode" -t "$threads" "$@" >/dev/null

echo
python3 "$here/canonicalise.py" "$out/py" "$out/rs"
