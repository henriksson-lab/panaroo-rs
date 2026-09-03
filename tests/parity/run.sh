#!/usr/bin/env bash
# Parity harness: run the Python reference and panaroo-rs on the same input, compare.
#
#   ./run.sh DATASET [-- PANAROO_ARGS...]
#
#   DATASET   ci | smoke | scale, or a path to an input list
#
# Options:
#   -r DIR    reference tree to use (default: tests/parity/build/reference)
#   -p        run the Python reference only (use while the Rust side is incomplete)
#   -m MODE   --clean-mode (default: strict)
#   -t N      threads (default: 8)
#
# Requires the pinned env:  conda activate panaroo-parity
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
# shellcheck source=tests/lib/safe_rm.sh
source "$repo/tests/lib/safe_rm.sh"
build="$repo/tests/parity/build"

ref="$build/reference"
python_only=0
mode=strict
threads=8

dataset="${1:-smoke}"; shift || true
while getopts "r:pm:t:" opt; do
  case "$opt" in
    r) ref="$OPTARG" ;;
    p) python_only=1 ;;
    m) mode="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    *) echo "usage: $0 DATASET [-r REF] [-p] [-m MODE] [-t N] [-- ARGS...]" >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))
[[ "${1:-}" == "--" ]] && shift

if [[ -f "$dataset" ]]; then
  input="$dataset"; name="$(basename "$(dirname "$dataset")")"
else
  input="$build/data/$dataset/input.txt"; name="$dataset"
fi
[[ -f "$input" ]] || { echo "no input list at $input -- run tests/parity/data/fetch.sh $name" >&2; exit 1; }

# The reference is run from its own tree, so the input list must hold absolute paths.
if grep -qv '^/' "$input"; then
  echo "input list $input contains relative paths; regenerate with fetch.sh" >&2
  exit 1
fi

command -v cd-hit >/dev/null || { echo "cd-hit not on PATH -- conda activate panaroo-parity" >&2; exit 1; }

out_py="$build/out/$name/py"
out_rs="$build/out/$name/rs"
safe_rm_rf "$repo/tests" "$out_py"; mkdir -p "$out_py"

echo "== reference ($ref) on $name, $(wc -l < "$input") genomes =="
( cd "$ref" && PYTHONHASHSEED=0 python -m panaroo -i "$input" -o "$out_py" \
    --clean-mode "$mode" -t "$threads" "$@" ) || { echo "reference run failed" >&2; exit 1; }

if [[ "$python_only" == 1 ]]; then
  echo; echo "reference output: $out_py"; exit 0
fi

safe_rm_rf "$repo/tests" "$out_rs"; mkdir -p "$out_rs"
cargo build --release --features cli --manifest-path "$repo/Cargo.toml" >/dev/null
echo; echo "== panaroo-rs on $name =="
# panaroo-rs takes the file list directly, as upstream does when given many -i arguments
mapfile -t infiles < "$input"
"$repo/target/release/panaroo" -i "${infiles[@]}" -o "$out_rs" --clean-mode "$mode" -t "$threads" "$@"

echo; echo "== comparison =="
python3 "$here/canonicalise.py" "$out_py" "$out_rs"
