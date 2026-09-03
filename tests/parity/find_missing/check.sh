#!/usr/bin/env bash
# Differential test: find_missing::{search_dna, translate_to_match} vs the reference Python.
#
# These two are intricate enough to be worth checking directly rather than only through the
# end-to-end run: search_dna pads with 'E', aligns both strands, picks the location closest
# to the centre using *signed* differences, and post-processes N-runs with two regexes.
#
# Requires: conda activate panaroo-parity
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

cargo build --release --example dump_find_missing --manifest-path "$repo/Cargo.toml" >/dev/null
python3 "$here/cases.py" > "$tmp/py.tsv"
"$repo/target/release/examples/dump_find_missing" "$tmp/py.tsv" > "$tmp/rs.tsv"

# compare the result columns (7 = sequence, 8 = loc); inputs are echoed back verbatim
if diff <(cut -f1,2,7,8 "$tmp/py.tsv") <(cut -f1,2,7,8 "$tmp/rs.tsv"); then
  echo "OK -- $(wc -l < "$tmp/py.tsv") cases match the reference"
else
  echo "FAIL"; exit 1
fi
