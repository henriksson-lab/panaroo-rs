#!/usr/bin/env bash
# Differential test: src/support/gff.rs vs the real gffutils, on real data.
#
#   tests/parity/gff/check.sh [FILE.gff ...]
#
# With no arguments, runs over every dataset fetched into tests/parity/build/data/.
# Requires: conda activate panaroo-parity
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

command -v python3 >/dev/null && python3 -c "import gffutils" 2>/dev/null \
  || { echo "gffutils not importable -- conda activate panaroo-parity" >&2; exit 1; }

cargo build --release --example dump_gff --manifest-path "$repo/Cargo.toml" >/dev/null
dumper="$repo/target/release/examples/dump_gff"

files=("$@")
if [[ ${#files[@]} -eq 0 ]]; then
  mapfile -t files < <(find "$repo/tests/parity/build/data" -name '*.gff' \
      -not -path '*/raw/*' -not -path '*__MACOSX*' -not -name '._*' | sort)
fi
[[ ${#files[@]} -gt 0 ]] || { echo "no .gff files found -- run tests/parity/data/fetch.sh" >&2; exit 1; }

bad=0 total=0
for f in "${files[@]}"; do
  n="$(basename "$f")"
  python3 "$here/dump_gffutils.py" "$f" > "$tmp/py.tsv" 2>/dev/null
  "$dumper" "$f" > "$tmp/rs.tsv"
  count=$(wc -l < "$tmp/py.tsv")
  total=$((total + count))
  if cmp -s "$tmp/py.tsv" "$tmp/rs.tsv"; then
    printf "  %-34s IDENTICAL (%s features)\n" "$n" "$count"
  else
    printf "  %-34s DIFFERS\n" "$n"
    diff "$tmp/py.tsv" "$tmp/rs.tsv" | head -6 | sed 's/^/      /'
    bad=$((bad + 1))
  fi
done

echo
if [[ $bad -eq 0 ]]; then
  echo "OK -- $total features matched gffutils exactly across ${#files[@]} file(s)"
else
  echo "FAIL -- $bad of ${#files[@]} file(s) differ"
  exit 1
fi
