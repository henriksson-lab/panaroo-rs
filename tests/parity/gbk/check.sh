#!/usr/bin/env bash
# Differential test: biocode_convert::convert_gbk_gff3 vs the reference Python.
#
# Exercises support::genbank (the Bio.SeqIO GenBank reader) and support::biocode_gff3
# (biocode's GFF3 writer) end to end on a real RefSeq GenBank file.
#
# Requires: conda activate panaroo-parity
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
gbk="${1:-$repo/tests/parity/build/data/gbk/GCF_000759575.2.gbff}"
[[ -f "$gbk" ]] || { echo "no GenBank fixture at $gbk -- run tests/parity/data/fetch.sh gbk" >&2; exit 1; }
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

cargo build --release --example dump_gbk --manifest-path "$repo/Cargo.toml" >/dev/null
( cd "$repo/tests/parity/build/reference" && python3 -c "
import sys
from panaroo.biocode_convert import convert_gbk_gff3
convert_gbk_gff3('$gbk', '$tmp/py.gff', True)
" >/dev/null 2>&1 )
"$repo/target/release/examples/dump_gbk" "$gbk" "$tmp/rs.gff" >/dev/null 2>&1

if diff -q "$tmp/py.gff" "$tmp/rs.gff" >/dev/null; then
  echo "OK -- $(wc -l < "$tmp/py.gff") GFF3 lines match the reference conversion"
else
  diff "$tmp/py.gff" "$tmp/rs.gff" | head -10
  echo "FAIL"; exit 1
fi
