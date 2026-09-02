#!/usr/bin/env bash
# Differential test: src/support/graph.rs vs real networkx.
#
# Replays tests/parity/graph/ops.txt through both and requires identical output. The ops
# file targets the order-sensitive behaviour that byte parity depends on: node insertion
# order, adjacency insertion order, the seen-dedup in edges()/edges(nbunch), BFS level
# order, and component order.
#
# Requires: conda activate panaroo-parity
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

python3 -c "import networkx" 2>/dev/null || { echo "networkx not importable -- conda activate panaroo-parity" >&2; exit 1; }
cargo build --release --example dump_graph --manifest-path "$repo/Cargo.toml" >/dev/null 2>&1

ops="${1:-$here/ops.txt}"
python3 "$here/dump_networkx.py" "$ops" > "$tmp/py.txt"
"$repo/target/release/examples/dump_graph" "$ops" > "$tmp/rs.txt"

if diff -u "$tmp/py.txt" "$tmp/rs.txt"; then
  echo "OK -- $(wc -l < "$tmp/py.txt") assertions match networkx $(python3 -c 'import networkx;print(networkx.__version__)')"
else
  echo "FAIL"; exit 1
fi
