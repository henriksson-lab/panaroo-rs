#!/usr/bin/env bash
# Build the parity reference: a copy of the pinned upstream Panaroo with the patches
# listed in enabled.txt applied.
#
#   ./apply.sh [-o OUTDIR] [-e ENABLED_FILE] [-w EXTRA_PATCH]...
#
# -w applies a patch that is not in enabled.txt, for one-off experiments
#    (e.g. measuring what a Tier B fix actually changes) without editing the file.
#
# The output tree is disposable; it is rebuilt from scratch on every run.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
# shellcheck source=tests/lib/safe_rm.sh
source "$repo/tests/lib/safe_rm.sh"
upstream="$repo/panaroo"
outdir="$repo/tests/parity/build/reference"
enabled="$here/enabled.txt"
extra=()

while getopts "o:e:w:" opt; do
  case "$opt" in
    o) outdir="$OPTARG" ;;
    e) enabled="$OPTARG" ;;
    w) extra+=("$OPTARG") ;;
    *) echo "usage: $0 [-o OUTDIR] [-e ENABLED_FILE] [-w EXTRA_PATCH]..." >&2; exit 2 ;;
  esac
done

[[ -d "$upstream/panaroo" ]] || { echo "no upstream checkout at $upstream" >&2; exit 1; }

# Refuse to build from a dirty upstream: the reference must be reproducible from a commit.
if ! git -C "$upstream" diff --quiet || ! git -C "$upstream" diff --cached --quiet; then
  echo "upstream checkout $upstream has uncommitted changes; refusing to build" >&2
  exit 1
fi
rev="$(git -C "$upstream" rev-parse HEAD)"

safe_rm_rf "$repo/tests" "$outdir"
mkdir -p "$(dirname "$outdir")"
git -C "$upstream" archive --format=tar HEAD | (mkdir -p "$outdir" && tar -x -C "$outdir")

applied=()
while IFS= read -r name; do
  name="${name%%#*}"; name="$(echo "$name" | xargs)"
  [[ -z "$name" ]] && continue
  applied+=("$here/patches/$name")
done < "$enabled"
applied+=("${extra[@]+"${extra[@]}"}")

# Use GNU patch, NOT `git apply`. `git apply` resolves paths against the enclosing git
# repository root, and $outdir lives inside this repo -- so it silently reports
# "Skipped patch" for every file and exits 0, producing an unpatched tree that looks
# successful. patch(1) has no repo awareness and honours -d.
for p in ${applied[@]+"${applied[@]}"}; do
  [[ -f "$p" ]] || { echo "missing patch: $p" >&2; exit 1; }
  if ! patch -p1 -d "$outdir" --forward --batch --silent < "$p"; then
    echo "FAILED to apply $(basename "$p")" >&2
    exit 1
  fi
  echo "applied $(basename "$p")"
done

# Verify rather than trust. A patch that applies to zero hunks is a silent no-op -- that
# is exactly how `git apply` failed here before. Re-extract a pristine tree and require
# that the set of files which actually differ equals the set the patches claim to touch.
pristine="$(mktemp -d)"
trap 'rm -rf "$pristine"' EXIT
git -C "$upstream" archive --format=tar HEAD | tar -x -C "$pristine"

expected="$(for p in ${applied[@]+"${applied[@]}"}; do
             grep '^+++ b/' "$p" | sed 's|^+++ b/||'
           done | sort -u)"
# `diff -rq` exits 1 when files differ, which is the expected case here -- guard it or
# `set -o pipefail` turns a successful verification into a script failure.
actual="$( { diff -rq "$pristine" "$outdir" 2>/dev/null || true; } \
          | sed -n 's|^Files '"$pristine"'/\(.*\) and .* differ$|\1|p' | sort -u)"

if [[ "$expected" != "$actual" ]]; then
  echo "patch verification FAILED: files changed on disk do not match files named in the patches" >&2
  echo "  expected:" >&2; echo "$expected" | sed 's/^/    /' >&2
  echo "  actual:"   >&2; echo "$actual"   | sed 's/^/    /' >&2
  exit 1
fi
n_expected=$(echo "$expected" | grep -c . || true)
echo "verified: $n_expected file(s) changed, exactly as the patches specify"

{
  echo "upstream_rev=$rev"
  for p in ${applied[@]+"${applied[@]}"}; do echo "patch=$(basename "$p")"; done
} > "$outdir/PARITY_REFERENCE"

echo
echo "reference built at $outdir"
echo "  upstream $rev"
echo "  ${#applied[@]} patch(es) applied"
