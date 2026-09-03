#!/usr/bin/env bash
# Fetch parity-testing datasets into tests/parity/build/data/.
#
#   ./fetch.sh ci       upstream Panaroo CI data: 4 real M. tuberculosis draft assemblies
#                       (Prokka GFF3 + FASTA, ~5100 CDS / 35 contigs each), plus the
#                       paralog-rich and alignment fixtures upstream's own tests use.
#                       This is the authoritative set: it is what upstream CI validates.
#
#   ./fetch.sh smoke    derived from `ci` -- the first N contigs of each genome, for fast
#                       iteration during the port. Deterministic, regenerable, not committed.
#
#   ./fetch.sh gbk      one RefSeq GenBank flat file (Campylobacter jejuni), the fixture
#                       for the biocode_convert GenBank-to-GFF3 path
#
#   ./fetch.sh scale N SPECIES_TAXID
#                       N RefSeq assemblies converted to Prokka-style GFF3, for exercising
#                       collapse_families at a realistic genome count.
#
# Data is never committed: it is large and it is reproducible from here.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
# shellcheck source=tests/lib/safe_rm.sh
source "$repo/tests/lib/safe_rm.sh"
out="$repo/tests/parity/build/data"

CI_URL="https://github.com/gtonkinhill/panaroo_test_data/releases/download/v0.0.2/travis_test_data.zip"
CI_SHA="25bde7222ff092700cf792c126443f59c6cc9ee799184ac2c24d3f4c8a7612b0"

fetch_ci() {
  mkdir -p "$out"
  local zip="$out/travis_test_data.zip"
  if [[ ! -f "$zip" ]]; then
    echo "downloading upstream CI data..."
    curl -fL --retry 3 -o "$zip" "$CI_URL"
  fi
  echo "$CI_SHA  $zip" | sha256sum -c - || { echo "checksum mismatch; refusing" >&2; exit 1; }
  safe_rm_rf "$repo/tests" "$out/ci" "$out/__MACOSX"
  unzip -q -o "$zip" -d "$out"
  # The release zip was built on macOS and carries a __MACOSX sidecar tree of AppleDouble
  # "._name" stubs. They have .gff names but are not GFF3, so anything globbing *.gff picks
  # them up and chokes.
  safe_rm_rf "$repo/tests" "$out/__MACOSX"
  find "$out/travis_test_data" -name '._*' -delete
  mv "$out/travis_test_data" "$out/ci"
  # aa1..aa4 are the four real M. tuberculosis draft assemblies. aln.gff (9 CDS) and
  # paralog.gff are single-purpose fixtures for upstream's own unit tests, not genomes to
  # build a pangenome from -- keep them out of the default input list.
  ls "$out"/ci/aa[1-4].gff | xargs -n1 readlink -f > "$out/ci/input.txt"
  ls "$out"/ci/*.gff | xargs -n1 readlink -f > "$out/ci/input_all.txt"
  echo "ci: $(wc -l < "$out/ci/input.txt") genomes -> $out/ci/input.txt"
  echo "    (all $(wc -l < "$out/ci/input_all.txt") files incl. fixtures: input_all.txt)"
}

# Keep only the first N contigs of each genome, and only the CDS on them. Produces a
# valid Prokka-style GFF3 an order of magnitude faster to run than the full set.
fetch_smoke() {
  local n="${1:-4}"
  [[ -d "$out/ci" ]] || fetch_ci
  safe_rm_rf "$repo/tests" "$out/smoke"; mkdir -p "$out/smoke"
  for f in "$out"/ci/aa[1-4].gff; do
    python3 "$here/subset_gff.py" "$f" "$out/smoke/$(basename "$f")" --contigs "$n"
  done
  ls "$out"/smoke/*.gff | xargs -n1 readlink -f > "$out/smoke/input.txt"
  echo "smoke: $(ls "$out"/smoke/*.gff | wc -l) gff files, first $n contigs each -> $out/smoke"
}

fetch_scale() {
  local n="${1:-20}" taxid="${2:-573}"   # 573 = Klebsiella pneumoniae
  safe_rm_rf "$repo/tests" "$out/scale"; mkdir -p "$out/scale/raw"
  python3 "$here/fetch_refseq.py" --taxid "$taxid" --n "$n" --outdir "$out/scale/raw"
  local conv="$repo/panaroo/scripts/convert_refseq_to_prokka_gff.py"
  for g in "$out"/scale/raw/*_genomic.gff; do
    local base="${g%_genomic.gff}"
    python3 "$conv" -g "$g" -f "${base}_genomic.fna" -o "$out/scale/$(basename "$base").gff"
  done
  ls "$out"/scale/*.gff | xargs -n1 readlink -f > "$out/scale/input.txt"
  echo "scale: $(ls "$out"/scale/*.gff | wc -l) gff files -> $out/scale"
}

fetch_gbk() {
  mkdir -p "$out/gbk"
  local dest="$out/gbk/GCF_000759575.2.gbff"
  [[ -f "$dest" ]] && { echo "gbk: already present"; return; }
  curl -fsSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/759/575/GCF_000759575.2_ASM75957v2/GCF_000759575.2_ASM75957v2_genomic.gbff.gz" \
    | gunzip > "$dest"
  echo "gbk: $dest"
}

case "${1:-ci}" in
  ci)    fetch_ci ;;
  gbk)   fetch_gbk ;;
  smoke) shift; fetch_smoke "$@" ;;
  scale) shift; fetch_scale "$@" ;;
  *)     echo "usage: $0 {ci|smoke [N_CONTIGS]|gbk|scale [N_GENOMES] [TAXID]}" >&2; exit 2 ;;
esac
