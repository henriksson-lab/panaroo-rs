#!/usr/bin/env python3
"""Download N RefSeq assemblies for a taxid, for the `scale` dataset.

Picks complete or chromosome-level assemblies deterministically (sorted by accession, so
re-running gives the same set) and fetches the _genomic.gff.gz and _genomic.fna.gz for
each. Panaroo's scripts/convert_refseq_to_prokka_gff.py then merges each pair into the
Prokka-style single file with an appended ##FASTA block that Panaroo expects.
"""
import argparse, gzip, os, shutil, sys, urllib.request

SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxid", type=int, required=True)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    cache = os.path.join(a.outdir, "assembly_summary.txt")
    if not os.path.exists(cache):
        print("fetching assembly summary (large, cached)...", file=sys.stderr)
        urllib.request.urlretrieve(SUMMARY, cache)

    rows = []
    with open(cache, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 20:
                continue
            # 5 = species_taxid, 11 = assembly_level, 19 = ftp_path
            if f[5] != str(a.taxid):
                continue
            if f[11] not in ("Complete Genome", "Chromosome"):
                continue
            if not f[19].startswith("http"):
                continue
            rows.append((f[0], f[19]))

    rows.sort()
    rows = rows[: a.n]
    if not rows:
        raise SystemExit(f"no assemblies found for taxid {a.taxid}")

    for acc, ftp in rows:
        stem = ftp.rstrip("/").rsplit("/", 1)[-1]
        for ext in ("genomic.gff", "genomic.fna"):
            dest = os.path.join(a.outdir, f"{acc}_{ext}")
            if os.path.exists(dest):
                continue
            url = f"{ftp}/{stem}_{ext}.gz"
            print(f"  {acc} {ext}", file=sys.stderr)
            tmp = dest + ".gz"
            urllib.request.urlretrieve(url, tmp)
            with gzip.open(tmp, "rb") as g, open(dest, "wb") as o:
                shutil.copyfileobj(g, o)
            os.remove(tmp)

    print(f"{len(rows)} assemblies in {a.outdir}")


if __name__ == "__main__":
    main()
