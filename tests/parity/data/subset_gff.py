#!/usr/bin/env python3
"""Cut a Prokka-style GFF3 down to its first N contigs, keeping the file valid.

Used to derive the fast `smoke` dataset from the real `ci` genomes. Deterministic: the
contig order is the order they appear in the ##FASTA section, so re-running gives the same
subset. Both the annotation block and the FASTA block are filtered consistently, and
##sequence-region lines for dropped contigs are removed.
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--contigs", type=int, default=4)
    a = ap.parse_args()

    text = open(a.infile).read()
    ann, _, fasta = text.partition("##FASTA")
    if not fasta:
        raise SystemExit(f"{a.infile}: no ##FASTA section")

    # Split the FASTA block into records, preserving the original text verbatim.
    records, cur = [], []
    for line in fasta.splitlines(keepends=True):
        if line.startswith(">"):
            if cur:
                records.append(cur)
            cur = [line]
        elif cur:
            cur.append(line)
    if cur:
        records.append(cur)

    keep_recs = records[: a.contigs]
    keep_ids = {r[0][1:].split()[0] for r in keep_recs}

    out_ann = []
    for line in ann.splitlines(keepends=True):
        if line.startswith("##sequence-region"):
            parts = line.split()
            if len(parts) > 1 and parts[1] not in keep_ids:
                continue
            out_ann.append(line)
        elif line.startswith("#"):
            out_ann.append(line)
        elif line.strip():
            if line.split("\t", 1)[0] in keep_ids:
                out_ann.append(line)

    with open(a.outfile, "w") as fh:
        fh.writelines(out_ann)
        fh.write("##FASTA\n")
        for r in keep_recs:
            fh.writelines(r)

    n_cds = sum(1 for l in out_ann if not l.startswith("#") and l.split("\t")[2:3] == ["CDS"])
    print(f"{a.outfile}: {len(keep_recs)} contigs, {n_cds} CDS")


if __name__ == "__main__":
    main()
