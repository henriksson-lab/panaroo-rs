#!/usr/bin/env python3
"""Dump a GFF3 file's parse via gffutils, in the same canonical form as
`cargo run --example dump_gff`. Diffing the two validates src/support/gff.rs.

Input handling mirrors prokka.get_gene_sequences: strip commas, split off ##FASTA, drop
##sequence-region lines.
"""
import sys, gffutils

path = sys.argv[1]
text = open(path).read().replace(",", "")
ann = text.split("##FASTA")[0]
ann = "\n".join(l for l in ann.splitlines() if "##sequence-region" not in l)

db = gffutils.create_db(ann, dbfn=":memory:", force=True, keep_order=True,
                        from_string=True, merge_strategy="create_unique")

out = sys.stdout
for e in db.all_features(featuretype=()):
    attrs = "\x02".join(f"{k}={chr(1).join(v)}" for k, v in e.attributes.items())
    out.write("\t".join([
        e.id, e.seqid, e.source, e.featuretype,
        str(e.start), str(e.stop), e.score, e.strand, e.frame, attrs,
    ]) + "\n")
