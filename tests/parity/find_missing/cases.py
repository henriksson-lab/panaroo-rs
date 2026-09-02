#!/usr/bin/env python3
"""Emit reference values for find_missing::{search_dna, translate_to_match, repl}.

Pairs with `cargo run --example dump_find_missing`.
"""
import os, sys, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build", "reference"))
from panaroo.find_missing import search_dna, translate_to_match

random.seed(7)
def rnd(n, alpha="ACGT"):
    return "".join(random.choice(alpha) for _ in range(n))

cases = []
gene = rnd(300)
cases.append(("exact", gene, "TTTT" + gene + "GGGG", 0.2, 0.95))
mut = list(gene); mut[10] = "A" if mut[10] != "A" else "C"; mut[150] = "T" if mut[150] != "T" else "G"
cases.append(("2 mismatches", "".join(mut), "AA" + gene + "CC", 0.2, 0.95))
cases.append(("truncated at contig end", gene, gene[100:], 0.2, 0.95))
cases.append(("absent", rnd(200), gene, 0.2, 0.95))
cases.append(("with N run", gene, gene[:50] + "N"*30 + gene[50:], 0.2, 0.95))
cases.append(("low threshold", rnd(120), gene, 0.05, 0.5))
cases.append(("revcomp", gene, "TT" + gene.translate(str.maketrans("ACGT","TGCA"))[::-1] + "GG", 0.2, 0.95))

# Randomised suite -- the hand-written cases above missed a real divergence on the smoke
# dataset, so generate a lot of shapes: long genes, indels, N-runs at both ends and in the
# middle, hits near contig boundaries, both strands.
def revc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

for i in range(400):
    L = random.choice([90, 150, 345, 512, 900, 1637])
    q = rnd(L)
    kind = i % 10
    flank_l, flank_r = rnd(random.randint(0, 60)), rnd(random.randint(0, 60))
    if kind == 0:
        db = flank_l + q + flank_r
    elif kind == 1:                     # substitutions
        m = list(q)
        for _ in range(max(1, L // 50)):
            j = random.randrange(L); m[j] = random.choice("ACGT")
        db = flank_l + "".join(m) + flank_r
    elif kind == 2:                     # deletion
        j = random.randrange(L - 10); db = flank_l + q[:j] + q[j + 9:] + flank_r
    elif kind == 3:                     # insertion
        j = random.randrange(L); db = flank_l + q[:j] + rnd(7) + q[j:] + flank_r
    elif kind == 4:                     # reverse complement
        db = flank_l + revc(q) + flank_r
    elif kind == 5:                     # truncated at the left contig end
        db = q[random.randrange(1, L // 2):] + flank_r
    elif kind == 6:                     # truncated at the right contig end
        db = flank_l + q[:random.randrange(L // 2, L)]
    elif kind == 7:                     # leading N run
        db = "N" * random.randint(20, 40) + q + flank_r
    elif kind == 8:                     # trailing N run
        db = flank_l + q + "N" * random.randint(20, 40)
    else:                               # internal N run
        j = random.randrange(L)
        db = flank_l + q[:j] + "N" * random.randint(20, 40) + q[j:] + flank_r
    cases.append((f"rand{i}-k{kind}-L{L}", q, db, 0.2, 0.95))

for name, q, db, pm, pid in cases:
    seq, loc = search_dna(db, q, pm, pid, True)
    # loc is [0, 0] on a miss and [start, end, strand] on a hit -- length is meaningful
    print(f"search_dna\t{name}\t{q}\t{db}\t{pm}\t{pid}\t{seq}\t{','.join(str(x) for x in loc)}")

for name, hit, prot in [
    ("fwd", "ATGAAATTTGGG", "MKFG"),
    ("empty", "", "MKFG"),
    ("revcomp", "CCCAAATTTCAT", "MKFG"),
    ("frameshift", "TATGAAATTTGGG", "MKFG"),
]:
    print(f"translate_to_match\t{name}\t{hit}\t{prot}\t{translate_to_match(hit, prot)}")
