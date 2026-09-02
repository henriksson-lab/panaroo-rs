# Bugs in upstream Panaroo

Findings from reading `gtonkinhill/panaroo` @ `e96b497` ("version bump", v1.8.0) closely
enough to translate it function-for-function into Rust.

These are recorded here because a faithful port has to decide, for each one, whether to
reproduce it or fix it — and because they are worth reporting upstream. What the port does
with them is in `PORTING_PLAN.md` §6.1 (Tier B patches, off by default) and §9.

**Confidence key**

| | meaning |
|---|---|
| **measured** | reproduced, patched, and the output difference quantified |
| **read** | derived from reading the source; consequence not yet measured |

| # | site | severity | confidence |
|---|---|---|---|
| [B1](#b1) | `find_missing.py:61` — `mem` read from a leaked loop variable | **affects results** | measured |
| [B2](#b2) | `find_missing.py:329` — validity check reads loop leftovers | **affects results** under `--refind-mode strict` | read |
| [B3](#b3) | `cdhit.py:535` — `pqid` typo | none (provably inert) | read |
| [B4](#b4) | `prokka.py:171` — contig scan never breaks | performance; latent correctness | read |
| [B5](#b5) | `__main__.py:190` — help text contradicts the default | documentation | read |
| [B6](#b6) | pipeline output is not reproducible run-to-run | reproducibility | measured |

Environment for all measurements: `tests/parity/env/environment.yml` (Python 3.11,
cd-hit 4.8.1, mafft 7.526, biopython 1.84, numpy 1.26.4, scipy 1.14.1), `--clean-mode
strict`, `PYTHONHASHSEED=0`.

Datasets: `ci` = the four real *M. tuberculosis* draft assemblies from Panaroo's own CI
release (`panaroo_test_data` v0.0.2; ~5100 CDS across ~35 contigs each). `smoke` = the first
three contigs of each of those.

---

<a name="b1"></a>
## B1 — `find_missing.py:61`: `mem` is read from a leaked loop variable

### The code

```python
    # identify nodes that have been merged at the protein level
    merged_ids = {}
    for node in G.nodes():
        if (len(G.nodes[node]['centroid']) > 1) or (G.nodes[node]['mergedDNA']):
            for sid in sorted(G.nodes[node]['seqIDs']):     # <- sid bound here
                merged_ids[sid] = node

    merged_nodes = defaultdict(dict)
    with open(gene_data_file, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.split(",")
            if line[2] in merged_ids:
                mem = int(sid.split("_")[0])               # <- and read here, line 61
                if merged_ids[line[2]] in merged_nodes[mem]:
                    merged_nodes[mem][merged_ids[line[2]]] = G.nodes[
                        merged_ids[line[2]]]["dna"][G.nodes[merged_ids[
                            line[2]]]['maxLenId']]
                else:
                    merged_nodes[mem][merged_ids[line[2]]] = line[5]
```

`sid` is not defined in the second loop. It survives from the loop that ends at line 53, so
it holds the lexicographically last seqID of the last node that satisfied the merged test.
`mem` is therefore **one arbitrary, constant genome index for every row of
`gene_data.csv`**.

### Consequence

All merged-node DNA piles into a single `merged_nodes[mem]` bucket. Every other genome gets
an empty dict from the `defaultdict`.

`merged_nodes` is consumed one genome at a time:

```python
delayed(search_gff)(search_list[member], conflicts[member], gff_handle,
                    merged_nodes=merged_nodes[member], ...)
```

and inside `search_gff` it gates the recovery of a merged gene's true extent:

```python
if node in merged_nodes:
    # re-align to find the full merged gene in a +/- search_radius window
    node_locs[node] = [gene[0], loc]
else:
    node_locs[node] = [gene[0], [start - 1, end]]   # just the annotated fragment
```

When several genes are merged into one node, the gene annotated in a given genome is only a
*fragment* of the merged gene, so its true extent has to be recovered by alignment.
`node_locs` then drives the overlap check in `find_missing`, where accumulated
`seq_coverage` decides whether to call `remove_member_from_node`.

**So for every genome except one, that recovery silently never runs.** Merged nodes keep
their short annotated bounds, less coverage is marked, fewer overlaps are detected, and
members Panaroo intended to drop survive into `gene_presence_absence.csv`.

For the one genome that does get the bucket, the recovery runs but with entries built from
every genome's rows: a node appearing in ≥2 rows anywhere ends up with the node's
`maxLenId` DNA rather than that genome's own sequence. That happens to be close to what the
correct code produces, which is plausibly why this went unnoticed.

### Why it reads as unintentional

The strongest argument is not the leaked variable, it's that **the data structure
contradicts the assignment**. `merged_nodes` is a `defaultdict(dict)` keyed by member and
consumed as `merged_nodes[member]` per genome — a shape that is pointless unless `mem`
varied. A flat `dict` would have been written instead.

Supporting: every other read in the block uses `line[2]`, and `line[2]` is exactly the
clustering ID whose first field is the genome index. The same idiom appears at line 77
(`member = int(sid.split("_")[0])`). The intended expression is almost certainly
`int(line[2].split("_")[0])` — a one-token edit.

Cutting slightly the other way: the author wrote `sorted(...)` in the loop just above, so
they were being deliberate about determinism right there.

### Fix

```diff
-                mem = int(sid.split("_")[0])
+                mem = int(line[2].split("_")[0])
```

`tests/parity/reference/patches/b01-merged-nodes-mem.patch`.

### Measured effect

Stock → patched:

| | `smoke` | `ci` |
|---|---|---|
| total gene clusters | 3732 → 3731 | 5118 → 5116 |
| core genes (99–100%) | 746 → 746 | 5097 → 5096 |
| shell genes (15–95%) | 2986 → 2985 | 21 → 20 |
| presence calls | 9838 → 9837 | 20431 → 20428 |
| refound genes | 18 → 18 | 56 → 56 |
| wall clock | 53.9 s → 53.1 s | 65.3 s → 66.0 s |

`gene_presence_absence.Rtab` differs, so this is a structural change to the pangenome, not a
rendering difference. `pre_filt_graph.gml` is byte-identical on both datasets, as expected —
`find_missing` runs after it is written. On `ci`, `gene_data.csv` and `combined_*_CDS.fasta`
also differ, because `find_missing` appends refound records to them.

Direction matches the analysis: more coverage marked → more overlaps detected → a node loses
its last member and is deleted. The magnitude is small but real: 1–2 clusters out of
3700–5100.

Runtime is unchanged. I had predicted a material slowdown, on the reasoning that the fix
un-masks a search that is currently skipped for all but one genome. That did not show up at
either scale (<1% either way). Worth re-measuring at higher genome counts before assuming it
stays free.

---

<a name="b2"></a>
## B2 — `find_missing.py:329`: validity check reads loop leftovers

### The code

```python
    for node in node_search_dict:
        best_hit = ""
        best_loc = None
        for search in node_search_dict[node]:
            ...
            hit, loc = search_dna(db_seq, search[0], prop_match,
                                  pairwise_id_thresh, refind=True)
            ...
            if len(hit) > len(best_hit):
                best_hit = hit
                best_loc = [gene[0], loc]

        if only_valid_genes:
            if not is_valid_gene(hit, translate(search[0])):   # <- line 330
                continue

        hits.append((node, best_hit))
```

The `if only_valid_genes` block sits at the outer loop level, so `hit` and `search` are
whatever the **last** inner iteration left behind, not the best candidate the loop selected.

Two distinct problems:

1. **Wrong candidate.** The node is accepted or rejected based on the last search's `hit`,
   even though `best_hit` is what actually gets appended. A node whose best hit is excellent
   is dropped if its last hit happened to be poor, and vice versa.

2. **Wrong protein.** `translate(search[0])` translates the *query* — `search[0]` is
   `G.nodes[node]["dna"][maxLenId]`, the node's own representative sequence — not the hit.
   So `is_valid_gene` is asked whether the thing we were searching *for* is a valid gene,
   which it essentially always is, since it came from a real annotation. The check is close
   to inert in the accept direction; its real effect is the length-divisible-by-3 test on
   `hit`, applied to the wrong `hit`.

### Scope

`only_valid_genes` is `True` only under `--refind-mode strict`:

```python
if args.refind_mode == "strict": only_valid_genes = True
else:                            only_valid_genes = False
```

The default is `--refind-mode default`, so **this bug does not affect a default run.**

### Correction to an earlier claim

I initially wrote that the `continue` desynchronises indices, because it skips
`hits.append` while a later loop indexes `hits_trans_dict[member][i]` positionally. That is
wrong. Both `hits_trans_dict[member]` and the consuming loop iterate the same `hits` list in
the same order:

```python
hits_trans_dict[member] = Parallel(...)(delayed(translate_to_match)(hit[1], ...) for hit in hits)
...
for node, dna_hit in hits:
    i += 1
    hit_protein = hits_trans_dict[member][i]
```

A skipped node is absent from `hits` entirely, so the two stay aligned. The indexing is
fine; only the candidate selection and the protein argument are wrong.

### Suggested fix

Not written — it needs a decision about intent that belongs to the Panaroo authors, since
`best_hit`'s own protein is not computed at that point in the function.

Status: **read only**, not measured, no patch. Would need a `--refind-mode strict` run to
quantify.

---

<a name="b3"></a>
## B3 — `cdhit.py:535`: `pqid` typo, provably no effect

```python
    if dna:
        pwid = 0.0
        for sA in [seqA, str(Seq(seqA).reverse_complement())]:
            aln = edlib.align(sA, seqB, mode="HW", task='distance',
                              k=0.5 * len(seqA), ...)
            if aln['editDistance'] == -1:
                pqid = max(pwid, 0.0)          # <- assigns to pqid, never read
            else:
                pwid = max(pwid, 1.0 - aln['editDistance'] / float(len(seqA)))
```

`pqid` is written and never read; `pwid` is left unchanged on the no-alignment branch.

**This is behaviour-neutral, and I want to be clear about that** because the obvious "fix"
is a no-op too. `pwid` starts at `0.0` and is only ever assigned `max(pwid, ...)`, so it is
monotonically non-decreasing and never negative. The intended `pwid = max(pwid, 0.0)` is
therefore also a no-op. (The `else` branch cannot produce a negative either: edlib is called
with `k = 0.5 * len(seqA)`, so any reported `editDistance` satisfies
`1 - editDistance/len(seqA) >= 0.5`.)

So: a real typo, zero observable consequence today. It is a latent landmine — if the
initialisation or the loop ever changed, the branch would start silently doing nothing — and
it is a trap for anyone porting the function, who may "correct" it into something that also
does nothing and conclude the port is faithful for the wrong reason.

Status: **read only**, no patch, nothing to measure.

---

<a name="b4"></a>
## B4 — `prokka.py:171`: contig scan never breaks

```python
    for entry in parsed_gff.all_features(featuretype=()):
        if "CDS" not in entry.featuretype:
            continue
        scaffold_id = None
        gene_sequence = None
        for sequence_index in range(len(sequences)):
            scaffold_id = sequences[sequence_index].id
            if scaffold_id == entry.seqid:
                ...
                scaffold_genes[scaffold_id].append(gene_record)
                # no break
```

Every CDS scans every contig, and keeps scanning after finding its match.

- **Performance:** O(n_genes × n_contigs) per genome. On a fragmented draft assembly — 5127
  CDS across 35 contigs in the `ci` data, and far more contigs in a poor assembly — this is
  a meaningful share of input parsing.
- **Latent correctness:** if two records in the `##FASTA` block shared an ID, the gene would
  be appended twice. Prokka does not emit duplicate contig names, so this does not fire in
  practice — but Panaroo accepts GFF3 from other sources, including a user-supplied
  `gff<TAB>fasta` pair, where the guarantee does not hold.

Note the missing `break` is also what makes the `if gene_sequence is None:` check below the
loop reachable at all, so adding a `break` requires restructuring rather than a one-line
edit.

Status: **read only**. Not a Tier B candidate — it is primarily a performance issue, tracked
in `PORTING_PLAN.md` §9 item 3.

---

<a name="b5"></a>
## B5 — `__main__.py:190`: help text contradicts the default

```python
        "--length_outlier_support_proportion",
        help=("... (default=0.01). Genes failing this test will be re-annotated ..."),
        type=float,
        default=0.1)
```

The help says `0.01`; the actual default is `0.1`, a factor of ten. The parameter controls
when a length-outlier gene is re-annotated at the shorter length, so a user following the
documentation gets ten times more aggressive re-annotation than they expect.

Documentation-only — but the kind that silently changes results for anyone who reads the
help and decides the default is fine.

---

<a name="b6"></a>
## B6 — output is not reproducible run to run

Not a coding mistake so much as an unintended property, but it has the same practical effect
and is worth reporting.

Panaroo iterates Python `set`s at several points where the iteration order reaches the
output. For `set[str]` — the `seqIDs` sets — CPython's ordering depends on the
per-process hash seed, so **two runs of Panaroo on identical input produce different files**
unless `PYTHONHASHSEED` is pinned.

### Measured

Stock upstream on `smoke` under `PYTHONHASHSEED` 0, 1 and 2:

| file | across seeds |
|---|---|
| `gene_data.csv`, `combined_{DNA,protein}_CDS.fasta`, `combined_protein_cdhit_out.txt{,.clstr}` | byte-identical |
| `pan_genome_reference.fa`, `gene_presence_absence.Rtab`, `struct_presence_absence.Rtab`, `summary_statistics.txt` | byte-identical |
| `pre_filt_graph.gml`, `final_graph.gml` | differ — **only** in `seqIDs`/`geneIDs` element order; every other key, node and edge identical |
| `gene_presence_absence{,_roary}.csv` | differ — only in `;`-token order within cells, and row order |

The affected sites, as finally audited (the first pass found ten; the rest surfaced later,
each caught by a specific test rather than by re-reading):

| site | container | reaches |
|---|---|---|
| `__main__.py:394,553` `";".join(seqIDs)` | `set[str]` | `geneIDs` in both GML files |
| `__main__.py:556` `seqIDs = list(seqIDs)` | `set[str]` | `final_graph.gml` |
| `generate_output.py:108` `for seq in seqIDs` | `set[str]` | presence/absence cell contents |
| `generate_output.py:68` `for node in component` | `set[int]` | `entry_count` → row order |
| `isvalid.py:172` `custom_stringizer`'s `set` branch | `set[str]` | `pre_filt_graph.gml`, where `seqIDs` is still a set |
| `clean_network.py:170` `list(search_space)` | `set[int]` | **structural** — seeds the BFS in `collapse_families` |
| `clean_network.py:431,467` `merge_node_cluster(G, <set>, …)` | `set[int]` | **structural** — feeds a stable `sorted(key=size)` |
| `clean_network.py:77` `list(set(labels[…]))` | `set[int]` | node-id assignment |
| `merge_nodes.py:126,145` `list(set([...]))` | `set[int]` | **structural** — adjacency insertion order |
| `find_missing.py:150,168` `for node in bad_nodes` | `set[int]` | **structural** — `delete_node` synthesises edges |
| `generate_output.py:405` `os.listdir(alignments_dir)` | **filesystem order** | column order of `core_gene_alignment.aln` |
| `generate_output.py:436,458` `for iso in isolates` | `set[str]` | row order of `core_gene_alignment.aln` |
| `generate_alignments.py:855,885` `os.listdir(...)` | **filesystem order** | the returned alignment list |

The two `os.listdir` sites are arguably worse than the sets: directory order is not stable
across machines or filesystems, so the same input on two computers can produce different
column order in the core alignment.

Two of these were missed on a first audit and found only by testing:

- `isvalid.py:172` — surfaced when the cross-seed acceptance test still failed on
  `pre_filt_graph.gml` after twelve of thirteen files had been fixed. `seqIDs` is still a
  Python `set` at that point in the run, so the GML stringizer renders it in set order.
- the alignment-path sites — invisible to a default-pipeline test, because they are only
  reached with `--alignment`.

Good news for anyone relying on Panaroo's science: `gene_presence_absence.Rtab` is
byte-stable across seeds, so **the pangenome partition itself does not move** — only its
rendering. All four differing files become identical under
`tests/parity/canonicalise.py`.

### Caveat

`PYTHONHASHSEED` perturbs `str` hashing only. Three further sites iterate `set[int]`, whose
order CPython never randomises, and those *are* structurally load-bearing:

| site | effect |
|---|---|
| `clean_network.py:170` `temp_node_list = list(search_space)` | seeds the BFS in `collapse_families`, so it decides which merges are attempted |
| `clean_network.py:431` `merge_node_cluster(G, list(cluster_dict[c]), …)` | feeds a stable `sorted(key=size)`, so ties resolve by set order |
| `find_missing.py:150,168` `for node in bad_nodes: delete_node(...)` | `delete_node` synthesises edges, so deletion order changes the edge set |

These are reproducible for a given input (no randomisation), so this is not a
non-determinism bug. But it does mean the experiment above says nothing about whether
`collapse_families` is confluent — testing that needs a patched build that sorts those
three, compared against stock. Not yet done.

### Suggested fix upstream

Insert `sorted()` at the sites above, and `sorted(os.listdir(...))` at the two directory
listings. Panaroo already does this in two places (`find_missing.py:52,77`), so the pattern
is established; it just is not applied consistently. The complete patch is
`tests/parity/reference/patches/d01-determinism.patch` in this repository — 22 one-line
edits across seven files — and `tests/parity/determinism.sh` verifies the result is
byte-identical across `PYTHONHASHSEED` values.

---

## Reproducing any of this

```sh
tests/parity/env/create.sh && conda activate panaroo-parity
tests/parity/data/fetch.sh ci
tests/parity/data/fetch.sh smoke 3
tests/parity/reference/apply.sh                       # stock (enabled.txt has Tier B off)

# with B1 applied, into a separate tree
tests/parity/reference/apply.sh -o /tmp/ref-b01 \
    -w "$PWD/tests/parity/reference/patches/b01-merged-nodes-mem.patch"

IN=$PWD/tests/parity/build/data/smoke/input.txt
(cd tests/parity/build/reference && PYTHONHASHSEED=0 python -m panaroo -i $IN -o /tmp/out_stock --clean-mode strict -t 8)
(cd /tmp/ref-b01                 && PYTHONHASHSEED=0 python -m panaroo -i $IN -o /tmp/out_b01   --clean-mode strict -t 8)

python3 tests/parity/canonicalise.py /tmp/out_stock /tmp/out_b01
diff /tmp/out_stock/summary_statistics.txt /tmp/out_b01/summary_statistics.txt
```

For B6, run the same input twice with different `PYTHONHASHSEED` values instead.

## Status in this port

Per `PORTING_PLAN.md` §6.1, all Tier B fixes are **off by default** while the translation is
in progress: the reference must stay as close to upstream as possible so that any parity
diff points at a translation error rather than a behaviour change we introduced. Enabling
B1 is a scientific decision to be made separately, and doing so makes `panaroo-rs` disagree
with published Panaroo results.

B6 is different: `PORTING_PLAN.md` §6.1 Tier D applies the `sorted()` fixes unconditionally,
because without them there is no stable reference to compare against at all.
