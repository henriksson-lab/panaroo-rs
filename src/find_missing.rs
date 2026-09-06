//! Translation of `panaroo/panaroo/find_missing.py` — the gene refinding stage.

use crate::support::graph::Graph;
use crate::support::pydict::PyDict;
use std::collections::{BTreeSet, HashMap};

/// `find_missing.py::blosum50`
///
/// Declared in the Python module but never read by any function in it. Kept so
/// `ccc-rs constants-diff` does not report a divergence.
pub const BLOSUM50: [[i8; 25]; 25] = [[0; 25]; 25];

/// What one `search_gff` call returns: `[hits, node_locs, max_seq_len]`.
pub struct SearchGffResult {
    /// `hits` — `(node, dna_hit)` pairs, in `node_search_dict` iteration order.
    pub hits: Vec<(usize, String)>,
    /// `node_locs` — node -> `(contig_id, loc)`. `loc` is `[start, end]` or
    /// `[start, end, strand_flag]`; see [`search_dna`].
    pub node_locs: PyDict<usize, (String, Vec<i64>)>,
    /// `max_seq_len` — longest contig in this genome.
    pub max_seq_len: usize,
}

/// `find_missing.py::find_missing`
///
/// Builds a per-genome search list of accessory genes that a neighbour has but this genome
/// lacks, hands each genome to [`search_gff`], then folds the hits back into the graph and
/// appends them to `combined_DNA_CDS.fasta`, `combined_protein_CDS.fasta` and
/// `gene_data.csv` as `{member}_refound_{n}` records.
///
/// // UPSTREAM BUG (`find_missing.py:61`): in the `merged_nodes` loop,
/// // `mem = int(sid.split("_")[0])` reads `sid`, which leaks from the loop ending at line
/// // 53, so `mem` is one arbitrary constant genome index for every row. All merged-node DNA
/// // lands in a single bucket and every other genome gets an empty dict, disabling the
/// // merged-gene extent recovery in `search_gff` for them. Reproduce exactly.
/// // See ORIGINAL_CODE_BUG.md B1 — measured at 1-2 gene clusters out of 3700-5100.
///
/// PORTING_PLAN.md §9 item 9: a fresh joblib pool per genome for the translation step.
#[allow(clippy::too_many_arguments)]
pub fn find_missing(
    g: &mut Graph,
    gff_file_handles: &[String],
    dna_seq_file: &str,
    prot_seq_file: &str,
    gene_data_file: &str,
    merge_id_thresh: f64,
    search_radius: i64,
    prop_match: f64,
    pairwise_id_thresh: f64,
    n_cpu: i64,
    remove_by_consensus: bool,
    only_valid_genes: bool,
    verbose: bool,
) {
    use crate::merge_nodes::{delete_node, remove_member_from_node};
    use crate::support::pydict::PyDict as Dict;
    use std::collections::{HashMap as Map, HashSet};
    use std::io::Write;

    // generate mapping between internal nodes and gff ids
    let gene_data = std::fs::read_to_string(gene_data_file)
        .unwrap_or_else(|e| panic!("could not read {gene_data_file}: {e}"));
    let mut id_to_gff: Map<String, String> = Map::new();
    for line in gene_data.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() < 4 {
            continue;
        }
        if id_to_gff.contains_key(f[2]) {
            panic!("NameError: Duplicate internal ids!");
        }
        id_to_gff.insert(f[2].to_string(), f[3].to_string());
    }

    // identify nodes that have been merged at the protein level
    let mut merged_ids: Map<String, usize> = Map::new();
    let mut last_sid: Option<String> = None;
    for node in g.nodes() {
        if g.node(node).centroid.len() > 1 || g.node(node).merged_dna {
            for sid in g.node(node).seq_ids.iter() {
                merged_ids.insert(sid.clone(), node);
                last_sid = Some(sid.clone());
            }
        }
    }

    // UPSTREAM BUG (find_missing.py:61): `mem` is computed from `sid`, which leaks from the
    // loop above -- so it is one arbitrary constant genome index for every row, not the
    // row's own genome. All merged-node DNA lands in a single bucket and every other genome
    // gets an empty dict, which disables the merged-gene extent recovery in search_gff for
    // them. Tier B patch b01 fixes it; with the patch enabled the reference computes
    // `int(line[2].split("_")[0])` instead. See ORIGINAL_CODE_BUG.md B1.
    //
    // The reference this port targets has b01 ENABLED (PORTING_PLAN.md §6.1), so the fixed
    // form is what we translate. `last_sid` is kept only to document the original.
    let _ = &last_sid;
    let mut merged_nodes: Map<usize, Map<usize, String>> = Map::new();
    for line in gene_data.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() < 6 {
            continue;
        }
        if let Some(&node) = merged_ids.get(f[2]) {
            let mem: usize = f[2].split('_').next().unwrap().parse().unwrap();
            let bucket = merged_nodes.entry(mem).or_default();
            if bucket.contains_key(&node) {
                let n = g.node(node);
                bucket.insert(node, n.dna[n.max_len_id].clone());
            } else {
                bucket.insert(node, f[5].to_string());
            }
        }
    }

    // iterate through nodes to identify accessory genes for searching: nodes missing a
    // member that at least one neighbour has
    let mut n_searches = 0usize;
    let mut search_list: Map<usize, Dict<usize, BTreeSet<(String, String)>>> = Map::new();
    let mut conflicts: Map<usize, BTreeSet<(usize, String)>> = Map::new();
    for node in g.nodes() {
        for neigh in g.neighbors(node) {
            for sid in g.node(neigh).seq_ids.iter() {
                let member: usize = sid.split('_').next().unwrap().parse().unwrap();
                conflicts
                    .entry(member)
                    .or_default()
                    .insert((neigh, id_to_gff[sid].clone()));
                if !g.node(node).members.contains(member) {
                    let n = g.node(node);
                    if n.dna[n.max_len_id].is_empty() {
                        panic!("NameError: Problem!");
                    }
                    search_list
                        .entry(member)
                        .or_default()
                        .entry_or_default(node)
                        .insert((n.dna[n.max_len_id].clone(), id_to_gff[sid].clone()));
                    n_searches += 1;
                }
            }
        }
    }

    if verbose {
        println!("Number of searches to perform:  {n_searches}");
        println!("Searching...");
    }

    let jobs: Vec<usize> = (0..gff_file_handles.len()).collect();
    let results = crate::support::parallel::parallel_map(n_cpu, jobs, |member| {
        let empty_search: Dict<usize, BTreeSet<(String, String)>> = Dict::new();
        let empty_conflicts: BTreeSet<(usize, String)> = BTreeSet::new();
        let empty_merged: Map<usize, String> = Map::new();
        let sl = search_list.get(&member).unwrap_or(&empty_search);
        let cf = conflicts.get(&member).unwrap_or(&empty_conflicts);
        let mn = merged_nodes.get(&member).unwrap_or(&empty_merged);
        search_gff(
            sl,
            cf,
            &gff_file_handles[member],
            mn,
            search_radius,
            prop_match,
            pairwise_id_thresh,
            merge_id_thresh,
            only_valid_genes,
            1,
        )
    });

    if verbose {
        println!("translating hits...");
    }

    let mut hits_trans_dict: Vec<Vec<String>> = Vec::with_capacity(results.len());
    for r in &results {
        let jobs: Vec<usize> = (0..r.hits.len()).collect();
        let trans = crate::support::parallel::parallel_map(n_cpu, jobs, |i| {
            let (node, hit) = &r.hits[i];
            translate_to_match(hit, &g.node(*node).protein[0])
        });
        hits_trans_dict.push(trans);
    }

    // remove nodes that conflict (overlap)
    let mut nodes_by_size: Vec<(usize, usize)> =
        g.nodes().into_iter().map(|n| (g.node(n).size, n)).collect();
    // `sorted(..., reverse=True)` on (size, node) tuples
    nodes_by_size.sort_by(|a, b| b.cmp(a));
    let nodes_by_size: Vec<usize> = nodes_by_size.into_iter().map(|(_, n)| n).collect();

    let mut bad_node_mem_pairs: HashSet<(usize, usize)> = HashSet::new();
    let mut bad_nodes: BTreeSet<usize> = BTreeSet::new();
    for (member, r) in results.iter().enumerate() {
        let mut seq_coverage: Map<String, Vec<bool>> = Map::new();
        for &node in &nodes_by_size {
            if bad_nodes.contains(&node) {
                continue;
            }
            let Some((contig_id, loc)) = r.node_locs.get(&node) else {
                continue;
            };
            let cov = seq_coverage
                .entry(contig_id.clone())
                .or_insert_with(|| vec![false; r.max_seq_len + 2]);
            let lo = loc[0].max(0) as usize;
            let hi = (loc[1].max(0) as usize).min(cov.len());
            let covered = if lo < hi {
                cov[lo..hi].iter().filter(|&&b| b).count()
            } else {
                0
            };
            let max_len = *g.node(node).lengths.iter().max().unwrap_or(&0);
            if covered as f64 >= 0.5 * max_len as f64 {
                if g.node(node).members.contains(member) {
                    remove_member_from_node(g, node, member);
                }
                bad_node_mem_pairs.insert((node, member));
            } else if lo < hi {
                for b in &mut cov[lo..hi] {
                    *b = true;
                }
            }
        }
    }

    for node in g.nodes() {
        if g.node(node).members.is_empty() {
            bad_nodes.insert(node);
        }
    }
    // Tier D: sorted(bad_nodes) -- delete_node synthesises edges, so order matters.
    for node in bad_nodes.clone() {
        if g.has_node(node) {
            delete_node(g, node);
        }
    }

    // remove by consensus
    if remove_by_consensus {
        if verbose {
            println!("removing by consensus...");
        }
        let mut node_hit_counter: Map<usize, usize> = Map::new();
        for (member, r) in results.iter().enumerate() {
            for (node, dna_hit) in &r.hits {
                if dna_hit.is_empty() || bad_nodes.contains(node) {
                    continue;
                }
                if bad_node_mem_pairs.contains(&(*node, member)) {
                    continue;
                }
                *node_hit_counter.entry(*node).or_insert(0) += 1;
            }
        }
        for node in g.nodes() {
            if node_hit_counter.get(&node).copied().unwrap_or(0) > g.node(node).size {
                bad_nodes.insert(node);
            }
        }
        for node in bad_nodes.clone() {
            if g.has_node(node) {
                delete_node(g, node);
            }
        }
    }

    if verbose {
        println!("Updating output...");
    }

    let mut n_found = 0usize;
    let mut dna_out = std::fs::OpenOptions::new()
        .append(true)
        .open(dna_seq_file)
        .expect("append dna");
    let mut prot_out = std::fs::OpenOptions::new()
        .append(true)
        .open(prot_seq_file)
        .expect("append prot");
    let mut data_out = std::fs::OpenOptions::new()
        .append(true)
        .open(gene_data_file)
        .expect("append csv");

    for (member, r) in results.iter().enumerate() {
        for (i, (node, dna_hit)) in r.hits.iter().enumerate() {
            if dna_hit.is_empty() || bad_nodes.contains(node) {
                continue;
            }
            if bad_node_mem_pairs.contains(&(*node, member)) {
                continue;
            }

            let hit_protein = &hits_trans_dict[member][i];
            let (contig_id, loc) = r.node_locs.get(node).expect("node_locs entry");
            let hit_strand = if loc[2] == 0 { '+' } else { '-' };
            let refound_id = format!("{member}_refound_{n_found}");

            {
                let n = g.node_mut(*node);
                n.members.add(member);
                n.size += 1;
                n.dna.push(dna_hit.clone());
                crate::isvalid::del_dups(&mut n.dna);
                n.protein.push(hit_protein.clone());
                crate::isvalid::del_dups(&mut n.protein);
                n.seq_ids.insert(refound_id.clone());
            }

            writeln!(dna_out, ">{refound_id}\n{dna_hit}").expect("write dna");
            writeln!(prot_out, ">{refound_id}\n{hit_protein}").expect("write prot");
            let gff_name = file_stem(&gff_file_handles[member]);
            writeln!(
                data_out,
                "{gff_name},{contig_id},{refound_id},{refound_id},{hit_protein},{dna_hit},,location:{}-{};strand:{hit_strand}",
                loc[0],
                loc[1],
            )
            .expect("write csv");

            n_found += 1;
        }
    }

    if verbose {
        println!("Number of refound genes:  {n_found}");
    }
}

/// `os.path.splitext(os.path.basename(p))[0]`
fn file_stem(p: &str) -> String {
    let base = p.rsplit('/').next().unwrap_or(p);
    match base.rfind('.') {
        Some(i) if i > 0 => base[..i].to_string(),
        _ => base.to_string(),
    }
}

fn filtered_gff_annotation(text: &str) -> String {
    let mut ann = String::new();
    for line in text.lines() {
        if line.contains("##sequence-region") {
            continue;
        }
        if !ann.is_empty() {
            ann.push('\n');
        }
        ann.push_str(line);
    }
    ann
}

/// `find_missing.py::search_gff`
///
/// // UPSTREAM BUG (`find_missing.py:329`): the `only_valid_genes` check reads `hit` and
/// // `search`, both leftovers from the inner loop, so it validates the *last* candidate
/// // rather than `best_hit`, and translates the *query* (`search[0]` is the node's own
/// // DNA) rather than the hit. Only active under `--refind-mode strict`. Preserve it.
/// // See ORIGINAL_CODE_BUG.md B2. The `continue` does NOT desynchronise indices: a skipped
/// // node is absent from `hits`, and `hits_trans_dict[member]` is built by iterating the
/// // same list, so the positional lookup in `find_missing` stays aligned.
///
/// PORTING_PLAN.md §9 item 5: `contigs` is built as a numpy char array and converted
/// straight back to a `String`, because the masking write that justified it is commented
/// out. Translate the round trip as a plain `String` — it is observationally identical —
/// but note the deviation.
#[allow(clippy::too_many_arguments)]
pub fn search_gff(
    node_search_dict: &PyDict<usize, BTreeSet<(String, String)>>,
    conflicts: &BTreeSet<(usize, String)>,
    gff_handle_name: &str,
    merged_nodes: &HashMap<usize, String>,
    search_radius: i64,
    prop_match: f64,
    pairwise_id_thresh: f64,
    merge_id_thresh: f64,
    only_valid_genes: bool,
    _n_cpu: i64,
) -> SearchGffResult {
    use crate::support::pydict::PyDict as Dict;
    use std::collections::HashSet;

    let raw = std::fs::read_to_string(gff_handle_name)
        .unwrap_or_else(|e| panic!("could not read {gff_handle_name}: {e}"));
    let text = raw.replace(',', "");
    let split: Vec<&str> = text.split("##FASTA\n").collect();
    let mut node_locs: Dict<usize, (String, Vec<i64>)> = Dict::new();

    if split.len() != 2 {
        panic!("NameError: File does not appear to be in GFF3 format!");
    }

    // load fasta
    //
    // PORTING_PLAN.md §9 item 5: the Python stores each contig as a numpy char array and
    // converts it straight back to a string below, because the masking write that justified
    // the array form is commented out (find_missing.py:296). Kept as a String -- that is
    // observationally identical, and the deviation is noted rather than hidden.
    let mut contigs: Dict<String, String> = Dict::new();
    let mut max_seq_len = 0usize;
    for record in crate::support::seqio::parse_fasta(split[1]) {
        max_seq_len = max_seq_len.max(record.seq.len());
        contigs.insert(record.id, record.seq);
    }

    // load gff annotation
    let ann = filtered_gff_annotation(split[0]);
    let parsed_gff = crate::support::gff::GffDb::create_db(&ann)
        .unwrap_or_else(|e| panic!("NameError: File does not appear to be in GFF3 format! {e}"));

    // mask regions that already have genes and convert back to string
    for (node, geneid) in conflicts {
        let gene = parsed_gff.get(geneid);
        let start = gene.start.min(gene.stop);
        let end = gene.start.max(gene.stop);

        if let Some(merged) = merged_nodes.get(node) {
            let contig = contigs
                .get(&gene.seqid)
                .unwrap_or_else(|| panic!("KeyError: {}", gene.seqid));
            let lo = 0.max(start - search_radius) as usize;
            let hi = ((end + search_radius) as usize).min(contig.len());
            let db_seq = if lo < hi {
                contig[lo..hi].to_string()
            } else {
                String::new()
            };

            let (_hit, mut loc) = search_dna(
                &db_seq,
                merged,
                (end - start) as f64 / merged.len() as f64,
                merge_id_thresh,
                false,
            );

            // update location
            loc[0] += 0.max(start - search_radius);
            loc[1] += 0.max(start - search_radius);
            node_locs.insert(*node, (gene.seqid.clone(), loc));
        } else {
            node_locs.insert(*node, (gene.seqid.clone(), vec![start - 1, end]));
        }
    }

    // Duplicate-entry check. The masking write it once guarded is commented out upstream,
    // so this loop now only raises on duplicates.
    let mut seen: HashSet<(String, i64, i64)> = HashSet::new();
    for (_node, geneid) in conflicts {
        let gene = parsed_gff.get(geneid);
        let start = gene.start.min(gene.stop);
        let end = gene.start.max(gene.stop);
        if !seen.insert((gene.seqid.clone(), start - 1, end)) {
            panic!("NameError: Duplicate entry!!!");
        }
    }

    // search for matches
    let mut hits: Vec<(usize, String)> = Vec::new();
    for (node, searches) in node_search_dict.items() {
        let mut best_hit = String::new();
        let mut best_loc: Option<(String, Vec<i64>)> = None;
        // `hit` and `search` leak out of this loop and are read below -- see the bug note.
        let mut last_hit = String::new();
        let mut last_search: Option<&(String, String)> = None;

        for search in searches {
            let gene = parsed_gff.get(&search.1);
            let start = gene.start.min(gene.stop);
            let end = gene.start.max(gene.stop);
            let contig = contigs
                .get(&gene.seqid)
                .unwrap_or_else(|| panic!("KeyError: {}", gene.seqid));
            let lo = 0.max(start - search_radius) as usize;
            let hi = ((end + search_radius) as usize).min(contig.len());
            let db_seq = if lo < hi {
                contig[lo..hi].to_string()
            } else {
                String::new()
            };

            let (hit, mut loc) =
                search_dna(&db_seq, &search.0, prop_match, pairwise_id_thresh, true);
            // update location
            loc[0] += 0.max(start - search_radius);
            loc[1] += 0.max(start - search_radius);

            if hit.len() > best_hit.len() {
                best_hit = hit.clone();
                best_loc = Some((gene.seqid.clone(), loc));
            }
            last_hit = hit;
            last_search = Some(search);
        }

        if only_valid_genes {
            // UPSTREAM BUG (find_missing.py:329): this reads `hit` and `search`, the LAST
            // inner-loop values, not `best_hit` and the search that produced it -- and it
            // translates the *query* (`search[0]` is the node's own DNA) rather than the
            // hit. Only active under --refind-mode strict. See ORIGINAL_CODE_BUG.md B2.
            // The `continue` does not desynchronise anything: the node is simply absent
            // from `hits`, and hits_trans_dict is built by iterating the same list.
            if let Some(sr) = last_search {
                if !crate::isvalid::is_valid_gene(&last_hit, &crate::support::seq::translate(&sr.0))
                {
                    continue;
                }
            }
        }

        hits.push((*node, best_hit.clone()));
        if let Some(bl) = best_loc {
            if !best_hit.is_empty() {
                node_locs.insert(*node, bl);
            }
        }
    }

    SearchGffResult {
        hits,
        node_locs,
        max_seq_len,
    }
}

/// `find_missing.py::repl` — the `re.sub` replacement callback; returns `'X' * len(match)`.
pub fn repl(matched: &str) -> String {
    "X".repeat(matched.len())
}

/// `find_missing.py::search_dna`
///
/// Returns `(seq, loc)`. **`loc` is a variable-length list**: it starts as `[0, 0]` and is
/// only replaced by the three-element `[start, end, strand_flag]` when a hit is accepted.
/// So a miss returns two elements, and `find_missing`'s `node_locs[node][1][2]` strand read
/// would raise `IndexError` on it — which is why that read is guarded by `best_hit != ""`.
/// Modelled as a `Vec<i64>` rather than `[i64; 3]` to keep that distinction.
///
/// Pads the database sequence with `len(search_sequence) / 2` `E` characters at each end and
/// declares `E` equal to every base, so a gene running off a contig end still aligns. Runs
/// edlib in `HW` mode on both strands, then picks the hit whose location is closest to the
/// centre.
///
/// Note `min(aln['locations'], key=lambda x: min(centre - x[0], centre - x[1]))` uses
/// **signed** differences, not absolute values. Preserve that; it is load-bearing.
pub fn search_dna(
    db_seq: &str,
    search_sequence: &str,
    prop_match: f64,
    pairwise_id_thresh: f64,
    _refind: bool,
) -> (String, Vec<i64>) {
    use std::borrow::Cow;

    use crate::support::edlib::{align, Mode, Task};
    use crate::support::seq::reverse_complement;

    let mut found_dna = String::new();
    let mut max_hit = 0.0f64;
    let mut loc: Vec<i64> = vec![0, 0];

    let added_e_len = search_sequence.len() / 2;
    let rc = reverse_complement(db_seq);

    for (i, db_raw) in [db_seq, rc.as_str()].into_iter().enumerate() {
        // add some Ns at the start and end to deal with fragments at the end of contigs
        let mut db = String::with_capacity(db_raw.len() + 2 * added_e_len);
        db.extend(std::iter::repeat('E').take(added_e_len));
        db.push_str(db_raw);
        db.extend(std::iter::repeat('E').take(added_e_len));

        let aln = align(
            search_sequence,
            &db,
            Mode::Hw,
            Task::Path,
            10 * search_sequence.len() as i64,
            &E_AND_N_EQUALITIES,
        );

        // remove terminal inserts
        let edit_distance =
            aln.edit_distance - terminal_insert_adjustment(aln.cigar.as_deref().unwrap_or(""));

        let (start, end, tloc);
        if edit_distance == -1 || aln.locations.is_empty() {
            continue;
        } else {
            // take hit that is closest to the centre of the neighbouring gene.
            // NOTE: `min(..., key=lambda x: min(centre - x[0], centre - x[1]))` uses SIGNED
            // differences, not absolute values. Preserved -- it is load-bearing.
            let centre = db.len() as f64 / 2.0;
            let mut best = aln.locations[0];
            let key = |l: (Option<i64>, i64)| -> f64 {
                let a = centre - l.0.unwrap_or(0) as f64;
                let b = centre - l.1 as f64;
                a.min(b)
            };
            let mut best_k = key(best);
            for &l in &aln.locations[1..] {
                let k = key(l);
                if k < best_k {
                    best_k = k;
                    best = l;
                }
            }
            tloc = best;
            start = tloc.0.unwrap_or(0);
            end = tloc.1 + 1;
        }

        let mut possible_dbs: Vec<Cow<'_, str>> = Vec::with_capacity(3);
        possible_dbs.push(Cow::Borrowed(db.as_str()));
        if db.contains("NNNNNNNNNNNNNNNNNNNN") {
            possible_dbs.push(Cow::Owned(sub_leading_n_run(&db)));
            possible_dbs.push(Cow::Owned(sub_trailing_n_run(&db)));
        }

        for posdb in possible_dbs {
            let posdb = posdb.as_ref();
            let seg = &posdb[start as usize..(end as usize).min(posdb.len())];
            let (n_x, n_e, acgt) = segment_counts(seg);
            let n_x = n_x as f64;
            let n_e = n_e as f64;

            let aln_length = (end - start) as f64 - n_x - n_e;
            if aln_length / search_sequence.len() as f64 <= prop_match {
                continue;
            }
            if acgt as f64 / search_sequence.len() as f64 <= prop_match {
                continue;
            }

            // determine an approximate percentage identity
            let pid = 1.0 - (edit_distance as f64 - n_x) / aln_length;
            if pid <= pairwise_id_thresh {
                continue;
            }

            if max_hit < pid * aln_length {
                found_dna = seg.to_string();
                max_hit = pid * aln_length;
                let l = if i == 0 {
                    [start, end]
                } else {
                    [
                        posdb.len() as i64 - tloc.1 - 1,
                        posdb.len() as i64 - tloc.0.unwrap_or(0),
                    ]
                };
                loc = vec![
                    0.max(l[0].min(l[1]) - added_e_len as i64),
                    (l[0].max(l[1]) - added_e_len as i64).min(db_seq.len() as i64),
                    i as i64,
                ];
            }
        }
    }

    (normalise_found_dna(&found_dna), loc)
}

/// `additionalEqualities` for `search_dna`: `N` and the padding character `E` both match
/// any base, so a gene running off a contig end still aligns.
const E_AND_N_EQUALITIES: [(char, char); 8] = [
    ('A', 'N'),
    ('C', 'N'),
    ('G', 'N'),
    ('T', 'N'),
    ('A', 'E'),
    ('C', 'E'),
    ('G', 'E'),
    ('T', 'E'),
];

fn normalise_found_dna(seq: &str) -> String {
    let bytes = seq.as_bytes();
    let mut start = 0usize;
    while start < bytes.len() && matches!(bytes[start], b'N' | b'X' | b'E') {
        start += 1;
    }
    let mut end = bytes.len();
    while end > start && matches!(bytes[end - 1], b'N' | b'X' | b'E') {
        end -= 1;
    }
    let mut out = String::with_capacity(end - start);
    for &b in &bytes[start..end] {
        match b {
            b'X' | b'E' => out.push('N'),
            _ => out.push(b as char),
        }
    }
    out
}

fn segment_counts(seg: &str) -> (usize, usize, usize) {
    let mut n_x = 0usize;
    let mut n_e = 0usize;
    let mut acgt = 0usize;
    for b in seg.bytes() {
        match b {
            b'X' => n_x += 1,
            b'E' => n_e += 1,
            b'A' | b'C' | b'G' | b'T' => acgt += 1,
            _ => {}
        }
    }
    (n_x, n_e, acgt)
}

fn terminal_insert_adjustment(cigar: &str) -> i64 {
    let bytes = cigar.as_bytes();
    let mut i = 0usize;
    let mut leading_insert = 0i64;
    let mut trailing_insert = 0i64;
    let mut first_op = true;

    while i < bytes.len() {
        if !bytes[i].is_ascii_digit() {
            i += 1;
            continue;
        }
        let mut n = 0i64;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            n = n * 10 + (bytes[i] - b'0') as i64;
            i += 1;
        }
        let op = bytes.get(i).copied();
        if first_op {
            first_op = false;
            if op == Some(b'I') {
                leading_insert = n;
            }
        }
        trailing_insert = if op == Some(b'I') && i + 1 == bytes.len() {
            n
        } else {
            0
        };
        if op.is_some() {
            i += 1;
        }
    }

    leading_insert + trailing_insert
}

/// `re.split(r'(\d+)', cigar)[1:]`
///
/// Python's `re.split` with a **capturing group** returns
/// `[text_before_first_match, group1, text_between, group2, ..., text_after_last_match]`.
/// Both the first and last elements are surrounding text, and either can be `""`:
///
/// ```text
/// "3374=100I" -> ['', '3374', '=', '100', 'I']   then [1:] -> ['3374', '=', '100', 'I']
/// "4="        -> ['', '4', '=']                  then [1:] -> ['4', '=']
/// "4"         -> ['', '4', '']                   then [1:] -> ['4', '']
/// "=4"        -> ['=', '4', '']                  then [1:] -> ['4', '']
/// ""          -> ['']                            then [1:] -> []
/// ```
///
/// So a trailing `""` appears only when the cigar **ends with a digit**, which a
/// well-formed cigar never does. Getting this backwards silently disables the
/// trailing-insert correction in [`search_dna`] -- `cig[-1] == "I"` never matches -- and
/// that correction is exactly what rescues a gene running off the end of a contig. On the
/// smoke dataset the inverted version cost three refound genes.
#[cfg(test)]
fn split_cigar(cigar: &str) -> Vec<String> {
    let b = cigar.as_bytes();
    let mut full: Vec<String> = Vec::new();
    let mut text = String::new();
    let mut i = 0;
    while i < b.len() {
        if b[i].is_ascii_digit() {
            full.push(std::mem::take(&mut text)); // text before this match
            let start = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            full.push(cigar[start..i].to_string());
        } else {
            text.push(b[i] as char);
            i += 1;
        }
    }
    full.push(text); // text after the last match
    full.remove(0); // the `[1:]`
    full
}

/// `re.sub("^[ACGTEX]{0,}NNNNNNNNNNNNNNNNNNNN", repl, db, 1)`
fn sub_leading_n_run(db: &str) -> String {
    const RUN: &str = "NNNNNNNNNNNNNNNNNNNN";
    let b = db.as_bytes();
    let mut i = 0;
    while i < b.len() && matches!(b[i], b'A' | b'C' | b'G' | b'T' | b'E' | b'X') {
        i += 1;
    }
    // greedy prefix, then backtrack to the first position where the run matches
    let mut j = i;
    loop {
        if db[j..].starts_with(RUN) {
            let end = j + RUN.len();
            return format!("{}{}", repl(&db[..end]), &db[end..]);
        }
        if j == 0 {
            return db.to_string();
        }
        j -= 1;
    }
}

/// `re.sub("NNNNNNNNNNNNNNNNNNNN[ACGTEX]{0,}$", repl, db, 1)`
fn sub_trailing_n_run(db: &str) -> String {
    const RUN: &str = "NNNNNNNNNNNNNNNNNNNN";
    match db.find(RUN) {
        Some(mut start) => {
            // leftmost match wins in re.sub; the tail must be all [ACGTEX] to the end
            loop {
                let after = start + RUN.len();
                if db[after..]
                    .bytes()
                    .all(|c| matches!(c, b'A' | b'C' | b'G' | b'T' | b'E' | b'X'))
                {
                    return format!("{}{}", &db[..start], repl(&db[start..]));
                }
                match db[start + 1..].find(RUN) {
                    Some(k) => start = start + 1 + k,
                    None => return db.to_string(),
                }
            }
        }
        None => db.to_string(),
    }
}

/// `find_missing.py::translate_to_match`
///
/// Translates the hit in all six frames and picks the frame sharing the most 3-mers with
/// the target protein. `max(alignments, key=lambda x: x[1])` returns the **first** maximum
/// (PORTING_PLAN.md §6.4).
///
/// Note the frame list is built as `[... for i in range(3) for s in dna_seqs]`, so the order
/// is `(frame 0, fwd), (frame 0, rev), (frame 1, fwd), ...` — that ordering decides ties.
pub fn translate_to_match(hit: &str, target_prot: &str) -> String {
    use crate::support::seq::reverse_complement;
    use std::collections::HashSet;

    if hit.is_empty() {
        return String::new();
    }

    // translate in all 6 frames splitting on unknown
    let rc = reverse_complement(hit);
    let dna_seqs = [hit, rc.as_str()];

    // `[... for i in range(3) for s in dna_seqs]` -- frame is the OUTER loop, so the order
    // is (frame0,fwd), (frame0,rev), (frame1,fwd), ... That ordering decides ties below.
    let search_set: HashSet<&str> = (0..target_prot.len().saturating_sub(2))
        .map(|i| &target_prot[i..i + 3])
        .collect();

    // `max(alignments, key=lambda x: x[1])` returns the FIRST maximum.
    let mut best = String::new();
    let mut best_n: i64 = -1;
    for i in 0..3 {
        for s in &dna_seqs {
            let target_sequence = translate_padded_frame(s, i);
            let query_set: HashSet<&str> = (0..target_sequence.len().saturating_sub(2))
                .map(|i| &target_sequence[i..i + 3])
                .collect();
            let n = search_set.intersection(&query_set).count() as i64;
            if n > best_n {
                best_n = n;
                best = target_sequence;
            }
        }
    }
    best
}

fn translate_padded_frame(seq: &str, frame: usize) -> String {
    use crate::support::seq::translate;

    let sub = if frame < seq.len() { &seq[frame..] } else { "" };
    // `s[i:].ljust(len + (3 - len % 3), 'N')` -- note when len % 3 == 0 this pads
    // by a further 3 Ns, which translate to X. Preserved.
    let pad = 3 - sub.len() % 3;
    let mut padded = String::with_capacity(sub.len() + pad);
    padded.push_str(sub);
    padded.extend(std::iter::repeat('N').take(pad));
    translate(&padded)
}

#[cfg(test)]
fn translate_to_match_collecting(hit: &str, target_prot: &str) -> String {
    use crate::support::seq::{reverse_complement, translate};
    use std::collections::HashSet;

    if hit.is_empty() {
        return String::new();
    }

    let dna_seqs = [hit.to_string(), reverse_complement(hit)];
    let mut proteins: Vec<String> = Vec::with_capacity(6);
    for i in 0..3 {
        for s in &dna_seqs {
            let sub = if i < s.len() { &s[i..] } else { "" };
            let pad = 3 - sub.len() % 3;
            let mut padded = sub.to_string();
            padded.push_str(&"N".repeat(pad));
            proteins.push(translate(&padded));
        }
    }

    let search_set: HashSet<&str> = (0..target_prot.len().saturating_sub(2))
        .map(|i| &target_prot[i..i + 3])
        .collect();
    let mut best = String::new();
    let mut best_n: i64 = -1;
    for target_sequence in proteins {
        let query_set: HashSet<&str> = (0..target_sequence.len().saturating_sub(2))
            .map(|i| &target_sequence[i..i + 3])
            .collect();
        let n = search_set.intersection(&query_set).count() as i64;
        if n > best_n {
            best_n = n;
            best = target_sequence;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Expected values from CPython's `re.split(r'(\d+)', s)[1:]`.
    #[test]
    fn split_cigar_matches_python_re_split() {
        assert_eq!(split_cigar("3374=100I"), ["3374", "=", "100", "I"]);
        assert_eq!(split_cigar("4="), ["4", "="]);
        assert_eq!(split_cigar("4"), ["4", ""]);
        assert!(split_cigar("").is_empty());
        // `[1:]` drops the leading text, so the "=" is not present
        assert_eq!(split_cigar("=4"), ["4", ""]);
        assert_eq!(split_cigar("10=2X3I"), ["10", "=", "2", "X", "3", "I"]);
    }

    /// The trailing-insert correction fires only when the last token is `"I"`, which needs
    /// the cigar to end with an operator. This is the case that rescues a gene running off
    /// a contig end.
    #[test]
    fn trailing_insert_is_detected() {
        let cig = split_cigar("3374=100I");
        assert_eq!(cig.last().map(|s| s.as_str()), Some("I"));
        assert_eq!(cig[cig.len() - 2], "100");
    }

    #[test]
    fn terminal_insert_adjustment_matches_split_cigar_logic() {
        for cigar in ["3374=100I", "4=", "4", "=4", "", "10I5=2I", "2I3=4I"] {
            let cig = split_cigar(cigar);
            let mut expected = 0i64;
            if cig.last().map(|s| s.as_str()) == Some("I") {
                expected += cig[cig.len() - 2].parse::<i64>().unwrap_or(0);
            }
            if cig.get(1).map(|s| s.as_str()) == Some("I") {
                expected += cig[0].parse::<i64>().unwrap_or(0);
            }
            assert_eq!(terminal_insert_adjustment(cigar), expected, "{cigar}");
        }
    }

    #[test]
    fn normalise_found_dna_matches_replace_then_trim() {
        for seq in ["", "NNN", "EXN", "NEXACXTEGNXE", "ACGT", "XXACEEGTXX"] {
            let expected = seq.replace(['X', 'E'], "N").trim_matches('N').to_string();
            assert_eq!(normalise_found_dna(seq), expected);
        }
    }

    #[test]
    fn segment_counts_match_individual_matches() {
        let seg = "ACGTNXEacgtXXEE";
        assert_eq!(
            segment_counts(seg),
            (
                seg.matches('X').count(),
                seg.matches('E').count(),
                seg.matches('A').count()
                    + seg.matches('C').count()
                    + seg.matches('G').count()
                    + seg.matches('T').count()
            )
        );
    }

    #[test]
    fn filtered_gff_annotation_matches_collect_join() {
        for text in [
            "",
            "##sequence-region ctg 1 10",
            "a\n##sequence-region ctg 1 10\nb\n",
            "a\nb",
            "a\n##sequence-region\n##sequence-region x\nb",
        ] {
            let expected = text
                .lines()
                .filter(|l| !l.contains("##sequence-region"))
                .collect::<Vec<_>>()
                .join("\n");
            assert_eq!(filtered_gff_annotation(text), expected, "{text:?}");
        }
    }

    #[test]
    fn translate_padded_frame_preserves_ljust_padding_rule() {
        assert_eq!(translate_padded_frame("ATGAAA", 0), "MKX");
        assert_eq!(translate_padded_frame("ATGAA", 0), "MX");
        assert_eq!(translate_padded_frame("ATGAAA", 1), "*X");
        assert_eq!(translate_padded_frame("ATGAAA", 3), "KX");
    }

    #[test]
    fn translate_to_match_matches_collect_then_first_max() {
        for (hit, target) in [
            ("ATGAAATAA", "MK*"),
            ("TTTATGAAATAA", "MK*"),
            ("ATGCCCAAATTT", "PF"),
            ("NNNATGAAACCC", "MKP"),
            ("ATGAAA", ""),
        ] {
            assert_eq!(
                translate_to_match(hit, target),
                translate_to_match_collecting(hit, target),
                "{hit} {target}"
            );
        }
    }
}
