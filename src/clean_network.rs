//! Translation of `panaroo/panaroo/clean_network.py`.
//!
//! This is the heart of the pipeline and the most parity-sensitive module: every merge
//! decision here changes node IDs, which appear in both GML files and drive the row order
//! of `gene_presence_absence.csv`. Translate the control flow literally, including the
//! order in which `search_space` is drained.

use crate::support::graph::Graph;
use crate::support::intbitset::IntBitSet;
use crate::support::sparse::CsrMatrix;
use std::collections::HashMap;

/// `clean_network.py::trim_low_support_trailing_ends`
///
/// Note the `removed` flag is set inside the removal loop, so a pass that finds no bad
/// nodes breaks out — but a pass that finds some always continues to the next iteration.
pub fn trim_low_support_trailing_ends(g: &mut Graph, min_support: usize, max_recursive: usize) {
    // fix trailing
    for _ in 0..max_recursive {
        let mut bad_nodes = Vec::new();
        let mut removed = false;
        for (node, val) in g.degrees() {
            if val <= 1 {
                // trailing node
                if g.node(node).size < min_support {
                    bad_nodes.push(node);
                }
            }
        }
        for node in bad_nodes {
            g.remove_node(node);
            removed = true;
        }

        if !removed {
            break;
        }
    }
}

/// `clean_network.py::mod_bfs_edges`
///
/// **Origin: networkx's `generic_bfs_edges`** (BSD 3-Clause), copied into Panaroo and
/// modified to also yield the remaining depth. Panaroo's own docstring says so. See
/// `NOTICE.md`. Yields `(parent, child, depth_now)`.
///
/// This is the *older* networkx implementation, which holds live neighbour iterators in a
/// deque and advances one per outer step, rather than networkx 3.x's level-order rewrite.
/// Both produce the same order; the transcription follows Panaroo's version.
/// Differentially tested against it in `tests/parity/graph/check.sh`.
pub fn mod_bfs_edges(
    g: &Graph,
    source: usize,
    depth_limit: Option<usize>,
) -> Vec<(usize, usize, usize)> {
    let depth_limit = depth_limit.unwrap_or_else(|| g.number_of_nodes());

    let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
    visited.insert(source);

    // The Python holds live neighbour *iterators* in the queue and advances one step per
    // outer iteration, so a parent stays at the front until its neighbours are exhausted.
    // An index into the neighbour list reproduces that exactly.
    let mut queue: std::collections::VecDeque<(usize, usize, Vec<usize>, usize)> =
        std::collections::VecDeque::new();
    queue.push_back((source, depth_limit, g.neighbors(source), 0));

    let mut out = Vec::new();
    while let Some(front) = queue.front_mut() {
        if front.3 >= front.2.len() {
            queue.pop_front();
            continue;
        }
        let parent = front.0;
        let depth_now = front.1;
        let child = front.2[front.3];
        front.3 += 1;

        if visited.insert(child) {
            out.push((parent, child, depth_now));
            if depth_now > 1 {
                let nbrs = g.neighbors(child);
                queue.push_back((child, depth_now - 1, nbrs, 0));
            }
        }
    }
    out
}

/// `clean_network.py::single_linkage`
///
/// Slices `distances_bwtn_centroids` down to the centroids of `neighbours`, runs
/// `connected_components`, merges labels that share a neighbour node, then groups.
///
/// Two parity notes:
///  - the label numbering from [`crate::support::sparse::connected_components`] must match
///    scipy's exactly (PORTING_PLAN.md §6.5), since `np.unique(labels)` fixes cluster order;
///  - `l = list(set(labels[neigh_array == neigh]))` iterates a Python `set` of numpy ints,
///    so `l[0]` — the label everything else is rewritten to — depends on set order.
///
/// PORTING_PLAN.md §9 item 6: the double fancy-index is the pipeline's hottest line. It
/// stays as-is until parity is signed off.
pub fn single_linkage(
    g: &Graph,
    distances_bwtn_centroids: &CsrMatrix,
    centroid_to_index: &HashMap<String, usize>,
    neighbours: &[usize],
) -> Vec<Vec<usize>> {
    let mut index: Vec<usize> = Vec::new();
    let mut neigh_offsets = Vec::with_capacity(neighbours.len() + 1);
    neigh_offsets.push(0);
    for &neigh in neighbours {
        for sid in &g.node(neigh).centroid {
            index.push(centroid_to_index[sid]);
        }
        neigh_offsets.push(index.len());
    }

    let mut labels = undirected_component_labels_for_index(distances_bwtn_centroids, &index);

    let mut label_sets: Vec<std::collections::BTreeSet<i64>> =
        vec![std::collections::BTreeSet::new(); neighbours.len()];
    for (neigh_i, window) in neigh_offsets.windows(2).enumerate() {
        for &label in &labels[window[0]..window[1]] {
            label_sets[neigh_i].insert(label);
        }
    }
    let max_label = labels.iter().copied().max().unwrap_or(-1);
    let mut label_merges = DisjointSet::new((max_label + 1) as usize);
    for labels_for_neigh in &label_sets {
        // Tier D: `sorted(set(labels[neigh_array == neigh]))`, so l[0] is the smallest
        // label rather than an arbitrary one. Every other label in the group is rewritten
        // to it, which is what merges components that share a node.
        if labels_for_neigh.len() > 1 {
            let mut labels_for_neigh = labels_for_neigh.iter().copied();
            if let Some(target) = labels_for_neigh.next() {
                for i in labels_for_neigh {
                    label_merges.union_min(target as usize, i as usize);
                }
            }
        }
    }
    for label in &mut labels {
        *label = label_merges.find(*label as usize) as i64;
    }

    // `np.unique(labels)` is sorted, so cluster order follows label value.
    let mut groups: std::collections::BTreeMap<i64, Vec<usize>> = std::collections::BTreeMap::new();
    for (neigh_i, window) in neigh_offsets.windows(2).enumerate() {
        let neigh = neighbours[neigh_i];
        for &label in &labels[window[0]..window[1]] {
            groups.entry(label).or_default().push(neigh);
        }
    }
    groups
        .into_values()
        .map(|mut v| {
            // `del_dups` from isvalid (the later import wins in clean_network.py), which
            // dedups in place preserving first-occurrence order.
            crate::isvalid::del_dups(&mut v);
            v
        })
        .collect()
}

/// Equivalent to `connected_components(distances[index][:, index], directed=False)[1]`.
///
/// `single_linkage` calls this for every candidate node. Building the temporary CSR
/// submatrix is observable only through SciPy's component label numbering, which is the
/// component order by the smallest local vertex. A union-find over the selected positions
/// preserves that numbering while avoiding the per-call sparse allocation and transpose.
fn undirected_component_labels_for_index(distances: &CsrMatrix, index: &[usize]) -> Vec<i64> {
    let n = index.len();
    let mut col_pos: HashMap<usize, Vec<usize>> = HashMap::new();
    for (p, &c) in index.iter().enumerate() {
        col_pos.entry(c).or_default().push(p);
    }

    let mut components = DisjointSet::new(n);
    for (row_pos, &row) in index.iter().enumerate() {
        for k in distances.indptr[row]..distances.indptr[row + 1] {
            if distances.data[k] == 0 {
                continue;
            }
            if let Some(positions) = col_pos.get(&distances.indices[k]) {
                for &col in positions {
                    components.union_min(row_pos, col);
                }
            }
        }
    }

    let mut root_to_label: HashMap<usize, i64> = HashMap::new();
    let mut next_label = 0i64;
    let mut labels = Vec::with_capacity(n);
    for v in 0..n {
        let root = components.find(v);
        let label = *root_to_label.entry(root).or_insert_with(|| {
            let label = next_label;
            next_label += 1;
            label
        });
        labels.push(label);
    }
    labels
}

#[derive(Debug, Clone)]
struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        while self.parent[x] != x {
            let parent = self.parent[x];
            self.parent[x] = root;
            x = parent;
        }
        root
    }

    fn union_min(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if ra < rb {
            self.parent[rb] = ra;
        } else {
            self.parent[ra] = rb;
        }
    }
}

/// `clean_network.py::collapse_families`
///
/// Returns `(G, distances_bwtn_centroids, centroid_to_index)`.
///
/// Called three times from `main` with different arguments:
///  1. `correct_mistranslations=True` — DNA-level, thresholds `[0.99, 0.98, 0.95, 0.9]`
///  2. `correct_mistranslations=False` — protein-level, `[0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]`
///  3. again after `find_missing`, reusing the matrix from (2)
///
/// `node_count` starts at `max(G.nodes()) + 10` and increments per merge; those IDs are
/// output-visible.
#[allow(clippy::too_many_arguments)]
pub fn collapse_families(
    g: &mut Graph,
    seqid_to_centroid: &HashMap<String, String>,
    outdir: &str,
    family_threshold: f64,
    dna_error_threshold: f64,
    family_len_dif_percent: Option<f64>,
    correct_mistranslations: bool,
    length_outlier_support_proportion: f64,
    n_cpu: i64,
    quiet: bool,
    distances_bwtn_centroids: Option<CsrMatrix>,
    centroid_to_index: Option<HashMap<String, usize>>,
    depths: &[usize],
    search_genome_ids: Option<&IntBitSet>,
) -> (CsrMatrix, HashMap<String, usize>) {
    use crate::cdhit::{iterative_cdhit, pwdist_edlib};
    use crate::merge_nodes::merge_node_cluster;
    use crate::support::pyfmt::PyNum;
    use std::collections::{BTreeSet, HashMap as Map, HashSet};

    let mut node_count = g.nodes().iter().copied().max().unwrap() + 10;

    let threshold: Vec<f64> = if correct_mistranslations {
        vec![0.99, 0.98, 0.95, 0.9]
    } else {
        vec![0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]
    };

    // precluster for speed
    let (distances_bwtn_centroids, centroid_to_index) = if correct_mistranslations {
        let cdhit_clusters = iterative_cdhit(
            g,
            outdir,
            true,
            PyNum::from_opt_f64(family_len_dif_percent),
            PyNum::Float(0.0),
            99999999,
            PyNum::Float(0.0),
            99999999,
            false,
            false,
            1,
            true,
            Some(7),
            &threshold,
            n_cpu,
        );
        pwdist_edlib(g, &cdhit_clusters, dna_error_threshold, true, n_cpu)
    } else if distances_bwtn_centroids.is_none() {
        let cdhit_clusters = iterative_cdhit(
            g,
            outdir,
            false,
            PyNum::from_opt_f64(family_len_dif_percent),
            PyNum::Float(0.0),
            99999999,
            PyNum::Float(0.0),
            99999999,
            true,
            false,
            1,
            true,
            None,
            &threshold,
            n_cpu,
        );
        pwdist_edlib(g, &cdhit_clusters, family_threshold, false, n_cpu)
    } else {
        (
            distances_bwtn_centroids.unwrap(),
            centroid_to_index.unwrap(),
        )
    };

    // keep track of centroids for each sequence. Need this to resolve clashes
    let mut seqid_to_index: Map<String, usize> = Map::new();
    for node in g.nodes() {
        for sid in g.node(node).seq_ids.iter() {
            if sid.contains("refound") {
                seqid_to_index.insert(
                    sid.clone(),
                    centroid_to_index[&g.node(node).long_centroid_id.1],
                );
            } else {
                seqid_to_index.insert(sid.clone(), centroid_to_index[&seqid_to_centroid[sid]]);
            }
        }
    }

    let nz = distances_bwtn_centroids.nonzero();
    let nonzero_dist: HashSet<(usize, usize)> =
        nz.0.iter().copied().zip(nz.1.iter().copied()).collect();

    // node -> genome -> set of centroid indices
    let mut node_mem_index: Map<usize, Map<usize, BTreeSet<usize>>> = Map::new();
    for n in g.nodes() {
        let mut per_genome: Map<usize, BTreeSet<usize>> = Map::new();
        for sid in g.node(n).seq_ids.iter() {
            let genome: usize = sid.split('_').next().unwrap().parse().unwrap();
            per_genome
                .entry(genome)
                .or_default()
                .insert(seqid_to_index[sid]);
        }
        node_mem_index.insert(n, per_genome);
    }

    for &depth in depths {
        if !quiet {
            println!("Processing depth:  {depth}");
        }
        let mut search_space: BTreeSet<usize> = match search_genome_ids {
            None => g.nodes().into_iter().collect(),
            Some(ids) => g
                .nodes()
                .into_iter()
                // `!a.intersection(b).is_empty()` materialised a whole IntBitSet to answer
                // a boolean. `isdisjoint` scans the same `min(len)` word range with the same
                // `&`, so it is the identical predicate -- but it short-circuits on the first
                // overlapping word and allocates nothing.
                .filter(|&n| !g.node(n).members.isdisjoint(ids))
                .collect(),
        };

        let mut iteration_num = 1;
        while !search_space.is_empty() {
            // look for nodes to merge. Tier D: sorted(search_space), which is what makes
            // the BFS seeding -- and hence the whole merge outcome -- deterministic.
            let temp_node_list: Vec<usize> = search_space.iter().copied().collect();
            let mut removed_nodes: HashSet<usize> = HashSet::new();
            if !quiet {
                println!("Iteration:  {iteration_num}");
            }
            iteration_num += 1;

            for node in temp_node_list {
                if removed_nodes.contains(&node) {
                    continue;
                }
                if g.degree(node) <= 2 {
                    search_space.remove(&node);
                    removed_nodes.insert(node);
                    continue;
                }

                // find neighbouring nodes and cluster their centroid with cdhit
                let mut neighbours: Vec<usize> =
                    crate::support::graph::bfs_edges(g, node, Some(depth))
                        .into_iter()
                        .map(|(_, v)| v)
                        .collect();
                neighbours.push(node);

                // find clusters
                let clusters = single_linkage(
                    g,
                    &distances_bwtn_centroids,
                    &centroid_to_index,
                    &neighbours,
                );

                for cluster in clusters {
                    // check if there are any to collapse
                    if cluster.len() <= 1 {
                        continue;
                    }

                    // check for conflicts
                    let mut seen = g.node(cluster[0]).members.copy();
                    let mut noconflict = true;
                    for &n in &cluster[1..] {
                        if !seen.isdisjoint(&g.node(n).members) {
                            noconflict = false;
                            break;
                        }
                        seen.union_update(&g.node(n).members);
                    }

                    if noconflict {
                        // no conflicts so merge
                        node_count += 1;
                        for &neig in &cluster {
                            removed_nodes.insert(neig);
                            search_space.remove(&neig);
                        }
                        merge_node_cluster(g, &cluster, node_count, !correct_mistranslations, true);
                        merge_mem_index(&mut node_mem_index, &cluster, node_count);
                        search_space.insert(node_count);
                    } else {
                        // merge if the centroids don't conflict and the nodes are adjacent
                        // in the conflicting genome -- a mistranslation / frame shift /
                        // premature stop where one gene was split in some genomes.
                        let mut cluster = cluster;
                        // sort by size, descending. `sorted(reverse=True)` is stable, so
                        // equal sizes keep their current order.
                        // `sorted(reverse=True)` is stable, so equal sizes keep their
                        // current order -- `sort_by_key` is stable too, so this is exact.
                        cluster.sort_by_key(|&a| std::cmp::Reverse(g.node(a).size));

                        let mut node_mem_count: Map<usize, usize> = Map::new();
                        for &n in &cluster {
                            for m in g.node(n).members.iter() {
                                *node_mem_count.entry(m).or_insert(0) += 1;
                            }
                        }
                        let mem_count: Vec<usize> = node_mem_count.values().copied().collect();
                        let ones = mem_count.iter().filter(|&&c| c == 1).count();
                        let merge_same_members = (ones as f64 / mem_count.len() as f64)
                            >= length_outlier_support_proportion;

                        while !cluster.is_empty() {
                            let mut sub_clust = vec![cluster[0]];
                            let n_a = cluster[0];
                            for &n_b in &cluster[1..] {
                                let mem_inter: Vec<usize> = g
                                    .node(n_a)
                                    .members
                                    .intersection(&g.node(n_b).members)
                                    .iter()
                                    .collect();
                                if !mem_inter.is_empty() {
                                    if merge_same_members {
                                        let mut shouldmerge = true;
                                        let ca: HashSet<&String> =
                                            g.node(n_a).centroid.iter().collect();
                                        if g.node(n_b).centroid.iter().any(|c| ca.contains(c)) {
                                            shouldmerge = false;
                                        }

                                        if shouldmerge {
                                            let mut edge_mem_count: Map<usize, usize> = Map::new();
                                            'outer: for members in
                                                crate::merge_nodes::gen_edge_iterables(
                                                    g,
                                                    &g.edges_of(&[n_a, n_b]),
                                                )
                                            {
                                                for e in members.iter() {
                                                    let c = edge_mem_count.entry(e).or_insert(0);
                                                    *c += 1;
                                                    if *c > 3 {
                                                        shouldmerge = false;
                                                        break 'outer;
                                                    }
                                                }
                                            }
                                        }

                                        if shouldmerge {
                                            'chk: for imem in &mem_inter {
                                                let a = node_mem_index[&n_a].get(imem);
                                                let b = node_mem_index[&n_b].get(imem);
                                                if let (Some(a), Some(b)) = (a, b) {
                                                    for &sid_a in a {
                                                        for &sid_b in b {
                                                            if nonzero_dist
                                                                .contains(&(sid_a, sid_b))
                                                                || nonzero_dist
                                                                    .contains(&(sid_b, sid_a))
                                                            {
                                                                shouldmerge = false;
                                                                break 'chk;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        if shouldmerge {
                                            sub_clust.push(n_b);
                                        }
                                    }
                                } else {
                                    sub_clust.push(n_b);
                                }
                            }

                            if sub_clust.len() > 1 {
                                let clique_clusters = single_linkage(
                                    g,
                                    &distances_bwtn_centroids,
                                    &centroid_to_index,
                                    &sub_clust,
                                );
                                for clust in clique_clusters {
                                    if clust.len() <= 1 {
                                        continue;
                                    }
                                    node_count += 1;
                                    for &neig in &clust {
                                        removed_nodes.insert(neig);
                                        search_space.remove(&neig);
                                    }
                                    merge_node_cluster(
                                        g,
                                        &clust,
                                        node_count,
                                        !correct_mistranslations,
                                        false,
                                    );
                                    merge_mem_index(&mut node_mem_index, &clust, node_count);
                                    search_space.insert(node_count);
                                }
                            }

                            let sub: HashSet<usize> = sub_clust.iter().copied().collect();
                            cluster.retain(|n| !sub.contains(n));
                        }
                    }
                }

                search_space.remove(&node);
            }
        }
    }

    (distances_bwtn_centroids, centroid_to_index)
}

/// Fold the merged nodes' `node_mem_index` entries into the new node's.
///
/// Not a Python function — the same eight lines appear twice in `collapse_families`, and
/// they mutate a structure whose contents are load-bearing. One copy so the two merge
/// paths cannot drift.
fn merge_mem_index(
    node_mem_index: &mut HashMap<usize, HashMap<usize, std::collections::BTreeSet<usize>>>,
    cluster: &[usize],
    node_count: usize,
) {
    let mut acc = node_mem_index.get(&cluster[0]).cloned().unwrap_or_default();
    for n in &cluster[1..] {
        if let Some(other) = node_mem_index.get(n) {
            for (m, v) in other {
                acc.entry(*m).or_default().extend(v.iter().copied());
            }
        }
        node_mem_index.remove(n);
    }
    node_mem_index.insert(node_count, acc);
}

/// `clean_network.py::collapse_paralogs`
///
/// PORTING_PLAN.md §9 item 7: the `shortest_path_length` call is uncapped and runs a full
/// BFS per (paralog, reference) pair; `ref_context` is recomputed inside the inner loop.
/// Both stay until parity is signed off.
///
/// `max(member_paralogs.items(), key=lambda x: len(x[1]))` returns the **first** maximum —
/// fold with a strict `>`, do not use `max_by_key` (PORTING_PLAN.md §6.4).
pub fn collapse_paralogs(
    g: &mut Graph,
    centroid_contexts: &mut crate::support::pydict::PyDict<String, Vec<(usize, usize)>>,
    max_context: usize,
    _quiet: bool,
) {
    use crate::merge_nodes::merge_node_cluster;
    use crate::support::pydict::PyDict;
    use std::collections::BTreeSet;

    let mut node_count = g.nodes().iter().copied().max().unwrap() + 10;

    // first sort by context length, context dist to ensure ties are broken the same way.
    // Python sorts a list of [node, genome] lists lexicographically.
    for v in centroid_contexts.values_mut() {
        v.sort();
    }

    // set up for context search
    let mut centroid_to_index: HashMap<String, usize> = HashMap::new();
    let mut ncentroids: i64 = -1;
    for node in g.nodes() {
        let centroid = g.node(node).centroid[0].clone();
        if !centroid_to_index.contains_key(&centroid) {
            ncentroids += 1;
            centroid_to_index.insert(centroid, ncentroids as usize);
        }
    }
    let ncentroids = (ncentroids + 1) as usize;

    // `centroid_contexts` is a defaultdict built by generate_network in gene order, and
    // upstream iterates it in that insertion order. A dict, not a set -- so this is
    // deterministic upstream and is NOT a Tier D site: the insertion order must be
    // preserved, not replaced with a sort. Hence PyDict.
    let centroids: Vec<String> = centroid_contexts.keys().cloned().collect();

    for centroid in centroids {
        // calculate distance
        let mut member_paralogs: PyDict<usize, Vec<(usize, usize)>> = PyDict::new();
        for para in centroid_contexts.get(&centroid).unwrap() {
            if !member_paralogs.contains_key(&para.1) {
                member_paralogs.insert(para.1, Vec::new());
            }
            member_paralogs.get_mut(&para.1).unwrap().push(*para);
        }

        // `max(..., key=len)` returns the FIRST maximum -- fold with a strict `>`.
        let mut ref_paralogs: Option<&Vec<(usize, usize)>> = None;
        for v in member_paralogs.values() {
            if ref_paralogs.is_none_or(|best| v.len() > best.len()) {
                ref_paralogs = Some(v);
            }
        }
        let ref_paralogs = ref_paralogs.expect("member_paralogs is not empty");

        // for each paralog find its closest reference paralog
        let mut cluster_dict: PyDict<usize, BTreeSet<usize>> = PyDict::new();
        let mut cluster_mems: PyDict<usize, BTreeSet<usize>> = PyDict::new();
        for (c, r) in ref_paralogs.iter().enumerate() {
            cluster_dict.insert(c, BTreeSet::from([r.0]));
            cluster_mems.insert(c, BTreeSet::from([r.1]));
        }

        for para in centroid_contexts.get(&centroid).unwrap() {
            let mut d_max = usize::MAX;
            let mut best_cluster: Option<usize> = None;

            if para.1 == ref_paralogs[0].1 {
                // this is the reference so skip
                continue;
            }

            // first attempt by shortest path
            for (c, r) in ref_paralogs.iter().enumerate() {
                if cluster_mems.get(&c).unwrap().contains(&para.1) {
                    // dont match paralogs of the same isolate
                    continue;
                }
                // PORTING_PLAN.md §9 item 7: uncapped full-graph BFS per pair. Stays.
                if let Some(d) = crate::support::graph::shortest_path_length(g, r.0, para.0) {
                    if d < d_max {
                        d_max = d;
                        best_cluster = Some(c);
                    }
                }
            }

            // if this fails use context
            if d_max == usize::MAX {
                best_cluster = Some(0);
                let mut s_max = f64::NEG_INFINITY;
                let mut para_context = vec![0f64; ncentroids];
                for (_, node, depth) in mod_bfs_edges(g, para.0, Some(max_context)) {
                    para_context[centroid_to_index[&g.node(node).centroid[0]]] = depth as f64;
                }
                for (c, r) in ref_paralogs.iter().enumerate() {
                    if cluster_mems.get(&c).unwrap().contains(&para.1) {
                        continue;
                    }
                    let mut ref_context = vec![0f64; ncentroids];
                    for (_, node, depth) in mod_bfs_edges(g, r.0, Some(max_context)) {
                        ref_context[centroid_to_index[&g.node(node).centroid[0]]] = depth as f64;
                    }
                    // sum(1 / (1 + abs((para - ref)[para*ref != 0])))
                    let terms: Vec<f64> = (0..ncentroids)
                        .filter(|&i| para_context[i] * ref_context[i] != 0.0)
                        .map(|i| 1.0 / (1.0 + (para_context[i] - ref_context[i]).abs()))
                        .collect();
                    let sc = crate::support::npmath::np_sum_f64(&terms);
                    if sc > s_max {
                        s_max = sc;
                        best_cluster = Some(c);
                    }
                }
            }

            let bc = best_cluster.expect("best_cluster set");
            cluster_dict.get_mut(&bc).unwrap().insert(para.0);
            cluster_mems.get_mut(&bc).unwrap().insert(para.1);
        }

        // merge
        for (_cluster, members) in cluster_dict.items() {
            if members.len() < 2 {
                continue;
            }
            node_count += 1;
            // Tier D: sorted(cluster_dict[cluster]) -- a BTreeSet is already sorted.
            let nodes: Vec<usize> = members.iter().copied().collect();
            merge_node_cluster(g, &nodes, node_count, true, true);
        }
    }
}

/// `clean_network.py::merge_paralogs`
///
/// PORTING_PLAN.md §9 item 8: hand-rolled connected components by repeated list copy.
/// Preserve the algorithm — the order in which `merge_clusters` is built determines the
/// `node_count` assignment.
pub fn merge_paralogs(g: &mut Graph) {
    use crate::merge_nodes::merge_node_cluster;
    use crate::support::pydict::PyDict;
    use std::collections::BTreeSet;

    let mut node_count = g.nodes().iter().copied().max().unwrap() + 10;

    // group paralog nodes by centroid
    let mut paralog_centroids: PyDict<String, Vec<usize>> = PyDict::new();
    for node in g.nodes() {
        if g.node(node).paralog {
            for centroid in &g.node(node).centroid {
                if !paralog_centroids.contains_key(&centroid) {
                    paralog_centroids.insert(centroid.clone(), Vec::new());
                }
                paralog_centroids.get_mut(&centroid).unwrap().push(node);
            }
        }
    }

    // find nodes that share common centroids
    //
    // PORTING_PLAN.md §9 item 8: the Python does this with a `while` loop that unpacks
    // `first, *rest`, copying the remaining list every round -- O(n^2) or worse. Preserved,
    // because the order in which merge_clusters is built decides node_count assignment.
    let mut rest: Vec<Vec<usize>> = paralog_centroids.values().cloned().collect();
    let mut merge_clusters: Vec<BTreeSet<usize>> = Vec::new();
    while !rest.is_empty() {
        let mut first: BTreeSet<usize> = rest[0].iter().copied().collect();
        let mut remaining: Vec<Vec<usize>> = rest[1..].to_vec();
        let mut lf = usize::MAX;
        while first.len() > lf || lf == usize::MAX {
            lf = first.len();
            let mut rest2: Vec<Vec<usize>> = Vec::new();
            for r in remaining {
                if r.iter().any(|n| first.contains(n)) {
                    first.extend(r.iter().copied());
                } else {
                    rest2.push(r);
                }
            }
            remaining = rest2;
            if first.len() == lf {
                break;
            }
        }
        merge_clusters.push(first);
        rest = remaining;
    }

    // merge paralog nodes that share the same centroid
    for temp_c in merge_clusters {
        if temp_c.len() > 1 {
            node_count += 1;
            // Tier D: sorted(temp_c) -- a BTreeSet is already sorted.
            let nodes: Vec<usize> = temp_c.into_iter().collect();
            merge_node_cluster(g, &nodes, node_count, true, false);
        }
    }
}

/// `clean_network.py::clean_misassembly_edges`
///
/// `bad_edges` is a Python `set` of `(node, neighbour)` tuples, but it is only iterated to
/// call `remove_edge`, and removal is idempotent under the `has_edge` guard — so set order
/// does not reach the output here. Verify this claim when translating.
pub fn clean_misassembly_edges(g: &mut Graph, edge_support_threshold: f64) {
    // `bad_edges` is a Python set of (node, neighbour) tuples, but it is only iterated to
    // call remove_edge under a has_edge guard, so its order is not observable. Audited;
    // deliberately not in the Tier D patch.
    let mut bad_edges: Vec<(usize, usize)> = Vec::new();
    let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();

    // remove edges with low support near contig ends
    for node in g.nodes() {
        for neigh in g.neighbors(node) {
            if g.node(neigh).has_end && (g.edge(node, neigh).size as f64) < edge_support_threshold {
                if seen.insert((node, neigh)) {
                    bad_edges.push((node, neigh));
                }
            }
        }
    }

    // remove edges that have much lower support than the nodes they connect
    for edge in g.edges() {
        let esize = g.edge(edge.0, edge.1).size as f64;
        let nmin = g.node(edge.0).size.min(g.node(edge.1).size) as f64;
        if esize < (0.05 * nmin) && esize < edge_support_threshold {
            if seen.insert(edge) {
                bad_edges.push(edge);
            }
        }
    }

    for edge in bad_edges {
        if g.has_edge(edge.0, edge.1) {
            g.remove_edge(edge.0, edge.1);
        }
    }
}

/// `clean_network.py::identify_possible_highly_variable` — not reachable from `main`.
/// Phase 6.
pub fn identify_possible_highly_variable(
    _g: &mut Graph,
    _cycle_threshold_max: usize,
    _cycle_threshold_min: usize,
    _size_diff_threshold: f64,
) {
    panic!("noimpl: clean_network::identify_possible_highly_variable")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn coo(data: Vec<i64>, r: Vec<usize>, c: Vec<usize>, n: usize) -> CsrMatrix {
        CsrMatrix::from_coo(data, r, c, (n, n))
    }

    #[test]
    fn direct_component_labels_match_sparse_slice() {
        let m = coo(vec![1, 1, 1, 1], vec![0, 1, 4, 3], vec![3, 2, 0, 1], 5);
        for index in [vec![3, 1, 2], vec![4, 0, 3, 1, 2], vec![3, 1, 2, 1]] {
            let sub = m.submatrix(&index);
            let (_, expected) = crate::support::sparse::connected_components(&sub, false);
            assert_eq!(undirected_component_labels_for_index(&m, &index), expected);
        }
    }
}
