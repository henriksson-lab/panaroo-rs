//! Translation of `panaroo/panaroo/generate_network.py`.

use crate::support::graph::Graph;
use std::collections::HashMap;

/// What `generate_network` returns: `(G, centroid_context, seqid_to_centroid)`.
pub struct NetworkResult {
    pub graph: Graph,
    /// `centroid_context` — centroid ID -> list of `[node, genome_id]`.
    ///
    /// A `PyDict`, not a `HashMap`: `collapse_paralogs` iterates it, and upstream's
    /// `defaultdict` is insertion-ordered. That order is deterministic in Python (a dict,
    /// not a set — so *not* a Tier D site), which means it has to be preserved here rather
    /// than replaced with a sort.
    pub centroid_context: crate::support::pydict::PyDict<String, Vec<(usize, usize)>>,
    /// `seqid_to_centroid` — clustering ID -> centroid clustering ID.
    pub seqid_to_centroid: HashMap<String, String>,
}

/// `generate_network.py::generate_network`
///
/// Walks the protein FASTA in file order, which is the order genes appear on each contig,
/// and threads a node per cd-hit cluster onto the graph, adding an edge from the previous
/// gene unless we are at the start of a contig (`loc[-1] == "0"`).
///
/// Details that carry into every output file:
///  - node IDs are the cd-hit cluster numbers; paralogs get fresh IDs counting up from
///    `len(cluster_members)`
///  - `prev` is *not* reset between genomes; the `loc[-1] == "0"` test on the first gene of
///    a contig is what closes the previous contig by setting `hasEnd`
///  - `paralogs` is a Python `set` of cluster ints, but is only membership-tested
///
/// PORTING_PLAN.md §9 item 11: the Python parses the protein FASTA through Biopython purely
/// to read record IDs. Left as-is.
pub fn generate_network(
    cluster_file: &str,
    data_file: &str,
    prot_seq_file: &str,
    all_dna: bool,
) -> NetworkResult {
    use crate::support::graph::EdgeAttrs;
    use crate::support::intbitset::IntBitSet;
    use crate::support::pydict::PyDict;
    use std::collections::HashSet;

    // associate sequences with their clusters
    let mut seq_to_cluster: HashMap<String, usize> = HashMap::new();
    let mut seqid_to_centroid: HashMap<String, String> = HashMap::new();
    let mut cluster_centroids: PyDict<usize, String> = PyDict::new();
    let mut cluster_members: PyDict<usize, Vec<String>> = PyDict::new();
    {
        let text = std::fs::read_to_string(cluster_file)
            .unwrap_or_else(|e| panic!("could not read {cluster_file}: {e}"));
        let mut cluster = 0usize;
        for line in text.lines() {
            if line.starts_with('>') {
                cluster = line.split_whitespace().last().unwrap().parse().unwrap();
            } else {
                let seq = line
                    .split('>')
                    .nth(1)
                    .unwrap()
                    .split("...")
                    .next()
                    .unwrap()
                    .to_string();
                seq_to_cluster.insert(seq.clone(), cluster);
                let genome = seq.split('_').next().unwrap().to_string();
                if !cluster_members.contains_key(&cluster) {
                    cluster_members.insert(cluster, Vec::new());
                }
                cluster_members.get_mut(&cluster).unwrap().push(genome);
                if line.split_whitespace().last() == Some("*") {
                    cluster_centroids.insert(cluster, seq);
                }
            }
        }
    }

    // determine paralogs if required
    let mut paralogs: HashSet<usize> = HashSet::new();
    for clust in cluster_members.keys().copied().collect::<Vec<_>>() {
        // `cluster_members[clust]` already holds the genome field, so the extra
        // `s.split("_")[0]` upstream is a no-op. Transcribed as the identity it is.
        let genomes = cluster_members.get(&clust).unwrap();
        let uniq: HashSet<&String> = genomes.iter().collect();
        if genomes.len() != uniq.len() {
            paralogs.insert(clust);
        }
    }

    // Load meta data such as sequence and annotation
    let mut cluster_centroid_data: HashMap<usize, CentroidData> = HashMap::new();
    let centroid_ids: HashSet<String> = cluster_centroids.values().cloned().collect();
    {
        let text = std::fs::read_to_string(data_file)
            .unwrap_or_else(|e| panic!("could not read {data_file}: {e}"));
        for line in text.lines().skip(1) {
            let f: Vec<&str> = line.trim_end().split(',').collect();
            if f.len() < 8 {
                continue;
            }
            if centroid_ids.contains(f[2]) {
                // this is a cluster centroid so keep it
                cluster_centroid_data.insert(
                    seq_to_cluster[f[2]],
                    CentroidData {
                        prot_sequence: f[4].to_string(),
                        dna_sequence: f[5].to_string(),
                        annotation: f[6].to_string(),
                        description: f[7].to_string(),
                    },
                );
            }
        }
    }

    // load headers which contain adjacency information
    //
    // PORTING_PLAN.md §9 item 11: the Python parses the whole protein FASTA through
    // Biopython just to read record ids. Left as-is.
    let seq_ids: Vec<String> = crate::support::seqio::parse_fasta_file(prot_seq_file)
        .into_iter()
        .map(|r| r.id)
        .collect();

    // build graph using adjacency information and optionally split paralogs
    let mut g = Graph::new();
    let mut centroid_context: PyDict<String, Vec<(usize, usize)>> = PyDict::new();
    let mut n_nodes = cluster_members.len();
    let mut prev: Option<usize> = None;

    for id in &seq_ids {
        let current_cluster = seq_to_cluster[id];
        seqid_to_centroid.insert(
            id.clone(),
            cluster_centroids.get(&current_cluster).unwrap().clone(),
        );
        let loc: Vec<&str> = id.split('_').collect();
        let genome_id: usize = loc[0].parse().unwrap();
        let centroid = cluster_centroids.get(&current_cluster).unwrap().clone();
        let data = cluster_centroid_data
            .get(&current_cluster)
            .unwrap_or_else(|| panic!("KeyError: cluster {current_cluster}"))
            .clone();

        if loc[loc.len() - 1] == "0" {
            // we're at the start of a contig
            if let Some(p) = prev {
                g.node_mut(p).has_end = true;
            }
            let mut cur = current_cluster;
            let cur_is_paralog = paralogs.contains(&current_cluster);
            if g.has_node(cur) && !cur_is_paralog {
                let n = g.node_mut(cur);
                n.size += 1;
                n.members.add(genome_id);
                n.seq_ids.insert(id.clone());
                n.has_end = true;
                n.lengths.push(data.dna_sequence.len());
                if all_dna {
                    n.dna.push(data.dna_sequence.clone());
                }
            } else {
                if cur_is_paralog {
                    // create a new paralog
                    n_nodes += 1;
                    cur = n_nodes;
                    if !centroid_context.contains_key(&centroid) {
                        centroid_context.insert(centroid.clone(), Vec::new());
                    }
                    centroid_context
                        .get_mut(&centroid)
                        .unwrap()
                        .push((cur, genome_id));
                }
                // add non paralog node
                g.add_node(
                    cur,
                    node_attrs(genome_id, id, true, &data, &centroid, cur_is_paralog),
                );
            }
            prev = Some(cur);
        } else {
            let is_paralog = paralogs.contains(&current_cluster);
            if is_paralog {
                // create a new paralog
                n_nodes += 1;
                let neighbour = n_nodes;
                if !centroid_context.contains_key(&centroid) {
                    centroid_context.insert(centroid.clone(), Vec::new());
                }
                centroid_context
                    .get_mut(&centroid)
                    .unwrap()
                    .push((neighbour, genome_id));
                g.add_node(
                    neighbour,
                    node_attrs(genome_id, id, false, &data, &centroid, true),
                );
                // add edge between nodes
                g.add_edge(
                    prev.expect("prev set"),
                    neighbour,
                    EdgeAttrs {
                        size: 1,
                        members: IntBitSet::from_iter_ints([genome_id]),
                        genome_ids: None,
                    },
                );
                prev = Some(neighbour);
            } else {
                if !g.has_node(current_cluster) {
                    // we need to add the gene in
                    g.add_node(
                        current_cluster,
                        node_attrs(genome_id, id, false, &data, &centroid, false),
                    );
                    // add edge between nodes
                    g.add_edge(
                        prev.expect("prev set"),
                        current_cluster,
                        EdgeAttrs {
                            size: 1,
                            members: IntBitSet::from_iter_ints([genome_id]),
                            genome_ids: None,
                        },
                    );
                } else {
                    {
                        let n = g.node_mut(current_cluster);
                        n.size += 1;
                        n.members.add(genome_id);
                        n.seq_ids.insert(id.clone());
                        n.lengths.push(data.dna_sequence.len());
                        if all_dna {
                            n.dna.push(data.dna_sequence.clone());
                        }
                    }
                    let p = prev.expect("prev set");
                    if g.has_edge(p, current_cluster) {
                        let e = g.edge_mut(p, current_cluster);
                        e.size += 1;
                        e.members.add(genome_id);
                    } else {
                        g.add_edge(
                            p,
                            current_cluster,
                            EdgeAttrs {
                                size: 1,
                                members: IntBitSet::from_iter_ints([genome_id]),
                                genome_ids: None,
                            },
                        );
                    }
                }
                prev = Some(current_cluster);
            }
        }
    }

    if let Some(p) = prev {
        g.node_mut(p).has_end = true;
    }

    NetworkResult {
        graph: g,
        centroid_context,
        seqid_to_centroid,
    }
}

/// The four `gene_data.csv` columns kept for a cluster centroid.
#[derive(Debug, Clone)]
struct CentroidData {
    prot_sequence: String,
    dna_sequence: String,
    annotation: String,
    description: String,
}

/// Build the `NodeAttrs` for a fresh node.
///
/// Not a Python function: the Python repeats the same 20-line `G.add_node(...)` kwarg block
/// three times. Rule 2 keeps duplicated *logic* duplicated, but this is a struct literal,
/// and three copies of it would be three places for a field to drift. The three call sites
/// differ only in `hasEnd` and `paralog`, which are parameters here.
#[allow(clippy::too_many_arguments)]
fn node_attrs(
    genome_id: usize,
    id: &str,
    has_end: bool,
    data: &CentroidData,
    centroid: &str,
    is_paralog: bool,
) -> crate::support::graph::NodeAttrs {
    use crate::support::intbitset::IntBitSet;
    use std::collections::BTreeSet;
    crate::support::graph::NodeAttrs {
        size: 1,
        centroid: vec![centroid.to_string()],
        max_len_id: 0,
        members: IntBitSet::from_iter_ints([genome_id]),
        seq_ids: BTreeSet::from([id.to_string()]),
        has_end,
        protein: vec![data.prot_sequence.clone()],
        dna: vec![data.dna_sequence.clone()],
        annotation: data.annotation.clone(),
        description: data.description.clone(),
        lengths: vec![data.dna_sequence.len()],
        long_centroid_id: (data.dna_sequence.len(), centroid.to_string()),
        paralog: is_paralog,
        merged_dna: false,
        prev_centroids: None,
        name: None,
        genome_ids: None,
        gene_ids: None,
        degrees: None,
        high_var: None,
        gml_late_attrs_before_name: false,
    }
}
