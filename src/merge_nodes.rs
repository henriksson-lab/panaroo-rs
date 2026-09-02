//! Translation of `panaroo/panaroo/merge_nodes.py`.

use crate::support::graph::{EdgeAttrs, Graph};
use crate::support::intbitset::IntBitSet;
use crate::support::pydict::PyDict;

/// Which node attribute `gen_node_iterables` is being asked for.
///
/// The Python indexes `G.nodes[n][feature]` with a runtime string. Rust needs the set of
/// legal keys spelled out; this enum is that set, not an extra function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeFeature {
    /// `'centroid'`
    Centroid,
    /// `'dna'`
    Dna,
    /// `'protein'`
    Protein,
    /// `'seqIDs'`
    SeqIds,
    /// `'annotation'` — always read with `split=";"`
    Annotation,
    /// `'description'` — always read with `split=";"`
    Description,
    /// `'prevCentroids'` — always read with `split=";"`
    PrevCentroids,
    /// `'members'`
    Members,
    /// `'lengths'`
    Lengths,
    /// `'hasEnd'`
    HasEnd,
    /// `'paralog'`
    Paralog,
    /// `'mergedDNA'`
    MergedDna,
    /// `'longCentroidID'`
    LongCentroidId,
}

/// One yielded value of `gen_node_iterables`.
///
/// Stands in for Python's dynamic typing: the generator yields whatever the attribute
/// happens to hold. Callers match on the variant they know they asked for.
#[derive(Debug, Clone)]
pub enum FeatureValue {
    StrList(Vec<String>),
    Members(IntBitSet),
    Lengths(Vec<usize>),
    Bool(bool),
    LongCentroidId((usize, String)),
}

/// `merge_nodes.py::gen_node_iterables`
///
/// ```python
/// def gen_node_iterables(G, nodes, feature, split=None):
///     for n in nodes:
///         if split is None:
///             yield G.nodes[n][feature]
///         else:
///             yield G.nodes[n][feature].split(split)
/// ```
///
/// `split` is only ever `";"`, and only for the three `str`-valued attributes.
///
/// Note Python's `"".split(";")` is `[""]`, not `[]` — an empty annotation contributes one
/// empty string, which `iter_del_dups` then dedups. Rust's `"".split(';')` agrees, so the
/// transcription is direct.
pub fn gen_node_iterables(
    g: &Graph,
    nodes: &[usize],
    feature: NodeFeature,
    split: Option<&str>,
) -> Vec<FeatureValue> {
    nodes
        .iter()
        .map(|&n| {
            let node = g.node(n);
            match feature {
                NodeFeature::Centroid => FeatureValue::StrList(node.centroid.clone()),
                NodeFeature::Dna => FeatureValue::StrList(node.dna.clone()),
                NodeFeature::Protein => FeatureValue::StrList(node.protein.clone()),
                NodeFeature::SeqIds => {
                    FeatureValue::StrList(node.seq_ids.iter().cloned().collect())
                }
                NodeFeature::Annotation => {
                    FeatureValue::StrList(split_or_whole(&node.annotation, split))
                }
                NodeFeature::Description => {
                    FeatureValue::StrList(split_or_whole(&node.description, split))
                }
                NodeFeature::PrevCentroids => FeatureValue::StrList(split_or_whole(
                    node.prev_centroids.as_deref().unwrap_or(""),
                    split,
                )),
                NodeFeature::Members => FeatureValue::Members(node.members.clone()),
                NodeFeature::Lengths => FeatureValue::Lengths(node.lengths.clone()),
                NodeFeature::HasEnd => FeatureValue::Bool(node.has_end),
                NodeFeature::Paralog => FeatureValue::Bool(node.paralog),
                NodeFeature::MergedDna => FeatureValue::Bool(node.merged_dna),
                NodeFeature::LongCentroidId => {
                    FeatureValue::LongCentroidId(node.long_centroid_id.clone())
                }
            }
        })
        .collect()
}

/// `value.split(split)` when `split` is given, else the value as a one-element list.
///
/// The `split is None` arm of `gen_node_iterables` yields the attribute itself. For a
/// `str`-valued attribute that is a bare string, but every call site that omits `split`
/// reads a list-valued attribute, so this arm is only reached through the `split=";"` path
/// in practice.
fn split_or_whole(value: &str, split: Option<&str>) -> Vec<String> {
    match split {
        Some(sep) => value.split(sep).map(|s| s.to_string()).collect(),
        None => vec![value.to_string()],
    }
}

/// `merge_nodes.py::gen_edge_iterables`
///
/// ```python
/// def gen_edge_iterables(G, edges, feature):
///     for e in edges:
///         yield G[e[0]][e[1]][feature]
/// ```
///
/// Only ever called with `feature='members'`, from `collapse_families`.
pub fn gen_edge_iterables(g: &Graph, edges: &[(usize, usize)]) -> Vec<IntBitSet> {
    edges
        .iter()
        .map(|&(u, v)| g.edge(u, v).members.clone())
        .collect()
}

/// `merge_nodes.py::iter_del_dups`
///
/// ```python
/// def iter_del_dups(iterable):
///     seen = {}
///     for f in itertools.chain.from_iterable(iterable):
///         seen[f] = None
///     return (list(seen.keys()))
/// ```
///
/// Flatten, then dedup preserving first-occurrence order. The Python uses a `dict` as an
/// ordered set, which is why the order is well defined — unlike a `set`, this one is not a
/// parity hazard.
pub fn iter_del_dups(iterable: Vec<Vec<String>>) -> Vec<String> {
    let mut seen: PyDict<String, ()> = PyDict::new();
    for inner in iterable {
        for f in inner {
            seen.insert(f, ());
        }
    }
    seen.keys().cloned().collect()
}

/// `merge_nodes.py::del_dups`
///
/// ```python
/// def del_dups(iterable):
///     seen = {}
///     for f in iterable:
///         seen[f] = None
///     return (list(seen.keys()))
/// ```
///
/// A different implementation from [`crate::isvalid::del_dups`], which compacts in place
/// and mutates its argument. Same observable result; both are kept per PORTING_PLAN.md §5.
///
/// Note `merge_nodes.py` shadows the `del_dups` it imports from `isvalid` with this
/// definition, so every unqualified `del_dups` call inside `merge_nodes` reaches *this*
/// one. `clean_network.py` imports both (`from panaroo.merge_nodes import *` then
/// `from panaroo.isvalid import del_dups`), and the later import wins — so
/// `clean_network.py:83` calls the **isvalid** version.
pub fn del_dups(iterable: Vec<String>) -> Vec<String> {
    let mut seen: PyDict<String, ()> = PyDict::new();
    for f in iterable {
        seen.insert(f, ());
    }
    seen.keys().cloned().collect()
}

/// `merge_nodes.py::merge_node_cluster`
///
/// Note two tie-breaking details that must survive translation:
///  - `nodes = sorted(nodes, key=lambda x: G.nodes[x]['size'])` — ascending, stable.
///  - the `maxLenId` loop uses `if len(s) >= max_l`, i.e. it deliberately keeps the
///    **last** longest sequence. Do not "fix" this to a `max_by_key`.
pub fn merge_node_cluster(
    g: &mut Graph,
    nodes: &[usize],
    new_node: usize,
    multi_centroid: bool,
    check_merge_mems: bool,
) {
    use crate::support::graph::{EdgeAttrs, NodeAttrs};
    use std::collections::{BTreeSet, HashSet};

    if check_merge_mems {
        let mut counts: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        for &n in nodes {
            for m in g.node(n).members.iter() {
                *counts.entry(m).or_insert(0) += 1;
            }
        }
        if counts.values().copied().max().unwrap_or(0) > 1 {
            panic!("ValueError: merging nodes with the same genome IDs!");
        }
    }

    // take node with most support as the 'consensus'.
    // `sorted(..., key=size)` is *stable*, so ties keep the caller's order -- which is why
    // the Tier D patch sorts the sets that reach here.
    let mut nodes: Vec<usize> = nodes.to_vec();
    nodes.sort_by_key(|&x| g.node(x).size);

    // First create a new node and combine the attributes
    let dna = iter_del_dups(nodes.iter().map(|&n| g.node(n).dna.clone()).collect());
    let mut max_len_id = 0usize;
    let mut max_l = 0usize;
    for (i, sq) in dna.iter().enumerate() {
        // `>=`, not `>`: the Python deliberately keeps the LAST longest sequence.
        if sq.len() >= max_l {
            max_l = sq.len();
            max_len_id = i;
        }
    }

    let mut members = g.node(nodes[0]).members.copy();
    for &n in &nodes[1..] {
        members.union_update(&g.node(n).members);
    }

    let merged_dna = if multi_centroid {
        nodes.iter().any(|&n| g.node(n).merged_dna)
    } else {
        true
    };

    let centroid = iter_del_dups(nodes.iter().map(|&n| g.node(n).centroid.clone()).collect());
    let seq_ids: BTreeSet<String> = iter_del_dups(
        nodes
            .iter()
            .map(|&n| g.node(n).seq_ids.iter().cloned().collect())
            .collect(),
    )
    .into_iter()
    .collect();
    let protein = iter_del_dups(nodes.iter().map(|&n| g.node(n).protein.clone()).collect());
    let annotation = iter_del_dups(
        nodes
            .iter()
            .map(|&n| split_semi(&g.node(n).annotation))
            .collect(),
    )
    .join(";");
    let description = iter_del_dups(
        nodes
            .iter()
            .map(|&n| split_semi(&g.node(n).description))
            .collect(),
    )
    .join(";");
    let lengths: Vec<usize> = nodes
        .iter()
        .flat_map(|&n| g.node(n).lengths.clone())
        .collect();
    // `max()` on a list of (int, str) tuples -- Python tuple comparison, so length first
    // then the centroid id lexicographically. Returns the FIRST maximum.
    let mut long_centroid_id = g.node(nodes[0]).long_centroid_id.clone();
    for &n in &nodes[1..] {
        let c = &g.node(n).long_centroid_id;
        if (c.0, c.1.as_str()) > (long_centroid_id.0, long_centroid_id.1.as_str()) {
            long_centroid_id = c.clone();
        }
    }
    let paralog = nodes.iter().any(|&n| g.node(n).paralog);
    let prev_centroids = if g.node(nodes[0]).prev_centroids.is_some() {
        // ";".join(set(iter_del_dups(...))) -- a Python set of str, so the join order is
        // set order. Only reachable for graphs read back from GML, which is out of scope
        // for phase 1; sorted for determinism if that path is ever taken.
        let mut v = iter_del_dups(
            nodes
                .iter()
                .map(|&n| split_semi(g.node(n).prev_centroids.as_deref().unwrap_or("")))
                .collect(),
        );
        v.sort();
        Some(v.join(";"))
    } else {
        None
    };

    g.add_node(
        new_node,
        NodeAttrs {
            size: members.len(),
            centroid,
            max_len_id,
            members,
            seq_ids,
            has_end: nodes.iter().any(|&n| g.node(n).has_end),
            protein,
            dna,
            annotation,
            description,
            lengths,
            long_centroid_id,
            paralog,
            merged_dna,
            prev_centroids,
            name: None,
            genome_ids: None,
            gene_ids: None,
            degrees: None,
            high_var: None,
            gml_late_attrs_before_name: false,
        },
    );

    // Now iterate through neighbours of each node and add them to the new node
    let merge_nodes: HashSet<usize> = nodes.iter().copied().collect();
    for &node in &nodes {
        for neighbour in g.neighbors(node) {
            if merge_nodes.contains(&neighbour) {
                continue;
            }
            if g.has_edge(new_node, neighbour) {
                let src = g.edge(node, neighbour).members.clone();
                let e = g.edge_mut(new_node, neighbour);
                e.members.union_update(&src);
                e.size = e.members.len();
            } else {
                let src = g.edge(node, neighbour).clone();
                g.add_edge(
                    new_node,
                    neighbour,
                    EdgeAttrs {
                        size: src.size,
                        members: src.members,
                        genome_ids: None,
                    },
                );
            }
        }
    }

    // remove old nodes from Graph
    g.remove_nodes_from(&nodes);
}

/// `value.split(";")` — Python's `"".split(";")` is `[""]`, and Rust agrees.
fn split_semi(v: &str) -> Vec<String> {
    v.split(';').map(|x| x.to_string()).collect()
}

/// `merge_nodes.py::delete_node`
///
/// Reconnects each member's neighbours pairwise before removing the node.
pub fn delete_node(g: &mut Graph, node: usize) {
    // add in new edges
    for mem in g.node(node).members.iter().collect::<Vec<_>>() {
        // Tier D: `sorted(set(...))`, not `list(set(...))`. The order decides adjacency
        // insertion order, which decides BFS order in collapse_families.
        let mut mem_edges: Vec<usize> = g
            .edges_of(&[node])
            .into_iter()
            .filter(|&e| g.edge(e.0, e.1).members.contains(mem))
            .map(|e| e.1)
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        mem_edges.dedup();
        if mem_edges.len() < 2 {
            continue;
        }
        for i in 0..mem_edges.len() {
            for j in (i + 1)..mem_edges.len() {
                let (n1, n2) = (mem_edges[i], mem_edges[j]);
                if g.has_edge(n1, n2) {
                    let e = g.edge_mut(n1, n2);
                    e.members.add(mem);
                    e.size = e.members.len();
                } else {
                    g.add_edge(
                        n1,
                        n2,
                        EdgeAttrs {
                            size: 1,
                            members: IntBitSet::from_iter_ints([mem]),
                            genome_ids: None,
                        },
                    );
                }
            }
        }
    }

    // now remove node
    g.remove_node(node);
}

/// `merge_nodes.py::remove_member_from_node`
pub fn remove_member_from_node(g: &mut Graph, node: usize, member: usize) {
    // add in replacement edges if required -- Tier D sorts this (see delete_node)
    let mem_edges: Vec<usize> = g
        .edges_of(&[node])
        .into_iter()
        .filter(|&e| g.edge(e.0, e.1).members.contains(member))
        .map(|e| e.1)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    if mem_edges.len() > 1 {
        for i in 0..mem_edges.len() {
            for j in (i + 1)..mem_edges.len() {
                let (n1, n2) = (mem_edges[i], mem_edges[j]);
                if g.has_edge(n1, n2) {
                    let e = g.edge_mut(n1, n2);
                    e.members.add(member);
                    e.size = e.members.len();
                } else {
                    g.add_edge(
                        n1,
                        n2,
                        EdgeAttrs {
                            size: 1,
                            members: IntBitSet::from_iter_ints([member]),
                            genome_ids: None,
                        },
                    );
                }
            }
        }
    }

    // remove member from node
    {
        let n = g.node_mut(node);
        n.members.discard(member);
        // `sid.split("_")[0] != str(member)` -- compare the genome field as a string, which
        // matters for refound IDs of the form "{member}_refound_{n}".
        let prefix = member.to_string();
        n.seq_ids = n
            .seq_ids
            .iter()
            .filter(|sid| sid.split('_').next() != Some(prefix.as_str()))
            .cloned()
            .collect();
        n.size -= 1;
    }

    // remove member from edges of node
    let mut edges_to_remove = Vec::new();
    for e in g.edges_of(&[node]) {
        if g.edge(e.0, e.1).members.contains(member) {
            if g.edge(e.0, e.1).members.len() == 1 {
                edges_to_remove.push(e);
            } else {
                let ea = g.edge_mut(e.0, e.1);
                ea.members.discard(member);
                ea.size = ea.members.len();
            }
        }
    }
    for e in edges_to_remove {
        g.remove_edge(e.0, e.1);
    }
}
