//! Replay tests/parity/graph/ops.txt against support::graph and dump the results.
//!
//! Pairs with `tests/parity/graph/dump_networkx.py`, which emits the identical format from
//! real networkx. Any difference is a bug in support::graph.
//!
//!     cargo run --release --example dump_graph -- tests/parity/graph/ops.txt
use panaroo::support::graph::{
    bfs_edges, connected_components, shortest_path_length, EdgeAttrs, Graph, NodeAttrs,
};
use panaroo::support::intbitset::IntBitSet;
use std::collections::BTreeSet;

fn node_attrs(size: usize) -> NodeAttrs {
    NodeAttrs {
        size,
        centroid: vec![],
        max_len_id: 0,
        members: IntBitSet::default(),
        seq_ids: BTreeSet::new(),
        has_end: false,
        protein: vec![],
        dna: vec![],
        annotation: String::new(),
        description: String::new(),
        lengths: vec![],
        long_centroid_id: (0, String::new()),
        paralog: false,
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

/// Render like Python's `repr` of a list, so the two dumps compare byte-for-byte.
fn pylist(v: &[usize]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    )
}
fn pypairs(v: &[(usize, usize)]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|(a, b)| format!("[{a}, {b}]"))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: dump_graph OPS.txt");
    let text = std::fs::read_to_string(&path).expect("read");
    let mut g = Graph::new();

    for raw in text.lines() {
        let line = raw.split('#').next().unwrap().trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        let nums = |from: usize| -> Vec<usize> {
            parts[from..].iter().map(|s| s.parse().unwrap()).collect()
        };
        match parts[0] {
            "add_node" => {
                let a = nums(1);
                g.add_node(a[0], node_attrs(a[0]));
            }
            "add_edge" => {
                let a = nums(1);
                g.add_edge(
                    a[0],
                    a[1],
                    EdgeAttrs {
                        size: 1,
                        members: IntBitSet::default(),
                        genome_ids: None,
                    },
                );
            }
            "remove_edge" => {
                let a = nums(1);
                g.remove_edge(a[0], a[1]);
            }
            "remove_node" => {
                let a = nums(1);
                g.remove_node(a[0]);
            }
            "dump" => {
                let a = nums(2);
                let r = match parts[1] {
                    "nodes" => pylist(&g.nodes()),
                    "edges" => pypairs(&g.edges()),
                    "neighbors" => pylist(&g.neighbors(a[0])),
                    "edges_of" => pypairs(&g.edges_of(&a)),
                    "degrees" => pypairs(&g.degrees()),
                    "bfs" => pypairs(&bfs_edges(&g, a[0], a.get(1).copied())),
                    "mod_bfs" => format!(
                        "[{}]",
                        panaroo::clean_network::mod_bfs_edges(&g, a[0], a.get(1).copied())
                            .iter()
                            .map(|(p, c, d)| format!("[{p}, {c}, {d}]"))
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                    "components" => format!(
                        "[{}]",
                        connected_components(&g)
                            .iter()
                            .map(|c| pylist(&c.iter().copied().collect::<Vec<_>>()))
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                    "spl" => match shortest_path_length(&g, a[0], a[1]) {
                        Some(d) => d.to_string(),
                        None => "None".to_string(),
                    },
                    other => panic!("unknown dump: {other}"),
                };
                println!("{line} -> {r}");
            }
            other => panic!("unknown op: {other}"),
        }
    }
}
