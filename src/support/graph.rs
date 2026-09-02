//! `networkx.Graph` — undirected, integer node keys, attribute dicts.
//!
//! Two ordering properties are load-bearing for byte parity and must not be traded away:
//!
//!  1. **Node insertion order.** `for node in G.nodes()` walks the node dict in insertion
//!     order, and that order reaches `pre_filt_graph.gml` / `final_graph.gml` directly.
//!  2. **Adjacency insertion order.** `G.neighbors(n)` walks `adj[n]` in insertion order,
//!     which fixes the traversal order of `nx.bfs_edges` in `collapse_families` and of
//!     `mod_bfs_edges` in `collapse_paralogs` — and therefore which nodes get merged, and
//!     therefore the new node IDs.
//!
//! Both are handled by [`PyDict`]; the rule is that removals go through `PyDict::pop`.
//!
//! Edge attributes are shared between `adj[u][v]` and `adj[v][u]` in networkx (the same
//! dict object). Here they live once in [`Graph::edge_attrs`] keyed by the normalised pair,
//! with `adj` carrying only neighbour order.

use super::intbitset::IntBitSet;
use super::pydict::PyDict;
use std::collections::BTreeSet;

/// Node attributes. Field names are snake_case; the **GML keys stay camelCase and are
/// emitted in the original insertion order** — see `PORTING_PLAN.md` §7.
///
/// Set in `generate_network::generate_network` and `merge_nodes::merge_node_cluster`.
#[derive(Debug, Clone)]
pub struct NodeAttrs {
    /// `size`
    pub size: usize,
    /// `centroid` — list of centroid sequence IDs
    pub centroid: Vec<String>,
    /// `maxLenId` — index into `dna` of the longest sequence
    pub max_len_id: usize,
    /// `members` — genome IDs
    pub members: IntBitSet,
    /// `seqIDs` — a Python `set`.
    ///
    /// A `BTreeSet` rather than a hash set: the Tier D reference wraps every observable
    /// iteration of this set in `sorted()` (PORTING_PLAN.md §6.1), so sorted order *is* the
    /// reference's order. `BTreeSet<String>` orders by UTF-8 bytes, which equals Python's
    /// code-point order for the ASCII sequence IDs Panaroo generates — so `"10_1_2"` sorts
    /// before `"9_1_2"`, matching Python. Do not substitute a natural sort.
    pub seq_ids: BTreeSet<String>,
    /// `hasEnd`
    pub has_end: bool,
    /// `protein`
    pub protein: Vec<String>,
    /// `dna`
    pub dna: Vec<String>,
    /// `annotation` — `;`-joined
    pub annotation: String,
    /// `description` — `;`-joined
    pub description: String,
    /// `lengths`
    pub lengths: Vec<usize>,
    /// `longCentroidID` — `(len(dna), centroid_id)`, compared as a Python tuple
    pub long_centroid_id: (usize, String),
    /// `paralog`
    pub paralog: bool,
    /// `mergedDNA`
    pub merged_dna: bool,

    // --- assigned later in the run; `None` until then -----------------------------------
    /// `prevCentroids` — only present on graphs read back from GML
    pub prev_centroids: Option<String>,
    /// `name` — set by `generate_roary_gene_presence_absence`
    pub name: Option<String>,
    /// `genomeIDs` — set in `__main__::main` before each GML write
    pub genome_ids: Option<String>,
    /// `geneIDs` — ditto
    pub gene_ids: Option<String>,
    /// `degrees` — ditto
    pub degrees: Option<usize>,
    /// `highVar` — set by `identify_possible_highly_variable` (not in the phase-1 path)
    pub high_var: Option<i64>,

    /// Whether `genomeIDs` / `geneIDs` / `degrees` were assigned to this node **before**
    /// `generate_roary_gene_presence_absence` set `name`.
    ///
    /// Not a Python attribute — it records Python dict *insertion order*, which
    /// `nx.write_gml` emits and which therefore has to be byte-reproduced. `main` assigns
    /// those three keys to every node before writing `pre_filt_graph.gml`; nodes created
    /// later by a merge do not have them, so when `main` assigns them again before the
    /// final write they land *after* `name`. Both orders occur in one file — on the `tiny`
    /// dataset, 1194 nodes have `name` last and 16 have it before `genomeIDs`.
    pub gml_late_attrs_before_name: bool,
}

/// Edge attributes: `size` and `members`, plus `genomeIDs` at GML-write time.
#[derive(Debug, Clone)]
pub struct EdgeAttrs {
    pub size: usize,
    pub members: IntBitSet,
    pub genome_ids: Option<String>,
}

/// `G.graph` — graph-level attributes.
#[derive(Debug, Clone, Default)]
pub struct GraphAttrs {
    /// `isolateNames`
    pub isolate_names: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct Graph {
    nodes: PyDict<usize, NodeAttrs>,
    /// node -> (neighbour -> ()), carrying neighbour insertion order only.
    adj: PyDict<usize, PyDict<usize, ()>>,
    /// Keyed by `(min(u,v), max(u,v))`.
    edge_attrs: PyDict<(usize, usize), EdgeAttrs>,
    pub graph: GraphAttrs,
}

impl Default for Graph {
    fn default() -> Self {
        Self::new()
    }
}

impl Graph {
    /// `nx.Graph()`
    pub fn new() -> Self {
        Graph {
            nodes: PyDict::new(),
            adj: PyDict::new(),
            edge_attrs: PyDict::new(),
            graph: GraphAttrs::default(),
        }
    }

    fn ekey(u: usize, v: usize) -> (usize, usize) {
        if u <= v {
            (u, v)
        } else {
            (v, u)
        }
    }

    /// `G.add_node(n, **attrs)`
    ///
    /// networkx *updates* an existing node's attribute dict and leaves its position in the
    /// node order alone. Our attributes are a struct, so "update" is a whole-struct
    /// replacement — which matches, because Panaroo only ever calls `add_node` with the
    /// complete attribute set. Position is preserved by `PyDict::insert`.
    pub fn add_node(&mut self, n: usize, attrs: NodeAttrs) {
        self.nodes.insert(n, attrs);
        if !self.adj.contains_key(&n) {
            self.adj.insert(n, PyDict::new());
        }
    }

    /// `G.has_node(n)`
    pub fn has_node(&self, n: usize) -> bool {
        self.nodes.contains_key(&n)
    }

    /// `G.nodes[n]`
    pub fn node(&self, n: usize) -> &NodeAttrs {
        self.nodes
            .get(&n)
            .unwrap_or_else(|| panic!("KeyError: node {n}"))
    }

    pub fn node_mut(&mut self, n: usize) -> &mut NodeAttrs {
        self.nodes
            .get_mut(&n)
            .unwrap_or_else(|| panic!("KeyError: node {n}"))
    }

    /// `G.nodes()` — insertion order.
    pub fn nodes(&self) -> Vec<usize> {
        self.nodes.keys().copied().collect()
    }

    /// `len(G)` / `len(G.nodes())`
    pub fn number_of_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// `G.remove_node(n)` — also drops it from every neighbour's adjacency.
    pub fn remove_node(&mut self, n: usize) {
        let nbrs: Vec<usize> = match self.adj.get(&n) {
            Some(a) => a.keys().copied().collect(),
            None => panic!("NetworkXError: node {n} is not in the graph"),
        };
        for v in nbrs {
            if let Some(av) = self.adj.get_mut(&v) {
                av.pop(&n);
            }
            self.edge_attrs.pop(&Self::ekey(n, v));
        }
        self.adj.pop(&n);
        self.nodes.pop(&n);
    }

    /// `G.remove_nodes_from(ns)` — silently ignores nodes that are not present.
    pub fn remove_nodes_from(&mut self, ns: &[usize]) {
        for &n in ns {
            if self.has_node(n) {
                self.remove_node(n);
            }
        }
    }

    /// `G.add_edge(u, v, **attrs)`
    ///
    /// networkx would create missing endpoints with empty attribute dicts. Panaroo never
    /// relies on that — every call site adds both nodes first — so a missing endpoint is a
    /// translation error and panics rather than fabricating a node with junk attributes.
    ///
    /// An existing edge keeps its position in both adjacency lists; only the attributes are
    /// replaced. A new edge is appended to `adj[u]` and `adj[v]`, in that order.
    pub fn add_edge(&mut self, u: usize, v: usize, attrs: EdgeAttrs) {
        assert!(
            self.has_node(u),
            "add_edge: node {u} not in graph (networkx would create it; Panaroo never does)"
        );
        assert!(
            self.has_node(v),
            "add_edge: node {v} not in graph (networkx would create it; Panaroo never does)"
        );
        self.adj.get_mut(&u).unwrap().insert(v, ());
        self.adj.get_mut(&v).unwrap().insert(u, ());
        self.edge_attrs.insert(Self::ekey(u, v), attrs);
    }

    /// `G.has_edge(u, v)`
    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        self.edge_attrs.contains_key(&Self::ekey(u, v))
    }

    /// `G[u][v]` / `G.edges[u, v]`
    pub fn edge(&self, u: usize, v: usize) -> &EdgeAttrs {
        self.edge_attrs
            .get(&Self::ekey(u, v))
            .unwrap_or_else(|| panic!("KeyError: edge ({u}, {v})"))
    }

    pub fn edge_mut(&mut self, u: usize, v: usize) -> &mut EdgeAttrs {
        self.edge_attrs
            .get_mut(&Self::ekey(u, v))
            .unwrap_or_else(|| panic!("KeyError: edge ({u}, {v})"))
    }

    /// `G.remove_edge(u, v)`
    pub fn remove_edge(&mut self, u: usize, v: usize) {
        if self.edge_attrs.pop(&Self::ekey(u, v)).is_none() {
            panic!("NetworkXError: edge ({u}, {v}) is not in the graph");
        }
        self.adj.get_mut(&u).unwrap().pop(&v);
        self.adj.get_mut(&v).unwrap().pop(&u);
    }

    /// `G.edges()`
    ///
    /// networkx's `EdgeView.__iter__`: walk nodes in insertion order, yield `(n, nbr)` for
    /// each neighbour not already emitted, then mark `n` seen. So each edge appears once,
    /// oriented from whichever endpoint comes first in node order. Verified against
    /// networkx 3.3.
    pub fn edges(&self) -> Vec<(usize, usize)> {
        self.edges_over(self.nodes.keys().copied())
    }

    /// `G.edges(nbunch)` — same walk, restricted to `ns`, with the `seen` set scoped to
    /// `ns`. Note this is order-sensitive: `G.edges([3,5])` and `G.edges([5,3])` differ.
    pub fn edges_of(&self, ns: &[usize]) -> Vec<(usize, usize)> {
        self.edges_over(ns.iter().copied())
    }

    fn edges_over(&self, it: impl Iterator<Item = usize>) -> Vec<(usize, usize)> {
        let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut out = Vec::new();
        for n in it {
            if let Some(nbrs) = self.adj.get(&n) {
                for &nbr in nbrs.keys() {
                    if !seen.contains(&nbr) {
                        out.push((n, nbr));
                    }
                }
            }
            seen.insert(n);
        }
        out
    }

    /// `G.neighbors(n)` — adjacency insertion order.
    pub fn neighbors(&self, n: usize) -> Vec<usize> {
        match self.adj.get(&n) {
            Some(a) => a.keys().copied().collect(),
            None => panic!("KeyError: node {n}"),
        }
    }

    /// `G.degree[n]`
    pub fn degree(&self, n: usize) -> usize {
        match self.adj.get(&n) {
            Some(a) => a.len(),
            None => panic!("KeyError: node {n}"),
        }
    }

    /// `G.degree()` — `(node, degree)` in node insertion order.
    pub fn degrees(&self) -> Vec<(usize, usize)> {
        self.nodes.keys().map(|&n| (n, self.degree(n))).collect()
    }

    /// `G.subgraph(nodes)` — a view in networkx; an owned copy here.
    ///
    /// The copy preserves the *original* node and adjacency order restricted to `nodes`,
    /// which is what the view does.
    pub fn subgraph(&self, nodes: &[usize]) -> Graph {
        let keep: std::collections::HashSet<usize> = nodes.iter().copied().collect();
        let mut g = Graph::new();
        for &n in self.nodes.keys() {
            if keep.contains(&n) {
                g.add_node(n, self.node(n).clone());
            }
        }
        for &n in self.nodes.keys() {
            if !keep.contains(&n) {
                continue;
            }
            for &nbr in self.adj.get(&n).unwrap().keys() {
                if keep.contains(&nbr) && !g.has_edge(n, nbr) {
                    g.add_edge(n, nbr, self.edge(n, nbr).clone());
                }
            }
        }
        g
    }
}

// --- networkx free functions ------------------------------------------------------------

/// `nx.connected_components(G)` — yields `set[usize]`. networkx builds each component with
/// a plain BFS from the first unvisited node in node insertion order. Callers iterate the
/// resulting set under `sorted()` in the Tier D reference, so a `BTreeSet` matches.
pub fn connected_components(g: &Graph) -> Vec<BTreeSet<usize>> {
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut out = Vec::new();
    for v in g.nodes() {
        if !seen.contains(&v) {
            let c = plain_bfs(g, v);
            seen.extend(c.iter().copied());
            out.push(c);
        }
    }
    out
}

/// networkx's `_plain_bfs` — level-order flood fill returning the component as a set.
fn plain_bfs(g: &Graph, source: usize) -> BTreeSet<usize> {
    let mut seen = BTreeSet::new();
    seen.insert(source);
    let mut nextlevel = vec![source];
    while !nextlevel.is_empty() {
        let thislevel = std::mem::take(&mut nextlevel);
        for v in thislevel {
            for w in g.neighbors(v) {
                if seen.insert(w) {
                    nextlevel.push(w);
                }
            }
        }
    }
    seen
}

/// `nx.is_connected(G)`
pub fn is_connected(g: &Graph) -> bool {
    if g.number_of_nodes() == 0 {
        panic!("NetworkXPointlessConcept: connectivity is undefined for the null graph");
    }
    let first = g.nodes()[0];
    plain_bfs(g, first).len() == g.number_of_nodes()
}

/// `nx.bfs_edges(G, source, depth_limit=None)` — yields `(parent, child)`.
///
/// networkx 3.3's `generic_bfs_edges`: strict level order, a child is marked seen the
/// moment it is yielded. `depth_limit` counts levels, so `depth_limit=1` yields only the
/// source's own neighbours. Verified against networkx 3.3.
pub fn bfs_edges(g: &Graph, source: usize, depth_limit: Option<usize>) -> Vec<(usize, usize)> {
    let depth_limit = depth_limit.unwrap_or_else(|| g.number_of_nodes());
    let n = g.number_of_nodes();
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    seen.insert(source);
    let mut out = Vec::new();
    let mut depth = 0;
    let mut next_parents = vec![source];
    while !next_parents.is_empty() && depth < depth_limit {
        let this_parents = std::mem::take(&mut next_parents);
        for parent in this_parents {
            for child in g.neighbors(parent) {
                if seen.insert(child) {
                    out.push((parent, child));
                    next_parents.push(child);
                }
            }
            if seen.len() == n {
                return out;
            }
        }
        depth += 1;
    }
    out
}

/// `nx.shortest_path_length(G, source, target)`; `None` for `NetworkXNoPath`.
pub fn shortest_path_length(g: &Graph, source: usize, target: usize) -> Option<usize> {
    if source == target {
        return Some(0);
    }
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    seen.insert(source);
    let mut level = 0usize;
    let mut frontier = vec![source];
    while !frontier.is_empty() {
        level += 1;
        let mut next = Vec::new();
        for v in frontier {
            for w in g.neighbors(v) {
                if seen.insert(w) {
                    if w == target {
                        return Some(level);
                    }
                    next.push(w);
                }
            }
        }
        frontier = next;
    }
    None
}

/// `nx.cycle_basis(G, root)`
pub fn cycle_basis(_g: &Graph, _root: usize) -> Vec<Vec<usize>> {
    panic!("noimpl: support::graph::cycle_basis")
}

/// networkx's `LIST_START_VALUE`.
pub const LIST_START_VALUE: &str = "_networkx_list_start";

/// A node's attributes for GML output: its label and its ordered key/value pairs.
pub type GmlNode = (String, Vec<(String, GmlAttr)>);

/// An edge's: source index, target index, and ordered key/value pairs.
pub type GmlEdge = (usize, usize, Vec<(String, GmlAttr)>);

/// One attribute value, in the shape networkx's `stringize` dispatches on.
///
/// Not a Python type — Python passes whatever the attribute holds and networkx
/// `isinstance`-dispatches. This enum is that set of cases.
#[derive(Debug, Clone)]
pub enum GmlAttr {
    Int(i64),
    Bool(bool),
    Float(f64),
    Str(String),
    /// A `list` or `tuple`. Rendered as repeated keys, with a `_networkx_list_start`
    /// marker first when it has exactly one element (so a 1-list is distinguishable from a
    /// scalar) and `"[]"` when empty.
    List(Vec<GmlAttr>),
    /// Anything else — a `set`, an `intbitset`. Only writable when a stringizer is given.
    Opaque(crate::isvalid::GmlValue),
}

/// `nx.write_gml(G, path, stringizer=...)`.
///
/// # Provenance
///
/// Reproduces the GML dialect of [NetworkX](https://networkx.org/) 3.3's
/// `generate_gml` (BSD 3-Clause). No NetworkX source was copied; the rules below were read
/// off its documented behaviour and verified against real output. See `NOTICE.md`.
///
/// # Rules
///
/// - node `id` is the node's **position in insertion order**, 0-based — not the node key.
///   `label` is the node key, always quoted.
/// - `int`/`bool`: `True`/`False` render as `1`/`0`; a value outside signed 32-bit range is
///   quoted; otherwise bare.
/// - `float`: `repr(v).upper()`, with `+` prefixed to `INF` and a `.` inserted before `E`
///   if the mantissa has none.
/// - `list`/`tuple` (and key != `"label"`, not already inside a list): empty renders as
///   `key "[]"`; a **single-element** list emits `key "_networkx_list_start"` first; then
///   each element is emitted under the same key.
/// - everything else goes through the stringizer, then is quoted and escaped.
///
/// **A `str` value falls into that last branch too** — there is no `isinstance(value, str)`
/// case before it. So with a stringizer supplied, strings are passed through it as well,
/// which is why `pre_filt_graph.gml` contains `annotation "'mmpL8_1'"` (repr-quoted inside
/// GML quotes) while `final_graph.gml`, written without a stringizer, has
/// `annotation ""`.
pub fn write_gml(
    path: &str,
    attrs: &[GmlNode],
    edges: &[GmlEdge],
    graph_attrs: &[(String, GmlAttr)],
    stringizer: Option<&dyn Fn(&crate::isvalid::GmlValue) -> String>,
) {
    use std::io::Write;
    let f = std::fs::File::create(path).unwrap_or_else(|e| panic!("create {path}: {e}"));
    let mut out = std::io::BufWriter::new(f);

    writeln!(out, "graph [").unwrap();
    for (k, v) in graph_attrs {
        for line in stringize_attr(k, v, "  ", false, stringizer) {
            writeln!(out, "{line}").unwrap();
        }
    }
    for (id, (label, kvs)) in attrs.iter().enumerate() {
        writeln!(out, "  node [").unwrap();
        writeln!(out, "    id {id}").unwrap();
        writeln!(out, "    label \"{}\"", escape(label)).unwrap();
        for (k, v) in kvs {
            for line in stringize_attr(k, v, "    ", false, stringizer) {
                writeln!(out, "{line}").unwrap();
            }
        }
        writeln!(out, "  ]").unwrap();
    }
    for (src, tgt, kvs) in edges {
        writeln!(out, "  edge [").unwrap();
        writeln!(out, "    source {src}").unwrap();
        writeln!(out, "    target {tgt}").unwrap();
        for (k, v) in kvs {
            for line in stringize_attr(k, v, "    ", false, stringizer) {
                writeln!(out, "{line}").unwrap();
            }
        }
        writeln!(out, "  ]").unwrap();
    }
    writeln!(out, "]").unwrap();
}

fn stringize_attr(
    key: &str,
    value: &GmlAttr,
    indent: &str,
    in_list: bool,
    stringizer: Option<&dyn Fn(&crate::isvalid::GmlValue) -> String>,
) -> Vec<String> {
    match value {
        GmlAttr::Bool(b) => vec![format!("{indent}{key} {}", if *b { 1 } else { 0 })],
        GmlAttr::Int(n) => {
            if *n < -(1i64 << 31) || *n >= (1i64 << 31) {
                vec![format!("{indent}{key} \"{n}\"")]
            } else {
                vec![format!("{indent}{key} {n}")]
            }
        }
        GmlAttr::Float(x) => {
            let mut text = crate::support::pyfmt::py_str_f64(*x).to_uppercase();
            if text == "INF" {
                text = format!("+{text}");
            } else if let Some(epos) = text.rfind('E') {
                if !text[..epos].contains('.') {
                    text = format!("{}.{}", &text[..epos], &text[epos..]);
                }
            }
            vec![format!("{indent}{key} {text}")]
        }
        GmlAttr::List(items) if !in_list => {
            let mut out = Vec::new();
            if items.is_empty() {
                out.push(format!("{indent}{key} \"[]\""));
            }
            if items.len() == 1 {
                out.push(format!("{indent}{key} \"{LIST_START_VALUE}\""));
            }
            for item in items {
                out.extend(stringize_attr(key, item, indent, true, stringizer));
            }
            out
        }
        other => {
            let text = match other {
                GmlAttr::Str(v) => match stringizer {
                    Some(f) => f(&crate::isvalid::GmlValue::Str(v.clone())),
                    None => v.clone(),
                },
                GmlAttr::Opaque(v) => match stringizer {
                    Some(f) => f(v),
                    None => panic!("NetworkXError: value is not a string"),
                },
                // a nested list reached with in_list=true is emitted element-wise above,
                // so this arm only sees a list when in_list is set -- treat it as opaque
                GmlAttr::List(_) => panic!("NetworkXError: nested list in GML"),
                _ => unreachable!(),
            };
            vec![format!("{indent}{key} \"{}\"", escape(&text))]
        }
    }
}

/// networkx's `escape`: XML character references for anything outside printable ASCII,
/// plus `&` and `"`. `re.sub('[^ -~]|[&"]', ...)`.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        let cp = c as u32;
        if !(0x20..=0x7e).contains(&cp) || c == '&' || c == '"' {
            out.push_str(&format!("&#{cp};"));
        } else {
            out.push(c);
        }
    }
    out
}
