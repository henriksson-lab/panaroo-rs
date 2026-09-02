#!/usr/bin/env python3
"""Replay tests/parity/graph/ops.txt against networkx and dump the results.

Pairs with `cargo run --example dump_graph`. Any difference is a bug in support::graph.
"""
import sys, os, networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build", "reference"))
from panaroo.clean_network import mod_bfs_edges

G = nx.Graph()
out = []
for raw in open(sys.argv[1]):
    line = raw.split("#")[0].strip()
    if not line:
        continue
    parts = line.split()
    op, args = parts[0], [int(x) for x in parts[2:]] if parts[0] == "dump" else [int(x) for x in parts[1:]]
    if op == "add_node":
        G.add_node(args[0], size=args[0])
    elif op == "add_edge":
        G.add_edge(args[0], args[1], size=1)
    elif op == "remove_edge":
        G.remove_edge(args[0], args[1])
    elif op == "remove_node":
        G.remove_node(args[0])
    elif op == "dump":
        what = parts[1]
        if what == "nodes":       r = list(G.nodes())
        elif what == "edges":     r = [list(e) for e in G.edges()]
        elif what == "neighbors": r = list(G.neighbors(args[0]))
        elif what == "edges_of":  r = [list(e) for e in G.edges(args)]
        elif what == "degrees":   r = [list(d) for d in G.degree()]
        elif what == "bfs":
            dl = args[1] if len(args) > 1 else None
            r = [list(e) for e in nx.bfs_edges(G, args[0], depth_limit=dl)]
        elif what == "mod_bfs":
            dl = args[1] if len(args) > 1 else None
            r = [list(e) for e in mod_bfs_edges(G, args[0], depth_limit=dl)]
        elif what == "components": r = [sorted(c) for c in nx.connected_components(G)]
        elif what == "spl":
            try: r = nx.shortest_path_length(G, args[0], args[1])
            except nx.NetworkXNoPath: r = None
        else: raise SystemExit(f"unknown dump: {what}")
        out.append(f"{line} -> {r}")
    else:
        raise SystemExit(f"unknown op: {op}")
print("\n".join(out))
