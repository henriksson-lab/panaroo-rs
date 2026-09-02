from intbitset import intbitset
from networkx.classes.graph import Graph
from networkx.classes.reportviews import EdgeDataView
from numpy import int64
from typing import (
    Any,
    Iterator,
    List,
    Optional,
)


def delete_node(G: Graph, node: int) -> Graph: ...


def gen_edge_iterables(
    G: Graph,
    edges: EdgeDataView,
    feature: str
) -> Iterator[intbitset]: ...


def gen_node_iterables(
    G: Graph,
    nodes: List[int64],
    feature: str,
    split: Optional[str] = ...
) -> Iterator[Any]: ...


def iter_del_dups(iterable: Iterator[Any]) -> List[str]: ...


def merge_node_cluster(
    G: Graph,
    nodes: List[int64],
    newNode: int,
    multi_centroid: bool = ...,
    check_merge_mems: bool = ...
) -> Graph: ...


def remove_member_from_node(G: Graph, node: int, member: int) -> Graph: ...
