from networkx.classes.graph import Graph
from typing import (
    Any,
    DefaultDict,
    Dict,
    Tuple,
)


def generate_network(
    cluster_file: str,
    data_file: str,
    prot_seq_file: str,
    all_dna: bool = ...
) -> Tuple[Graph, DefaultDict[Any, Any], Dict[str, str]]: ...
