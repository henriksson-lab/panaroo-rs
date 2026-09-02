from networkx.classes.graph import Graph
from typing import (
    Dict,
    Tuple,
)


def generate_common_struct_presence_absence(
    G: Graph,
    output_dir: str,
    mems_to_isolates: Dict[int, str],
    min_variant_support: int = ...
) -> None: ...


def generate_pan_genome_reference(
    G: Graph,
    output_dir: str,
    ids_len_stop: Dict[str, Tuple[int, bool, bool]],
    split_paralogs: bool = ...
) -> None: ...


def generate_roary_gene_presence_absence(
    G: Graph,
    mems_to_isolates: Dict[int, str],
    orig_ids: Dict[str, str],
    ids_len_stop: Dict[str, Tuple[int, bool, bool]],
    output_dir: str
) -> Graph: ...


def generate_summary_stats(output_dir: str) -> bool: ...
