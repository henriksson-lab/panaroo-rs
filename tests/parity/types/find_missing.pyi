from networkx.classes.graph import Graph
from typing import List


def find_missing(
    G: Graph,
    gff_file_handles: List[str],
    dna_seq_file: str,
    prot_seq_file: str,
    gene_data_file: str,
    merge_id_thresh: float,
    search_radius: int,
    prop_match: float,
    pairwise_id_thresh: float,
    n_cpu: int,
    remove_by_consensus: bool = ...,
    only_valid_genes: bool = ...,
    verbose: bool = ...
) -> Graph: ...
