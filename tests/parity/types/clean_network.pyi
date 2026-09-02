from networkx.classes.graph import Graph
from numpy import int64
from scipy.sparse._csr import csr_matrix
from typing import (
    Any,
    DefaultDict,
    Dict,
    List,
    Optional,
    Tuple,
    Union,
)


def clean_misassembly_edges(
    G: Graph,
    edge_support_threshold: int
) -> Graph: ...


def collapse_families(
    G: Graph,
    seqid_to_centroid: Dict[str, str],
    outdir: str,
    family_threshold: float = ...,
    dna_error_threshold: float = ...,
    family_len_dif_percent: None = ...,
    correct_mistranslations: bool = ...,
    length_outlier_support_proportion: float = ...,
    n_cpu: int = ...,
    quiet: bool = ...,
    distances_bwtn_centroids: Optional[csr_matrix] = ...,
    centroid_to_index: Optional[Dict[str, int]] = ...,
    depths: List[int] = ...,
    search_genome_ids: None = ...
) -> Tuple[Graph, csr_matrix, Dict[str, int]]: ...


def collapse_paralogs(
    G: Graph,
    centroid_contexts: DefaultDict[Any, Any],
    max_context: int = ...,
    quiet: bool = ...
) -> Graph: ...


def single_linkage(
    G: Graph,
    distances_bwtn_centroids: csr_matrix,
    centroid_to_index: Dict[str, int],
    neighbours: List[Union[int, int64]]
) -> List[List[int64]]: ...


def trim_low_support_trailing_ends(
    G: Graph,
    min_support: int = ...,
    max_recursive: int = ...
) -> Graph: ...
