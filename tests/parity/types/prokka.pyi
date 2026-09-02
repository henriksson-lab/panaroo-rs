from Bio.SeqRecord import SeqRecord
from io import TextIOWrapper
from collections import OrderedDict
from numpy import ndarray
from typing import (
    List,
    Set,
    Union,
)


def get_trans_table(table: int) -> List[Union[ndarray, Set[str]]]: ...


def output_files(
    dna_dictionary: OrderedDict,
    protien_list: List[SeqRecord],
    prot_handle: TextIOWrapper,
    dna_handle: TextIOWrapper,
    csv_handle: TextIOWrapper,
    gff_filename: str
) -> None: ...


def process_prokka_input(
    gff_list: List[str],
    output_dir: str,
    filter_seqs: bool,
    quiet: bool,
    n_cpu: int,
    table: int
) -> bool: ...
