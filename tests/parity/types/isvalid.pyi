from intbitset import intbitset
from numpy import int64
from typing import (
    List,
    Set,
    Union,
)


def conv_list(maybe_list: List[str]) -> List[str]: ...


def custom_stringizer(value: Union[str, intbitset, Set[str]]) -> str: ...


def del_dups(seq: List[Union[str, int64]]) -> List[Union[str, int64]]: ...


def is_valid_gene(dna: str, protein: str) -> bool: ...
