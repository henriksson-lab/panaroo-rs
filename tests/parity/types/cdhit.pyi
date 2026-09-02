from typing import (
    Optional,
    Union,
)


def check_cdhit_version(cdhit_exec: str = ...) -> float: ...


def run_cdhit(
    input_file: str,
    output_file: str,
    id: float = ...,
    n_cpu: int = ...,
    s: Optional[float] = ...,
    aL: float = ...,
    AL: int = ...,
    aS: Union[float, int] = ...,
    AS: int = ...,
    accurate: bool = ...,
    use_local: bool = ...,
    word_length: None = ...,
    min_length: None = ...,
    quiet: bool = ...
) -> None: ...


def run_cdhit_est(
    input_file: str,
    output_file: str,
    id: float = ...,
    n_cpu: int = ...,
    s: None = ...,
    aL: float = ...,
    AL: int = ...,
    aS: int = ...,
    AS: int = ...,
    accurate: bool = ...,
    use_local: bool = ...,
    strand: int = ...,
    print_aln: bool = ...,
    word_length: Optional[int] = ...,
    mask: bool = ...,
    quiet: bool = ...
) -> None: ...
