//! Infrastructure with **no Python counterpart**.
//!
//! Each module here stands in for a third-party library that Panaroo depends on. The
//! contract is behavioural parity with the Python library, not API beauty: these types
//! exist so that the translated code in the parent modules can be a line-by-line
//! transcription of the Python.
//!
//! | module | replaces |
//! |---|---|
//! | [`align_io`]   | `Bio.AlignIO` |
//! | [`biocode_gff3`]| `biocode.things.*.print_as(format='gff3')` and `biocode.utils` |
//! | [`codon_table`]| `Bio.Data.CodonTable` — NCBI genetic code **data**, generated from Biopython |
//! | [`edlib`]      | the `edlib` Python binding (FFI to the same C library) |
//! | [`genbank`]    | `Bio.SeqIO` GenBank flat-file reader |
//! | [`gff`]        | `gffutils` |
//! | [`graph`]      | `networkx.Graph` |
//! | [`intbitset`]  | `intbitset.intbitset` |
//! | [`npmath`]     | `numpy` reductions, with numpy's exact summation order |
//! | [`parallel`]   | `joblib.Parallel` / `joblib.delayed` |
//! | [`proc`]       | `subprocess.run` / `subprocess.Popen` |
//! | [`pydict`]     | CPython `dict` (insertion-ordered) |
//! | [`pyfmt`]      | CPython `str()` / `repr()` number formatting |
//! | [`seq`]        | `Bio.Seq` |
//! | [`seqio`]      | `Bio.SeqIO` |
//! | [`sparse`]     | `scipy.sparse.csr_matrix` + `scipy.sparse.csgraph` |

pub mod align_io;
pub mod biocode_gff3;
pub mod codon_table;
pub mod edlib;
pub mod genbank;
pub mod gff;
pub mod graph;
pub mod intbitset;
pub mod npmath;
pub mod parallel;
pub mod proc;
pub mod pydict;
pub mod pyfmt;
pub mod pytempfile;
pub mod seq;
pub mod seqio;
pub mod sparse;
