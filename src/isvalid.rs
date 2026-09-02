//! Translation of `panaroo/panaroo/isvalid.py`.
//!
//! # Vendored third-party code
//!
//! `custom_stringizer` / `stringize` are **not Panaroo's own**. Panaroo's docstring says so:
//! they are a copy of networkx's `literal_stringizer`
//! ([NetworkX](https://networkx.org/), BSD 3-Clause), extended with `set` and `intbitset`
//! branches. They live here because they live in `isvalid.py` and this port is
//! file-for-file, but the origin is recorded here and in `NOTICE.md` per the project rule
//! that third-party code carries its source.
//!
//! The rest of the Python module — the argparse path validators — belongs to the auxiliary
//! entry points and is out of scope.

use std::collections::HashSet;
use std::hash::Hash;

/// `isvalid.py::conv_list` — wraps a scalar in a list.
///
/// ```python
/// def conv_list(maybe_list):
///     if not isinstance(maybe_list, list):
///         maybe_list = [maybe_list]
///     return (maybe_list)
/// ```
///
/// On the `__main__::main` path (lines 548-549) it is only ever called with a `list[str]`
/// — MonkeyType records `List[str] -> List[str]` — so it is the identity there. The
/// scalar arm exists for graphs read back from GML, where `dna`/`protein` are `;`-joined
/// strings; see [`conv_list_scalar`].
///
/// Duplicated verbatim as [`crate::get_neighborhood::conv_list`]; both are kept per
/// PORTING_PLAN.md §5.
pub fn conv_list(maybe_list: Vec<String>) -> Vec<String> {
    maybe_list
}

/// The `not isinstance(maybe_list, list)` arm of `isvalid.py::conv_list`.
///
/// Rust cannot dispatch on runtime type, so the two arms are separate entry points. This is
/// a spelling difference, not a behavioural one — the Python picks the arm by type, and
/// each Rust call site already knows which it has.
pub fn conv_list_scalar(maybe_list: String) -> Vec<String> {
    vec![maybe_list]
}

/// `isvalid.py::del_dups` — deduplicate in place, keeping first-occurrence order.
///
/// ```python
/// def del_dups(seq):
///     seen = set()
///     pos = 0
///     for item in seq:
///         if item not in seen:
///             seen.add(item)
///             seq[pos] = item
///             pos += 1
///     del seq[pos:]
///     return (seq)
/// ```
///
/// Generic because it genuinely is: MonkeyType records `List[Union[str, int64]]`, and both
/// arms are live — `find_missing.py:192,196` passes `list[str]` (dna/protein), while
/// `clean_network.py:83` passes node IDs as numpy `int64`.
///
/// Distinct from [`crate::merge_nodes::del_dups`], which is a different implementation of
/// the same idea (it builds a dict and returns its keys). Both are kept.
pub fn del_dups<T: Eq + Hash + Clone>(seq: &mut Vec<T>) {
    let mut seen: HashSet<T> = HashSet::new();
    let mut pos = 0;
    for i in 0..seq.len() {
        if !seen.contains(&seq[i]) {
            seen.insert(seq[i].clone());
            seq.swap(pos, i);
            pos += 1;
        }
    }
    seq.truncate(pos);
}

/// `isvalid.py::is_valid_gene`
///
/// ```python
/// def is_valid_gene(dna, protein):
///     if len(dna) % 3 != 0:
///         return False
///     protein = protein.strip("X")
///     if protein[0] != 'M':
///         return False
///     if "*" in protein[:-1]:
///         return False
///     return True
/// ```
///
/// Note `protein[0]` after `strip("X")`: an all-`X` or empty protein raises `IndexError` in
/// Python rather than returning `False`. That is preserved as a panic — a caller relying on
/// `False` there would be relying on behaviour the Python does not have.
///
/// `len(dna)` counts bytes here because Panaroo's sequences are ASCII; Python counts code
/// points, which agrees for any input that reaches this function.
pub fn is_valid_gene(dna: &str, protein: &str) -> bool {
    // Check if sequence is divisible by 3
    if dna.len() % 3 != 0 {
        return false;
    }

    // Check for start codon
    let protein = protein.trim_matches('X');
    if protein.is_empty() {
        panic!("is_valid_gene: IndexError: string index out of range");
    }
    if !protein.starts_with('M') {
        return false;
    }

    // Check for premature stop
    let upto_last = &protein[..protein.len() - protein.chars().last().unwrap().len_utf8()];
    if upto_last.contains('*') {
        return false;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values from running the reference Python; see the oracle in the commit
    // message for this batch.

    #[test]
    fn del_dups_keeps_first_occurrence_order() {
        // python: del_dups(['b','a','b','c','a']) -> ['b', 'a', 'c']
        let mut v: Vec<String> = ["b", "a", "b", "c", "a"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        del_dups(&mut v);
        assert_eq!(v, ["b", "a", "c"]);

        // python: del_dups(['x']) -> ['x']
        let mut v: Vec<String> = vec!["x".to_string()];
        del_dups(&mut v);
        assert_eq!(v, ["x"]);

        // python: del_dups([]) -> []
        let mut v: Vec<String> = vec![];
        del_dups(&mut v);
        assert!(v.is_empty());
    }

    #[test]
    fn del_dups_works_on_integers_too() {
        // python: del_dups([3,1,3,2,1]) -> [3, 1, 2]
        // This arm is live: clean_network.py:83 passes numpy int64 node IDs.
        let mut v: Vec<i64> = vec![3, 1, 3, 2, 1];
        del_dups(&mut v);
        assert_eq!(v, [3, 1, 2]);
    }

    #[test]
    fn conv_list_is_identity_on_a_list() {
        // python: conv_list(['a','b']) -> ['a', 'b'];  conv_list('a') -> ['a']
        assert_eq!(conv_list(vec!["a".into(), "b".into()]), ["a", "b"]);
        assert_eq!(conv_list_scalar("a".into()), ["a"]);
    }

    #[test]
    fn is_valid_gene_matches_python() {
        // Each case run through the reference Python.
        assert!(is_valid_gene("ATGAAA", "MK")); // -> True
        assert!(!is_valid_gene("ATGA", "MK")); // -> False  (len % 3 != 0)
        assert!(is_valid_gene("ATGAAA", "XXMKXX")); // -> True   (X stripped both ends)
        assert!(!is_valid_gene("ATGAAA", "K")); // -> False  (no start codon)
        assert!(!is_valid_gene("ATGAAA", "M*K")); // -> False  (premature stop)
        assert!(is_valid_gene("ATGAAA", "MK*")); // -> True   (terminal stop is fine)
        assert!(is_valid_gene("ATGAAA", "M")); // -> True   (protein[:-1] is "")
    }

    #[test]
    #[should_panic(expected = "IndexError")]
    fn is_valid_gene_raises_on_all_x_protein() {
        // python: is_valid_gene("ATGAAA", "XXX") -> IndexError
        is_valid_gene("ATGAAA", "XXX");
    }
}

// --- vendored from networkx: literal_stringizer ------------------------------------------

/// `isvalid.py::custom_stringizer`
///
/// The `stringizer=` callback handed to `nx.write_gml` for `pre_filt_graph.gml`. Converts a
/// value to its GML literal representation.
///
/// **Origin: networkx's `literal_stringizer`** (BSD 3-Clause), copied into Panaroo and
/// extended with branches for `set` and `intbitset`. See the module header and `NOTICE.md`.
///
/// This was invisible to the first call-graph pass: `nx.write_gml(G, path,
/// stringizer=custom_stringizer)` passes it as a *value*, so there is no syntactic call.
/// `tools/port_order.py` now walks every call's arguments for bare identifiers naming known
/// functions, which is what surfaced it — the same class of miss as the joblib `delayed(f)`
/// case.
///
/// Tier D note: the `set` branch is patched to iterate `sorted(value)`. Without that,
/// `seqIDs` renders in CPython set order and `pre_filt_graph.gml` is not reproducible.
/// See `ORIGINAL_CODE_BUG.md` B6.
pub fn custom_stringizer(value: &GmlValue) -> String {
    let mut buf = String::new();
    stringize(value, &mut buf);
    buf
}

/// `isvalid.py::stringize` — the recursive worker nested inside `custom_stringizer`.
///
/// Kept as its own function because the Python defines it as one (rule 1), even though it
/// is a closure there.
pub fn stringize(value: &GmlValue, buf: &mut String) {
    match value {
        // `isinstance(value, (int, long, bool))` -- True/False become 1/0
        GmlValue::Bool(b) => buf.push(if *b { '1' } else { '0' }),
        GmlValue::Int(n) => buf.push_str(&crate::support::pyfmt::py_str_i64(*n)),
        GmlValue::None => buf.push_str("None"),
        // `repr(value)` for str
        GmlValue::Str(v) => buf.push_str(&py_repr_str(v)),
        GmlValue::Float(x) => buf.push_str(&crate::support::pyfmt::py_str_f64(*x)),
        // lists, sets and intbitsets all render with [] and no spaces after commas
        GmlValue::List(items) | GmlValue::Set(items) => {
            buf.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    buf.push(',');
                }
                stringize(item, buf);
            }
            buf.push(']');
        }
        GmlValue::IntBitSet(items) => {
            buf.push('[');
            for (i, n) in items.iter().enumerate() {
                if i > 0 {
                    buf.push(',');
                }
                buf.push_str(&crate::support::pyfmt::py_str_i64(*n));
            }
            buf.push(']');
        }
        GmlValue::Tuple(items) => {
            if items.len() > 1 {
                buf.push('(');
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        buf.push(',');
                    }
                    stringize(item, buf);
                }
                buf.push(')');
            } else if items.len() == 1 {
                buf.push('(');
                stringize(&items[0], buf);
                buf.push_str(",)");
            } else {
                buf.push_str("()");
            }
        }
    }
}

/// Python `repr()` of a `str`.
///
/// Single quotes unless the string contains a `'` and no `"`. Escapes `\`, the quote
/// character, and the control characters `\n \r \t`; other non-printables become `\xNN`.
/// Non-ASCII printable characters are kept verbatim in Python 3 — networkx's `escape` then
/// turns them into XML character references, so they never reach the file raw.
fn py_repr_str(s: &str) -> String {
    let has_single = s.contains('\'');
    let has_double = s.contains('"');
    let quote = if has_single && !has_double { '"' } else { '\'' };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(quote);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// The dynamically-typed value `stringize` dispatches on.
///
/// Not a Python type — Python passes whatever the attribute holds and the function
/// `isinstance`-dispatches. This enum is that set of cases.
#[derive(Debug, Clone)]
pub enum GmlValue {
    /// `int`, `long`, `bool` (written as 1/0) or `None`
    Int(i64),
    Bool(bool),
    None,
    Str(String),
    Float(f64),
    List(Vec<GmlValue>),
    /// A Python `set` — Tier D renders this sorted.
    Set(Vec<GmlValue>),
    /// An `intbitset` — always ascending.
    IntBitSet(Vec<i64>),
    Tuple(Vec<GmlValue>),
}
