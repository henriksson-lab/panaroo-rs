//! CPython number formatting and integer semantics.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [CPython](https://python.org/) 3.11 (`float.__repr__`, `int.__floordiv__`,
//! `int.__mod__`, `math.ceil`), Python Software Foundation License 2.0. No CPython source
//! was copied; the rules below were established by running the interpreter. See
//! `NOTICE.md`.
//!
//! # Why this exists
//!
//! Python `str(float)` and Rust `{}` on `f64` both emit shortest-round-trip digits, but
//! they disagree about when to switch to exponent form. Every float that reaches an output
//! file must go through [`py_str_f64`] — in Panaroo that is at minimum the
//! `Avg sequences per isolate` and `Avg group size nuc` columns of
//! `gene_presence_absence_roary.csv`, and the numeric flags in the cd-hit command lines.
//!
//! | value | Python `str()` | Rust `{}` |
//! |---|---|---|
//! | `0.00001` | `1e-05` | `0.00001` |
//! | `1e16` | `1e+16` | `10000000000000000` |
//! | `1e15` | `1000000000000000.0` | `1000000000000000` |
//! | `1.0` | `1.0` | `1` |

/// `str(x)` / `repr(x)` for a Python `float`.
///
/// CPython formats the shortest round-tripping digit string, then chooses a layout from the
/// decimal point position `decpt` (the value is `0.<digits> * 10^decpt`):
///
///   - **exponent form** iff `decpt <= -4 || decpt > 16`
///   - otherwise **fixed form**, always with at least one fractional digit
///
/// In exponent form the exponent carries a sign and at least two digits, and a single-digit
/// mantissa gets no `.0`. OBSERVED against CPython 3.11:
/// `1e16 -> "1e+16"`, `1.5e16 -> "1.5e+16"`, `1e-5 -> "1e-05"`, `1e-10 -> "1e-10"`,
/// `1e15 -> "1000000000000000.0"`, `5e-324 -> "5e-324"`, `-0.0 -> "-0.0"`.
pub fn py_str_f64(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }

    let sign = if x.is_sign_negative() { "-" } else { "" };
    let a = x.abs();
    if a == 0.0 {
        return format!("{sign}0.0");
    }

    // Rust's `{:e}` gives the shortest round-tripping mantissa in scientific form, which is
    // the same digit string CPython computes.
    let sci = format!("{:e}", a); // e.g. "1e16", "1.5e-5"
    let (mant, exp) = sci.split_once('e').expect("scientific form");
    let exp: i32 = exp.parse().expect("exponent");
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let decpt = exp + 1; // value == 0.<digits> * 10^decpt

    if decpt <= -4 || decpt > 16 {
        // exponent form
        let e = decpt - 1;
        let head = &digits[..1];
        let rest = &digits[1..];
        let mantissa = if rest.is_empty() {
            head.to_string()
        } else {
            format!("{head}.{rest}")
        };
        return format!(
            "{sign}{mantissa}e{}{:02}",
            if e < 0 { '-' } else { '+' },
            e.abs()
        );
    }

    // fixed form
    if decpt <= 0 {
        // 0.000ddd
        format!("{sign}0.{}{}", "0".repeat((-decpt) as usize), digits)
    } else if (decpt as usize) >= digits.len() {
        // ddd000.0
        format!(
            "{sign}{}{}.0",
            digits,
            "0".repeat(decpt as usize - digits.len())
        )
    } else {
        let (i, f) = digits.split_at(decpt as usize);
        format!("{sign}{i}.{f}")
    }
}

/// `str(x)` for a numpy scalar reaching `",".join(map(str, ...))`.
///
/// numpy float64 scalars stringify through the same shortest-repr path as Python floats, so
/// this delegates. It exists so call sites record which library produced the value.
pub fn py_str_np_float64(x: f64) -> String {
    py_str_f64(x)
}

/// `str(n)` for a Python `int`. Arbitrary precision upstream, but every value Panaroo
/// prints fits in `i64`.
pub fn py_str_i64(n: i64) -> String {
    n.to_string()
}

/// Python `%` on integers: the result takes the sign of the **divisor**.
///
/// Rust's `%` truncates instead, so `-7 % 2` is `-1` in Rust and `1` in Python.
/// OBSERVED: `-7 % 2 == 1`.
pub fn py_mod_i64(a: i64, b: i64) -> i64 {
    let r = a % b;
    if r != 0 && (r < 0) != (b < 0) {
        r + b
    } else {
        r
    }
}

/// Python `//` on integers: floor division. OBSERVED: `-7 // 2 == -4`.
pub fn py_floordiv_i64(a: i64, b: i64) -> i64 {
    let q = a / b;
    let r = a % b;
    if r != 0 && (r < 0) != (b < 0) {
        q - 1
    } else {
        q
    }
}

/// A Python number that may be an `int` or a `float`, because the value's *type* changes
/// how it prints and Panaroo relies on that.
///
/// `cdhit.py`'s `aL`/`aS` parameters default to the float `0.0`, but `iterative_cdhit`
/// calls `run_cdhit(..., aS=AS)` — passing the **int** `99999999` from the `AS` control
/// parameter (lines 333, 417, 431; almost certainly a slip, but it is what upstream does).
/// So the same parameter reaches `str()` as `"0.0"` in one call and `"99999999"` in
/// another. Formatting it as a float either way would emit `-aS 99999999.0`, a different
/// cd-hit flag, a different clustering, and a different pangenome.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PyNum {
    Int(i64),
    Float(f64),
    /// Python `None`, which `str()` renders as `"None"`.
    ///
    /// Reachable: `--family_len_dif_percent` has `type=float` and **no** `default=`, and
    /// `set_default_args` never fills it, so a default run really does invoke
    /// `cd-hit ... -s None`. (cd-hit's `atof` parses that as 0.0, so the numeric effect is
    /// nil — but the command string is what we are reproducing.)
    None,
}

impl PyNum {
    /// `str(x)` for whichever type it holds.
    pub fn to_py_string(self) -> String {
        match self {
            PyNum::Int(n) => py_str_i64(n),
            PyNum::Float(x) => py_str_f64(x),
            PyNum::None => "None".to_string(),
        }
    }

    /// `Option<f64>` as Python would pass it: the value, or `None`.
    pub fn from_opt_f64(v: Option<f64>) -> Self {
        match v {
            Some(x) => PyNum::Float(x),
            Option::None => PyNum::None,
        }
    }
}

impl std::fmt::Display for PyNum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.to_py_string())
    }
}

/// `math.ceil(x)` — returns an `int` in Python 3.
pub fn py_ceil(x: f64) -> i64 {
    x.ceil() as i64
}

/// `repr()` of a Python `bytes` object.
///
/// Needed because `cdhit::check_cdhit_version` and
/// `generate_alignments::check_aligner_install` regex over
/// `str(subprocess.run(...))`, whose text embeds `stdout=b'...'`. See
/// [`super::proc::run_shell_capture_repr`].
///
/// CPython prefers single quotes, escapes `\\ \n \r \t`, escapes `'` only when the string
/// also contains no `"`, and renders any other non-printable or non-ASCII byte as `\xNN`.
pub fn py_repr_bytes(b: &[u8]) -> String {
    let has_single = b.contains(&b'\'');
    let has_double = b.contains(&b'"');
    let quote = if has_single && !has_double { '"' } else { '\'' };
    let mut out = String::with_capacity(b.len() + 3);
    out.push('b');
    out.push(quote);
    for &c in b {
        match c {
            b'\\' => out.push_str("\\\\"),
            b'\n' => out.push_str("\\n"),
            b'\r' => out.push_str("\\r"),
            b'\t' => out.push_str("\\t"),
            c if c as char == quote => {
                out.push('\\');
                out.push(quote);
            }
            0x20..=0x7e => out.push(c as char),
            _ => out.push_str(&format!("\\x{c:02x}")),
        }
    }
    out.push(quote);
    out
}

/// `repr()` of a list of `str`, as `print` renders it.
///
/// Needed because `__main__.py::main` does `print("Problem reading input line: ", line)`
/// where `line` has already been rebound to `line.strip().split()` -- so the user sees the
/// LIST, `['a', 'b', 'c']`, not the original text.
///
/// Python's `repr` of a str prefers single quotes, switching to double quotes only when the
/// string contains a `'` and no `"`. Backslashes and control characters are escaped.
pub fn py_repr_str_list(items: &[&str]) -> String {
    let mut out = String::from("[");
    for (i, s) in items.iter().enumerate() {
        if i > 0 {
            out.push_str(", ");
        }
        out.push_str(&py_repr_str(s));
    }
    out.push(']');
    out
}

/// `repr()` of a single `str`.
pub fn py_repr_str(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
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
                out.push(c);
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every expectation is a value printed by CPython 3.11.

    #[test]
    fn float_repr_matches_python() {
        assert_eq!(py_str_f64(1.0), "1.0");
        assert_eq!(py_str_f64(0.1), "0.1");
        assert_eq!(py_str_f64(0.98), "0.98");
        assert_eq!(py_str_f64(100.0), "100.0");
        assert_eq!(py_str_f64(2.3333333333333335), "2.3333333333333335");
        assert_eq!(py_str_f64(0.0001), "0.0001");
        assert_eq!(py_str_f64(0.0), "0.0");
        assert_eq!(py_str_f64(-0.0), "-0.0");
        assert_eq!(py_str_f64(-1.5), "-1.5");
    }

    #[test]
    fn float_repr_switches_to_exponent_where_python_does() {
        assert_eq!(py_str_f64(1e-5), "1e-05");
        assert_eq!(py_str_f64(1.5e-5), "1.5e-05");
        assert_eq!(py_str_f64(1e-10), "1e-10");
        assert_eq!(py_str_f64(1e16), "1e+16");
        assert_eq!(py_str_f64(1.5e16), "1.5e+16");
        assert_eq!(py_str_f64(1e17), "1e+17");
        assert_eq!(py_str_f64(5e-324), "5e-324");
        assert_eq!(
            py_str_f64(1.7976931348623157e308),
            "1.7976931348623157e+308"
        );
        // just below the threshold -- stays fixed
        assert_eq!(py_str_f64(1e15), "1000000000000000.0");
        assert_eq!(py_str_f64(9999999999999998.0), "9999999999999998.0");
        assert_eq!(py_str_f64(1234567890123456.8), "1234567890123456.8");
    }

    #[test]
    fn float_repr_specials() {
        assert_eq!(py_str_f64(f64::NAN), "nan");
        assert_eq!(py_str_f64(f64::INFINITY), "inf");
        assert_eq!(py_str_f64(f64::NEG_INFINITY), "-inf");
    }

    #[test]
    fn integer_semantics_match_python() {
        // python: -7 // 2 == -4 ; -7 % 2 == 1
        assert_eq!(py_floordiv_i64(-7, 2), -4);
        assert_eq!(py_mod_i64(-7, 2), 1);
        assert_eq!(py_floordiv_i64(7, 2), 3);
        assert_eq!(py_mod_i64(7, 2), 1);
        // python: math.ceil(0.05*3) == 1
        assert_eq!(py_ceil(0.05 * 3.0), 1);
        assert_eq!(py_ceil(2.0), 2);
        // negative divisors: the result takes the sign of the divisor
        // python: 7 // -2 == -4 ; 7 % -2 == -1 ; -7 // -2 == 3 ; -7 % -2 == -1
        assert_eq!(py_floordiv_i64(7, -2), -4);
        assert_eq!(py_mod_i64(7, -2), -1);
        assert_eq!(py_floordiv_i64(-7, -2), 3);
        assert_eq!(py_mod_i64(-7, -2), -1);
    }

    #[test]
    fn pynum_prints_by_type_not_by_value() {
        // python: str(99999999) == '99999999' but str(99999999.0) == '99999999.0'
        assert_eq!(PyNum::Int(99999999).to_py_string(), "99999999");
        assert_eq!(PyNum::Float(99999999.0).to_py_string(), "99999999.0");
        assert_eq!(PyNum::Float(0.0).to_py_string(), "0.0");
        assert_eq!(PyNum::Float(0.98).to_py_string(), "0.98");
        // python: str(None) == 'None'; cd-hit is really invoked with `-s None` by default
        assert_eq!(PyNum::None.to_py_string(), "None");
        assert_eq!(PyNum::from_opt_f64(None).to_py_string(), "None");
        assert_eq!(PyNum::from_opt_f64(Some(0.5)).to_py_string(), "0.5");
    }

    #[test]
    fn repr_of_a_split_line() {
        // what `print("Problem reading input line: ", line)` shows for a 3-column line
        assert_eq!(
            py_repr_str_list(&["a.gff", "a.fna", "extra"]),
            "['a.gff', 'a.fna', 'extra']"
        );
        assert_eq!(py_repr_str_list(&[]), "[]");
        // repr switches to double quotes when the string holds an apostrophe
        assert_eq!(py_repr_str("it's"), "\"it's\"");
        assert_eq!(py_repr_str("a\tb"), "'a\\tb'");
    }

    #[test]
    fn bytes_repr_matches_python() {
        // python: repr(b'hi\n') == "b'hi\\n'"
        assert_eq!(py_repr_bytes(b"hi\n"), r"b'hi\n'");
        // python: repr(b'\t== CD-HIT version 4.8.1 ==') keeps the ASCII verbatim
        assert_eq!(
            py_repr_bytes(b"\tCD-HIT version 4.8.1"),
            r"b'\tCD-HIT version 4.8.1'"
        );
        // python: repr(b'\x00\xff') == "b'\\x00\\xff'"
        assert_eq!(py_repr_bytes(b"\x00\xff"), r"b'\x00\xff'");
    }
}
