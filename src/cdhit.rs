//! Translation of `panaroo/panaroo/cdhit.py`.
//!
//! cd-hit itself stays an external process (rule 3). The command strings are assembled
//! here byte-identically to the Python, because they are printed verbatim when
//! `--verbose` is on and because cd-hit's clustering depends on the exact flags.
//!
//! `cluster_nodes_cdhit`, `is_valid` and `align_dna_cdhit` are **not reachable** from
//! `__main__::main` in the current upstream tree (`align_dna_cdhit` is imported by
//! `find_missing.py` but never called). They are declared for completeness and scheduled
//! in phase 6.

use crate::support::graph::Graph;
use crate::support::pyfmt::{py_str_f64, PyNum};
use crate::support::sparse::CsrMatrix;
use std::collections::HashMap;

/// `cdhit.py::check_cdhit_version`
///
/// Regexes `CD-HIT version \d+\.\d+` out of the **repr of the CompletedProcess object**,
/// not out of stdout — see [`crate::support::proc::run_shell_capture_repr`]. Exits the
/// process with status 1 if cd-hit is not runnable.
pub fn check_cdhit_version(cdhit_exec: &str) -> f64 {
    let p = crate::support::proc::run_shell_capture_repr(&format!("{cdhit_exec} -h"));
    // re.search(r'CD-HIT version \d+\.\d+', p) -- note this matches only two components,
    // so "CD-HIT version 4.8.1" yields 4.8, not 4.8.1.
    let version = find_version(&p);
    match version {
        Some(v) => v,
        None => {
            eprintln!("Need cd-hit to be runnable through: {cdhit_exec}");
            std::process::exit(1);
        }
    }
}

/// `re.search(r'CD-HIT version \d+\.\d+', text)` then `float(match.split()[-1])`.
///
/// Hand-rolled rather than pulling in a regex crate for one pattern.
fn find_version(text: &str) -> Option<f64> {
    const PREFIX: &str = "CD-HIT version ";
    let mut from = 0;
    while let Some(i) = text[from..].find(PREFIX) {
        let start = from + i + PREFIX.len();
        let rest = &text[start..];
        let digits1: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        if !digits1.is_empty() && rest[digits1.len()..].starts_with('.') {
            let after = &rest[digits1.len() + 1..];
            let digits2: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
            if !digits2.is_empty() {
                return format!("{digits1}.{digits2}").parse().ok();
            }
        }
        from = start;
    }
    None
}

/// Dispatch one assembled cd-hit command.
///
/// Default: hand it to the shell, exactly as `subprocess.run(cmd, shell=True, check=True)`
/// does. With the `cdhit-embedded` feature: run the translated cd-hit in-process instead,
/// which removes a `fork`/`exec` per call (12 per run) but leaves every file it reads and
/// writes exactly where it was.
///
/// `est` selects the nucleotide front-end, which differs from the protein one in its
/// pre-parse defaults and alphabet -- see [`crate::cdhit_embedded`].
fn dispatch_cdhit(cmd: &str, #[allow(unused_variables)] est: bool) {
    #[cfg(feature = "cdhit-embedded")]
    {
        let argv = crate::cdhit_embedded::argv_from_command(cmd);
        if est {
            crate::cdhit_embedded::run_cd_hit_est_main(&argv);
        } else {
            crate::cdhit_embedded::run_cd_hit_main(&argv);
        }
    }
    #[cfg(not(feature = "cdhit-embedded"))]
    {
        crate::support::proc::run_shell_check(cmd);
    }
}

/// `cdhit.py::run_cdhit`
///
/// Emits the flags in exactly the Python's order. Numeric flags go through Python's
/// `str()`, which is why `s`/`aL`/`aS` are [`PyNum`] rather than `f64` — see its docs.
///
/// OBSERVED, by monkeypatching `subprocess` in the reference Python:
///
/// ```text
/// run_cdhit('IN','OUT',id=0.98,s=0.98,quiet=True,n_cpu=4)
///   cd-hit -T 4 -i IN -o OUT -c 0.98 -s 0.98 -aL 0.0 -AL 99999999 -aS 0.0 -AS 99999999 -M 0 -d 999 -g 1 -n 2 > /dev/null
///
/// run_cdhit('IN','OUT',id=0.99,s=0.0,aS=99999999,quiet=True,n_cpu=4)
///   cd-hit -T 4 -i IN -o OUT -c 0.99 -s 0.0 -aL 0.0 -AL 99999999 -aS 99999999 -AS 99999999 -M 0 -d 999 -g 1 -n 2 > /dev/null
/// ```
///
/// Note the second: `-aS 99999999`, not `99999999.0`. That is the `aS=AS` call shape from
/// `iterative_cdhit`.
#[allow(clippy::too_many_arguments)]
pub fn run_cdhit(
    input_file: &str,
    output_file: &str,
    id: f64,
    n_cpu: i64,
    s: PyNum,
    a_l: PyNum,
    big_a_l: i64,
    a_s: PyNum,
    big_a_s: i64,
    accurate: bool,
    use_local: bool,
    word_length: Option<i64>,
    min_length: Option<i64>,
    quiet: bool,
) {
    let mut cmd = String::from("cd-hit");
    cmd += &format!(" -T {n_cpu}");
    cmd += &format!(" -i {input_file}");
    cmd += &format!(" -o {output_file}");
    cmd += &format!(" -c {}", py_str_f64(id));
    cmd += &format!(" -s {s}");
    cmd += &format!(" -aL {a_l}");
    cmd += &format!(" -AL {big_a_l}");
    cmd += &format!(" -aS {a_s}");
    cmd += &format!(" -AS {big_a_s}");
    cmd += " -M 0 -d 999";

    if use_local {
        cmd += " -G 0";
    }
    if accurate {
        cmd += " -g 1 -n 2";
    }
    if let Some(w) = word_length {
        if !accurate {
            cmd += &format!(" -n {w}");
        }
    }
    if let Some(l) = min_length {
        cmd += &format!(" -l {l}");
    }

    if !quiet {
        println!("running cmd: {cmd}");
    } else {
        cmd += " > /dev/null";
    }

    dispatch_cdhit(&cmd, false);
}

/// `cdhit.py::run_cdhit_est`
///
/// OBSERVED:
/// ```text
/// run_cdhit_est('IN','OUT',id=0.99,s=0.0,aS=99999999,accurate=False,word_length=7,quiet=True,n_cpu=4)
///   cd-hit-est -T 4 -i IN -o OUT -c 0.99 -s 0.0 -aL 0.0 -AL 99999999 -aS 99999999 -AS 99999999 -r 1 -M 0 -d 999 -mask NX -n 7 > /dev/null
/// ```
#[allow(clippy::too_many_arguments)]
pub fn run_cdhit_est(
    input_file: &str,
    output_file: &str,
    id: f64,
    n_cpu: i64,
    s: PyNum,
    a_l: PyNum,
    big_a_l: i64,
    a_s: PyNum,
    big_a_s: i64,
    accurate: bool,
    use_local: bool,
    strand: i64,
    print_aln: bool,
    word_length: Option<i64>,
    mask: bool,
    quiet: bool,
) {
    let mut cmd = String::from("cd-hit-est");
    cmd += &format!(" -T {n_cpu}");
    cmd += &format!(" -i {input_file}");
    cmd += &format!(" -o {output_file}");
    cmd += &format!(" -c {}", py_str_f64(id));
    cmd += &format!(" -s {s}");
    cmd += &format!(" -aL {a_l}");
    cmd += &format!(" -AL {big_a_l}");
    cmd += &format!(" -aS {a_s}");
    cmd += &format!(" -AS {big_a_s}");
    cmd += &format!(" -r {strand}");
    cmd += " -M 0 -d 999";

    if mask {
        cmd += " -mask NX";
    }
    if use_local {
        cmd += " -G 0";
    }
    if accurate {
        cmd += " -g 1 -n 6";
    }
    if let Some(w) = word_length {
        if !accurate {
            cmd += &format!(" -n {w}");
        }
    }
    if print_aln {
        cmd += " -p 1";
    }

    if !quiet {
        println!("running cmd: {cmd}");
    } else {
        cmd += " > /dev/null";
    }

    dispatch_cdhit(&cmd, true);
}

/// `cdhit.py::iterative_cdhit`
///
/// Runs cd-hit at each threshold in turn, feeding the previous round's representative
/// FASTA into the next. Returns the accumulated clusters as centroid-ID lists.
///
/// The Python plays a trick with `temp_input_file.name = temp_output_file.name` — it
/// rebinds the *attribute* on the NamedTemporaryFile object so the next round reads the
/// previous round's output. Reproduce the resulting file-path sequence, not the object
/// mutation.
#[allow(clippy::too_many_arguments)]
pub fn iterative_cdhit(
    g: &Graph,
    outdir: &str,
    dna: bool,
    s: PyNum,
    a_l: PyNum,
    big_a_l: i64,
    _a_s: PyNum,
    big_a_s: i64,
    accurate: bool,
    use_local: bool,
    strand: i64,
    quiet: bool,
    word_length: Option<i64>,
    thresholds: &[f64],
    n_cpu: i64,
) -> Vec<Vec<String>> {
    use crate::support::pydict::PyDict;

    let temp_input_file = temp_name(outdir);
    let temp_output_file = temp_name(outdir);

    let centroid_to_seq = centroid_to_seq(g, dna);

    let mut clusters: Vec<Vec<String>> = Vec::new();
    {
        let mut out = String::new();
        for centroid in centroid_to_seq.keys() {
            clusters.push(vec![centroid.clone()]);
            out.push('>');
            out.push_str(centroid);
            out.push('\n');
            out.push_str(centroid_to_seq.get(centroid).unwrap());
            out.push('\n');
        }
        std::fs::write(&temp_input_file, out).expect("write cd-hit input");
    }

    // The Python rebinds `temp_input_file.name` / `temp_output_file.name` each round so the
    // next round reads the previous round's output, and appends "t{cid}" to the output name.
    // Reproduce the resulting *path sequence*, not the object mutation.
    let mut in_path = temp_input_file.clone();
    let mut out_path = temp_output_file.clone();

    for &cid in thresholds {
        if dna {
            run_cdhit_est(
                &in_path,
                &out_path,
                cid,
                n_cpu,
                s,
                a_l,
                big_a_l,
                // UPSTREAM: `aS=AS` -- the Python passes the *control* parameter's int into
                // the coverage-fraction parameter (cdhit.py:417). Almost certainly a slip,
                // but it is what upstream does, and the type matters: `str(99999999)` is
                // "99999999" while `str(99999999.0)` would be "99999999.0", a different
                // cd-hit flag. Hence PyNum. Not in the Tier B set because the effect has
                // not been measured.
                PyNum::Int(big_a_s),
                big_a_s,
                accurate,
                use_local,
                strand,
                false,
                word_length,
                true,
                quiet,
            );
        } else {
            run_cdhit(
                // `aS=AS` again (cdhit.py:431); see the est branch above.
                &in_path,
                &out_path,
                cid,
                n_cpu,
                s,
                a_l,
                big_a_l,
                PyNum::Int(big_a_s),
                big_a_s,
                accurate,
                use_local,
                word_length,
                None,
                quiet,
            );
        }

        // process the output
        let clstr = format!("{out_path}.clstr");
        let temp_clusters = parse_clstr_strings(&clstr);

        // collapse previously clustered
        let mut temp_clust_dict: HashMap<String, i64> = HashMap::new();
        for (c, clust) in temp_clusters.iter().enumerate() {
            for n in clust {
                temp_clust_dict.insert(n.clone(), c as i64);
            }
        }
        let mut clust_dict: PyDict<i64, Vec<String>> = PyDict::new();
        for clust in &clusters {
            let mut c: i64 = -1;
            for n in clust {
                if let Some(&v) = temp_clust_dict.get(n) {
                    c = v;
                }
            }
            for n in clust {
                if !clust_dict.contains_key(&c) {
                    clust_dict.insert(c, Vec::new());
                }
                clust_dict.get_mut(&c).unwrap().push(n.clone());
            }
        }
        clusters = clust_dict.values().cloned().collect();

        // cleanup and rename for next round
        let _ = std::fs::remove_file(&in_path);
        let _ = std::fs::remove_file(&clstr);
        in_path = out_path.clone();
        out_path = format!("{out_path}t{}", crate::support::pyfmt::py_str_f64(cid));
    }

    clusters
}

/// Parse a cd-hit `.clstr` file into clusters of sequence names.
///
/// Not a Python function — the same seven-line loop appears in `iterative_cdhit` and
/// `cluster_nodes_cdhit`, differing only in whether the name is parsed as an `int`. This is
/// the string form. Kept here because it reads a *tool's* output format, not Panaroo logic.
///
/// The Python appends the accumulating cluster on every `>` line and then drops the first
/// (empty) entry, which is equivalent to splitting on `>Cluster` headers.
fn parse_clstr_strings(path: &str) -> Vec<Vec<String>> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("could not read {path}: {e}"));
    let mut clusters: Vec<Vec<String>> = Vec::new();
    let mut c: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            clusters.push(std::mem::take(&mut c));
        } else if let Some(rest) = line.split_once('>') {
            let name = rest.1.split("...").next().unwrap_or("").to_string();
            c.push(name);
        }
    }
    clusters.push(c);
    clusters.remove(0);
    clusters
}

/// `tempfile.NamedTemporaryFile(delete=False, dir=outdir)` — a unique path in `outdir`.
///
/// The exact name never reaches the output, so only uniqueness matters.
fn temp_name(outdir: &str) -> String {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static N: AtomicUsize = AtomicUsize::new(0);
    let i = N.fetch_add(1, Ordering::Relaxed);
    let pid = std::process::id();
    format!("{outdir}tmp{pid}_{i}")
}

/// `cdhit.py::pwdist_edlib`
///
/// All-pairs identity within each cd-hit cluster, thresholded into a boolean CSR matrix.
/// Returns `(distances_bwtn_centroids, centroid_to_index)`.
///
/// PORTING_PLAN.md §9 item 2: the Python creates a fresh joblib pool per cluster. Result
/// order is unaffected, so [`crate::support::parallel::parallel_map`] may use one pool —
/// see the note there.
pub fn pwdist_edlib(
    g: &Graph,
    cdhit_clusters: &[Vec<String>],
    threshold: f64,
    dna: bool,
    n_cpu: i64,
) -> (CsrMatrix, HashMap<String, usize>) {
    // Prepare sequences
    let centroid_to_seq = centroid_to_seq(g, dna);
    let ncentroids = centroid_to_seq.len();

    // centroid to index -- insertion order of centroid_to_seq, i.e. node order then the
    // order of each node's centroid list
    let mut centroid_to_index: HashMap<String, usize> = HashMap::new();
    for (i, centroid) in centroid_to_seq.keys().enumerate() {
        centroid_to_index.insert(centroid.clone(), i);
    }

    // get pairwise id between sequences in the same cdhit clusters
    //
    // PORTING_PLAN.md §9 item 2: the Python creates a fresh joblib pool per cluster. Result
    // order is unaffected, so one pool over the flattened pair list is observationally
    // identical -- see support::parallel.
    let mut pairs: Vec<(String, String)> = Vec::with_capacity(
        cdhit_clusters
            .iter()
            .map(|c| c.len() * c.len().saturating_sub(1) / 2)
            .sum(),
    );
    for cluster in cdhit_clusters {
        for i in 0..cluster.len() {
            for j in (i + 1)..cluster.len() {
                pairs.push((cluster[i].clone(), cluster[j].clone()));
            }
        }
    }

    let all_distances = crate::support::parallel::parallel_map(n_cpu, pairs, |(c1, c2)| {
        run_pw_thresholded(
            centroid_to_seq
                .get(&c1)
                .unwrap_or_else(|| panic!("KeyError: {c1}")),
            centroid_to_seq
                .get(&c2)
                .unwrap_or_else(|| panic!("KeyError: {c2}")),
            centroid_to_index[&c1],
            centroid_to_index[&c2],
            dna,
            threshold,
        )
    });

    let mut data = Vec::new();
    let mut row_ind = Vec::new();
    let mut col_ind = Vec::new();
    for d in all_distances {
        if d.2 >= threshold {
            data.push(1);
            row_ind.push(d.0);
            col_ind.push(d.1);
        }
    }

    let distances_bwtn_centroids =
        CsrMatrix::from_coo(data, row_ind, col_ind, (ncentroids, ncentroids));

    (distances_bwtn_centroids, centroid_to_index)
}

/// The `centroid_to_seq` dict built identically at the top of `pwdist_edlib` and
/// `iterative_cdhit`.
///
/// Not a Python function — the Python inlines this loop in both places. Rule 2 says keep
/// duplicated logic duplicated, but this one is *support* for two translated functions
/// rather than a translation itself, and the insertion order it produces is load-bearing
/// (it fixes `centroid_to_index`, hence every matrix index). Factoring it guarantees the
/// two call sites cannot drift apart.
fn centroid_to_seq(g: &Graph, dna: bool) -> crate::support::pydict::PyDict<String, String> {
    let mut m = crate::support::pydict::PyDict::new();
    for node in g.nodes() {
        let n = g.node(node);
        let seqs = if dna { &n.dna } else { &n.protein };
        for (sid, seq) in n.centroid.iter().zip(seqs.iter()) {
            m.insert(sid.clone(), seq.clone());
        }
    }
    m
}

/// `cdhit.py::run_pw`
///
/// Pairwise identity via edlib in `HW` mode with `k = 0.5 * len(seqA)`. Swaps so `seqA` is
/// the shorter. For DNA it takes the max over the forward and reverse-complement strands.
///
/// The DNA branch assigns `pqid = max(pwid, 0.0)` (note the typo) when `editDistance == -1`.
/// This is a real typo but is **provably behaviour-neutral**: `pwid` starts at `0.0` and is
/// only ever assigned `max(pwid, ...)`, so the intended `pwid = max(pwid, 0.0)` is a no-op
/// too. Transcribe the branch as a no-op and do not "fix" it. See ORIGINAL_CODE_BUG.md B3.
pub fn run_pw(seq_a: &str, seq_b: &str, n1: usize, n2: usize, dna: bool) -> (usize, usize, f64) {
    // NEG_INFINITY disables the extra bound, so this is the Python's `k` exactly.
    run_pw_thresholded(seq_a, seq_b, n1, n2, dna, f64::NEG_INFINITY)
}

/// [`run_pw`] with the caller's identity threshold used as an additional edlib `k` bound.
///
/// The only consumer of `run_pw`'s float is `pwdist_edlib`'s `d.2 >= threshold`. An
/// alignment rejected by the tighter bound has `editDistance > (1 - threshold) * len(seqA)`,
/// hence `1 - ed/len < threshold`, so it could only ever have contributed a value BELOW the
/// threshold -- and dropping a value below `T` from a `max` cannot change `max >= T`. When
/// every strand is rejected the result is `0.0`, which is `< T` for any `T > 0`. So the
/// boolean the caller computes is unchanged, while edlib's Ukkonen band
/// (`ceil((k+1)/64)` blocks) shrinks by up to 8x on the DNA pass.
///
/// The returned float itself is NOT the same for rejected pairs, which is why `run_pw`
/// keeps the unbounded behaviour for any other caller.
fn run_pw_thresholded(
    seq_a: &str,
    seq_b: &str,
    n1: usize,
    n2: usize,
    dna: bool,
    threshold: f64,
) -> (usize, usize, f64) {
    use crate::support::edlib::{align, Mode, Task};

    let (seq_a, seq_b) = if seq_a.len() > seq_b.len() {
        (seq_b, seq_a)
    } else {
        (seq_a, seq_b)
    };

    // The Python's bound, unchanged.
    let k_python = (0.5 * seq_a.len() as f64) as i64;
    let k = if threshold.is_finite() && threshold > 0.0 {
        let k_thresh = (((1.0 - threshold) * seq_a.len() as f64).ceil() as i64 + 1).max(0);
        k_python.min(k_thresh)
    } else {
        k_python
    };

    let pwid = if dna {
        let mut acc = 0.0f64;
        let rc = crate::support::seq::reverse_complement(seq_a);
        for s_a in [seq_a, rc.as_str()] {
            let aln = align(s_a, seq_b, Mode::Hw, Task::Distance, k, &DNA_N_EQUALITIES);
            if aln.edit_distance == -1 {
                // UPSTREAM: the Python writes `pqid = max(pwid, 0.0)` here -- a typo, since
                // `pqid` is never read. It is provably behaviour-neutral: `pwid` starts at
                // 0.0 and only ever takes `max(pwid, ...)`, so the intended
                // `pwid = max(pwid, 0.0)` is a no-op too. Transcribed as a no-op; do not
                // "fix" it. See ORIGINAL_CODE_BUG.md B3.
            } else {
                acc = acc.max(1.0 - aln.edit_distance as f64 / seq_a.len() as f64);
            }
        }
        acc
    } else {
        let aln = align(
            seq_a,
            seq_b,
            Mode::Hw,
            Task::Distance,
            k,
            &PROTEIN_X_EQUALITIES,
        );
        if aln.edit_distance == -1 {
            0.0
        } else {
            1.0 - aln.edit_distance as f64 / seq_a.len() as f64
        }
    };

    (n1, n2, pwid)
}

/// `additionalEqualities` for the DNA branch of `run_pw`: `N` matches any base.
const DNA_N_EQUALITIES: [(char, char); 4] = [('A', 'N'), ('C', 'N'), ('G', 'N'), ('T', 'N')];

/// `additionalEqualities` for the protein branch of `run_pw`: `X` matches any residue, plus
/// the two IUPAC ambiguity groups `B` (D/N) and `Z` (E/Q). Order is copied verbatim from
/// the Python — edlib treats the list as a set, so order does not matter, but keeping it
/// makes the two readable side by side.
const PROTEIN_X_EQUALITIES: [(char, char); 28] = [
    ('*', 'X'),
    ('A', 'X'),
    ('C', 'X'),
    ('B', 'X'),
    ('E', 'X'),
    ('D', 'X'),
    ('G', 'X'),
    ('F', 'X'),
    ('I', 'X'),
    ('H', 'X'),
    ('K', 'X'),
    ('M', 'X'),
    ('L', 'X'),
    ('N', 'X'),
    ('Q', 'X'),
    ('P', 'X'),
    ('S', 'X'),
    ('R', 'X'),
    ('T', 'X'),
    ('W', 'X'),
    ('V', 'X'),
    ('Y', 'X'),
    ('X', 'X'),
    ('Z', 'X'),
    ('D', 'B'),
    ('N', 'B'),
    ('E', 'Z'),
    ('Q', 'Z'),
];

// --- not reachable from main; phase 6 ---------------------------------------------------

/// `cdhit.py::cluster_nodes_cdhit` — unreachable from `main`. Phase 6.
#[allow(clippy::too_many_arguments)]
pub fn cluster_nodes_cdhit(
    _g: &Graph,
    _nodes: &[usize],
    _outdir: &str,
    _id: f64,
    _dna: bool,
    _s: PyNum,
    _a_l: PyNum,
    _big_a_l: i64,
    _a_s: PyNum,
    _big_a_s: i64,
    _accurate: bool,
    _use_local: bool,
    _strand: i64,
    _quiet: bool,
    _prevent_para: bool,
    _n_cpu: i64,
) -> Vec<Vec<usize>> {
    panic!("noimpl: cdhit::cluster_nodes_cdhit")
}

/// `cdhit.py::is_valid` — unreachable from `main`. Phase 6.
pub fn is_valid(_g: &Graph, _node: usize, _cluster: &[usize]) -> bool {
    panic!("noimpl: cdhit::is_valid")
}

/// `cdhit.py::align_dna_cdhit` — imported by `find_missing.py` but never called. Phase 6.
#[allow(clippy::too_many_arguments)]
pub fn align_dna_cdhit(
    _query: &str,
    _target: &str,
    _temp_dir: &str,
    _id: f64,
    _n_cpu: i64,
    _s: PyNum,
    _a_l: PyNum,
    _big_a_l: i64,
    _a_s: PyNum,
    _big_a_s: i64,
    _accurate: bool,
    _use_local: bool,
    _strand: i64,
    _mask: bool,
    _quiet: bool,
) -> String {
    panic!("noimpl: cdhit::align_dna_cdhit")
}
