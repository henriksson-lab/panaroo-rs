//! Translation of `panaroo/panaroo/generate_alignments.py`.
//!
//! The largest module in scope: 36 functions covering resume-manifest bookkeeping, aligner
//! invocation, and the codon-alignment path. External aligners stay external (rule 3).

use crate::prokka::TransTable;
use crate::support::align_io::MultipleSeqAlignment;
use crate::support::graph::NodeAttrs;
use crate::support::seqio::SeqRecord;

/// `generate_alignments.py::get_trans_table`
///
/// Byte-identical to [`crate::prokka::get_trans_table`] upstream. Both are kept; see
/// PORTING_PLAN.md §5.
pub fn get_trans_table(table: i64) -> TransTable {
    // Byte-identical to prokka's; rule 2 keeps the duplication rather than delegating,
    // but the body is the same code.
    crate::prokka::get_trans_table(table)
}

/// `generate_alignments.py::translate`
///
/// Byte-identical to [`crate::prokka::translate`] upstream. Both are kept.
pub fn translate(seq: &str, translation_table: &TransTable) -> String {
    crate::prokka::translate(seq, translation_table)
}

// --- path construction ------------------------------------------------------------------
//
// The truncation limits (237/236 and 248) exist to stay under filesystem name limits. They
// are load-bearing: two genes whose names share a 236-character prefix collide, and the
// pipeline depends on that collision being resolved the same way it is in Python.

/// `generate_alignments.py::get_alignment_basename`
///
/// ```python
/// gene_name = node["name"]
/// if len(gene_name) >= 237:
///     return gene_name[:236]
/// return gene_name
/// ```
///
/// The 237/236 limit keeps the derived filenames under the ext4 255-byte cap. It is
/// load-bearing, not cosmetic: two genes whose names share a 236-character prefix collide
/// on disk, and the pipeline depends on that collision resolving the same way it does in
/// Python.
pub fn get_alignment_basename(node: &NodeAttrs) -> String {
    let gene_name = node.name.as_deref().expect("node['name'] not set");
    if gene_name.len() >= 237 {
        return gene_name[..236].to_string();
    }
    gene_name.to_string()
}

/// `generate_alignments.py::get_temp_dna_input_path`
///
/// ```python
/// outname = temp_directory + node["name"] + ".fasta"
/// if len(outname) >= 248:
///     outname = outname[:248] + ".fasta"
/// return outname
/// ```
///
/// Note the truncation appends `.fasta` to an already-248-character string, so the result
/// is 254 characters and ends `.fasta.fasta` minus the overlap. Verified against Python:
/// a 300-character name under `/t/` yields a 254-character path. Preserve it.
///
/// This is string concatenation, not `os.path.join`, so `temp_directory` must already end
/// in a separator — which is how `main` passes it.
pub fn get_temp_dna_input_path(node: &NodeAttrs, temp_directory: &str) -> String {
    let name = node.name.as_deref().expect("node['name'] not set");
    let outname = format!("{temp_directory}{name}.fasta");
    if outname.len() >= 248 {
        return format!("{}.fasta", &outname[..248]);
    }
    outname
}

/// `generate_alignments.py::get_expected_gene_alignment_path`
///
/// With `codons`, always the `.aln.fas` path. Without, a gene with a single sequence gets
/// `{node['name']}.fasta` — note that branch uses the **untruncated** name, unlike every
/// other path here.
pub fn get_expected_gene_alignment_path(
    node: &NodeAttrs,
    output_dir: &str,
    codons: bool,
) -> String {
    if codons {
        return join(
            output_dir,
            &[
                "aligned_gene_sequences",
                &format!("{}.aln.fas", get_alignment_basename(node)),
            ],
        );
    }

    if node.seq_ids.len() > 1 {
        return join(
            output_dir,
            &[
                "aligned_gene_sequences",
                &format!("{}.aln.fas", get_alignment_basename(node)),
            ],
        );
    }

    let name = node.name.as_deref().expect("node['name'] not set");
    join(
        output_dir,
        &["aligned_gene_sequences", &format!("{name}.fasta")],
    )
}

/// `generate_alignments.py::get_expected_protein_input_path`
pub fn get_expected_protein_input_path(node: &NodeAttrs, temp_directory: &str) -> String {
    join(
        temp_directory,
        &[&format!("{}.fasta", get_alignment_basename(node))],
    )
}

/// `generate_alignments.py::get_expected_protein_alignment_path`
pub fn get_expected_protein_alignment_path(node: &NodeAttrs, output_dir: &str) -> String {
    join(
        output_dir,
        &[
            "aligned_protein_sequences",
            &format!("{}.aln.fas", get_alignment_basename(node)),
        ],
    )
}

/// `generate_alignments.py::get_expected_unaligned_dna_path`
pub fn get_expected_unaligned_dna_path(node: &NodeAttrs, output_dir: &str) -> String {
    join(
        output_dir,
        &[
            "unaligned_dna_sequences",
            &format!("{}.fasta", get_alignment_basename(node)),
        ],
    )
}

/// `generate_alignments.py::get_resume_manifest_path`
pub fn get_resume_manifest_path(output_dir: &str) -> String {
    join(output_dir, &["alignment_resume_state.json"])
}

/// `os.path.join` for POSIX, which is what the paths above use.
///
/// Not a Python function — `os.path.join` has no direct Rust equivalent that matches its
/// behaviour on already-separator-terminated components, and `PathBuf::join` normalises in
/// ways that would change the emitted strings. Rules reproduced: an absolute later
/// component discards everything before it, and a separator is inserted only when the
/// accumulated path does not already end in one.
fn join(base: &str, parts: &[&str]) -> String {
    let mut out = base.to_string();
    for p in parts {
        if p.starts_with('/') {
            out = p.to_string();
        } else if out.is_empty() || out.ends_with('/') {
            out.push_str(p);
        } else {
            out.push('/');
            out.push_str(p);
        }
    }
    out
}

// --- resume manifest --------------------------------------------------------------------

/// The JSON object written by [`write_resume_manifest`].
///
/// Serialised with `json.dump(..., indent=2, sort_keys=True)` plus a trailing newline, so
/// keys are emitted in **alphabetical** order: `aligner, alignment, codons, core_threshold,
/// started_at, strict_codons, subset`.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ResumeManifest {
    /// `datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")`
    pub started_at: String,
    /// `"core"` or `"pan"`
    pub alignment: String,
    pub aligner: String,
    pub codons: bool,
    pub strict_codons: bool,
    pub core_threshold: f64,
    pub subset: Option<i64>,
}

/// `generate_alignments.py::load_resume_manifest` — `None` if absent.
pub fn load_resume_manifest(output_dir: &str) -> Option<ResumeManifest> {
    let manifest_path = get_resume_manifest_path(output_dir);
    if !std::path::Path::new(&manifest_path).is_file() {
        return None;
    }
    let text = std::fs::read_to_string(&manifest_path).ok()?;
    serde_json::from_str(&text).ok()
}

/// `generate_alignments.py::check_resume_manifest_collision`
pub fn check_resume_manifest_collision(output_dir: &str, resume: bool) {
    if resume {
        return;
    }
    let manifest_path = get_resume_manifest_path(output_dir);
    if !std::path::Path::new(&manifest_path).is_file() {
        return;
    }
    panic!(
        "RuntimeError: Found an existing gene-alignment resume manifest in {output_dir}. \
         Re-run with --resume to continue the previous alignment. To start the sequence \
         alignment again from scratch, delete: {manifest_path}, {}, {}, and {}.",
        join(output_dir, &["aligned_gene_sequences/"]),
        join(output_dir, &["aligned_protein_sequences/"]),
        join(output_dir, &["unaligned_dna_sequences/"]),
    );
}

/// `generate_alignments.py::write_resume_manifest`
#[allow(clippy::too_many_arguments)]
pub fn write_resume_manifest(
    output_dir: &str,
    alignment: &str,
    aligner: &str,
    codons: bool,
    strict_codons: bool,
    core_threshold: f64,
    subset: Option<i64>,
    resume: bool,
) -> ResumeManifest {
    let manifest_path = get_resume_manifest_path(output_dir);
    let existing = load_resume_manifest(output_dir);

    if resume {
        let Some(m) = existing else {
            panic!(
                "RuntimeError: Cannot resume panaroo-msa: no alignment resume manifest was found."
            );
        };
        for (field, matches) in [
            ("alignment", m.alignment == alignment),
            ("aligner", m.aligner == aligner),
            ("codons", m.codons == codons),
            ("strict_codons", m.strict_codons == strict_codons),
            ("core_threshold", m.core_threshold == core_threshold),
            ("subset", m.subset == subset),
        ] {
            if !matches {
                panic!(
                    "RuntimeError: Cannot resume panaroo-msa: current run does not match \
                     the existing manifest for '{field}'."
                );
            }
        }
        return m;
    }

    let manifest = ResumeManifest {
        started_at: utc_now_iso(),
        alignment: alignment.to_string(),
        aligner: aligner.to_string(),
        codons,
        strict_codons,
        core_threshold,
        subset,
    };

    // `json.dump(manifest, handle, indent=2, sort_keys=True)` then a trailing newline.
    // sort_keys means alphabetical: aligner, alignment, codons, core_threshold,
    // started_at, strict_codons, subset.
    let json = format!(
        "{{\n  \"aligner\": {},\n  \"alignment\": {},\n  \"codons\": {},\n  \"core_threshold\": {},\n  \"started_at\": {},\n  \"strict_codons\": {},\n  \"subset\": {}\n}}\n",
        serde_json::to_string(&manifest.aligner).unwrap(),
        serde_json::to_string(&manifest.alignment).unwrap(),
        manifest.codons,
        crate::support::pyfmt::py_str_f64(manifest.core_threshold),
        serde_json::to_string(&manifest.started_at).unwrap(),
        manifest.strict_codons,
        match manifest.subset { Some(n) => n.to_string(), None => "null".to_string() },
    );
    std::fs::write(&manifest_path, json).expect("write resume manifest");
    manifest
}

/// `datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")`
///
/// Not a Python function. Deliberately excluded from parity comparison — a timestamp
/// cannot match between two runs, so `tests/parity/canonicalise.py` must skip
/// `alignment_resume_state.json` (or its `started_at` field) when comparing.
fn utc_now_iso() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system clock before 1970")
        .as_secs() as i64;
    // civil-from-days, Howard Hinnant's algorithm
    let days = secs.div_euclid(86400);
    let rem = secs.rem_euclid(86400);
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z.rem_euclid(146097);
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60
    )
}

// --- resume predicates ------------------------------------------------------------------

/// `generate_alignments.py::is_valid_fasta`
pub fn is_valid_fasta(path: &str) -> bool {
    if !std::path::Path::new(path).is_file() {
        return false;
    }
    match std::fs::read_to_string(path) {
        Ok(text) => !crate::support::seqio::parse_fasta(&text).is_empty(),
        Err(_) => false,
    }
}

/// `generate_alignments.py::is_valid_alignment`
pub fn is_valid_alignment(path: &str) -> bool {
    if !std::path::Path::new(path).is_file() {
        return false;
    }
    // AlignIO.read raises unless every record has the same length; that check is what
    // distinguishes this from is_valid_fasta.
    match std::fs::read_to_string(path) {
        Ok(text) => {
            let recs = crate::support::seqio::parse_fasta(&text);
            !recs.is_empty() && recs.iter().all(|r| r.seq.len() == recs[0].seq.len())
        }
        Err(_) => false,
    }
}

/// `generate_alignments.py::gene_has_valid_final_output`
pub fn gene_has_valid_final_output(node: &NodeAttrs, output_dir: &str, codons: bool) -> bool {
    let output_path = get_expected_gene_alignment_path(node, output_dir, codons);
    if output_path.ends_with(".aln.fas") {
        is_valid_alignment(&output_path)
    } else {
        is_valid_fasta(&output_path)
    }
}

/// `generate_alignments.py::gene_has_valid_protein_output`
pub fn gene_has_valid_protein_output(node: &NodeAttrs, output_dir: &str) -> bool {
    is_valid_alignment(&get_expected_protein_alignment_path(node, output_dir))
}

/// `generate_alignments.py::node_requires_msa`
///
/// ```python
/// sequence_ids = node["seqIDs"]
/// if isinstance(sequence_ids, str):
///     sequence_ids = [sequence_ids]
/// return len(sequence_ids) > 1
/// ```
///
/// The `isinstance(..., str)` guard is for graphs read back from GML, where a single-member
/// node's `seqIDs` round-trips as a bare string; `len()` on it would count characters. Our
/// `seq_ids` is always a set, so the guard is structurally impossible here.
pub fn node_requires_msa(node: &NodeAttrs) -> bool {
    node.seq_ids.len() > 1
}

/// `generate_alignments.py::get_pending_gene_ids`
pub fn get_pending_gene_ids(
    nodes: &[(usize, &NodeAttrs)],
    output_dir: &str,
    codons: bool,
    resume: bool,
) -> Vec<usize> {
    let mut pending = Vec::new();
    for (node_id, node) in nodes {
        if resume && gene_has_valid_final_output(node, output_dir, codons) {
            continue;
        }
        pending.push(*node_id);
    }
    pending
}

/// `generate_alignments.py::get_pending_codon_gene_ids`
///
/// Returns `(protein_pending_gene_ids, reverse_translate_pending_gene_ids)`.
pub fn get_pending_codon_gene_ids(
    nodes: &[(usize, &NodeAttrs)],
    output_dir: &str,
    resume: bool,
) -> (Vec<usize>, Vec<usize>) {
    let mut protein_pending = Vec::new();
    let mut reverse_translate_pending = Vec::new();

    for (node_id, node) in nodes {
        if resume && gene_has_valid_final_output(node, output_dir, true) {
            continue;
        }
        if !node_requires_msa(node) {
            protein_pending.push(*node_id);
            continue;
        }
        reverse_translate_pending.push(*node_id);
        if resume && gene_has_valid_protein_output(node, output_dir) {
            continue;
        }
        protein_pending.push(*node_id);
    }

    (protein_pending, reverse_translate_pending)
}

/// `generate_alignments.py::get_codon_pending_files`
///
/// Returns `(protein_alignment_files, dna_sequence_files)`.
pub fn get_codon_pending_files(
    nodes: &[(usize, &NodeAttrs)],
    output_dir: &str,
    gene_ids: &[usize],
) -> (Vec<String>, Vec<String>) {
    let gene_id_set: std::collections::HashSet<usize> = gene_ids.iter().copied().collect();
    let selected: Vec<&NodeAttrs> = nodes
        .iter()
        .filter(|(id, _)| gene_id_set.contains(id))
        .map(|(_, n)| *n)
        .collect();

    let protein_alignment_files = selected
        .iter()
        .map(|n| get_expected_protein_alignment_path(n, output_dir))
        .collect();
    let dna_sequence_files = selected
        .iter()
        .map(|n| get_expected_unaligned_dna_path(n, output_dir))
        .collect();

    (protein_alignment_files, dna_sequence_files)
}

/// `generate_alignments.py::print_stage_progress`
pub fn print_stage_progress(
    stage_name: &str,
    completed: usize,
    remaining: usize,
    total: Option<usize>,
) {
    let total = total.unwrap_or(completed + remaining);
    println!(
        "{stage_name}: {completed} completed alignments found, \
         {remaining} to be aligned out of {total}."
    );
}

// --- aligner discovery ------------------------------------------------------------------

/// `generate_alignments.py::check_aligner_install`
///
/// Like `check_cdhit_version`, this greps the **repr of the CompletedProcess object**.
pub fn check_aligner_install(aligner: &str) -> bool {
    let command = match aligner {
        "clustal" => "clustalo --help",
        "prank" => "prank -help",
        "mafft" => "mafft --help",
        "muscle" | "muscle-super5" => "muscle -h",
        "famsa" => "famsa -h",
        "none" => return true,
        _ => {
            eprintln!("Incorrect aligner specification");
            std::process::exit(0); // bare sys.exit() -- status 0, matching upstream
        }
    };

    // As with check_cdhit_version, this regexes the *repr of the CompletedProcess object*.
    // Note this call captures stderr too, so a banner printed there is visible -- unlike
    // check_cdhit_version, which captures only stdout.
    let p = crate::support::proc::run_shell_capture_repr_both(command);

    let present = match aligner {
        "clustal" => find_pat(
            &p,
            "Clustal Omega - ",
            &[
                Pat::Digits,
                Pat::Lit("."),
                Pat::Digits,
                Pat::Lit("."),
                Pat::Digits,
            ],
        ),
        "prank" => find_pat(&p, "prank v.", &[Pat::Digits, Pat::Lit(".")]),
        "mafft" => find_pat(&p, "MAFFT v", &[Pat::Digits, Pat::Lit("."), Pat::Digits]),
        // case-insensitive `muscle\s+\d+\.\d+\.\S+`
        "muscle" | "muscle-super5" => find_muscle(&p),
        "famsa" => find_famsa(&p),
        _ => false,
    };

    if !present {
        eprintln!("Need specified aligner to be installed ");
        std::process::exit(1);
    }
    present
}

/// A tiny pattern element, so the version probes above can be written without a regex crate.
enum Pat {
    Digits,
    Lit(&'static str),
}

fn find_pat(text: &str, prefix: &str, rest: &[Pat]) -> bool {
    let mut from = 0;
    while let Some(i) = text[from..].find(prefix) {
        let start = from + i + prefix.len();
        let mut pos = start;
        let mut ok = true;
        for p in rest {
            match p {
                Pat::Digits => {
                    let n = text[pos..]
                        .chars()
                        .take_while(|c| c.is_ascii_digit())
                        .count();
                    if n == 0 {
                        ok = false;
                        break;
                    }
                    pos += n;
                }
                Pat::Lit(l) => {
                    if text[pos..].starts_with(l) {
                        pos += l.len();
                    } else {
                        ok = false;
                        break;
                    }
                }
            }
        }
        if ok {
            return true;
        }
        from = start;
    }
    false
}

/// `re.search(r"muscle\s+\d+\.\d+\.\S+", p, re.IGNORECASE)`
fn find_muscle(text: &str) -> bool {
    let lower = text.to_lowercase();
    let mut from = 0;
    while let Some(i) = lower[from..].find("muscle") {
        let mut pos = from + i + "muscle".len();
        let ws = lower[pos..]
            .chars()
            .take_while(|c| c.is_whitespace())
            .count();
        if ws > 0 {
            pos += ws;
            let a = lower[pos..]
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .count();
            if a > 0 && lower[pos + a..].starts_with('.') {
                let mut q = pos + a + 1;
                let b = lower[q..]
                    .chars()
                    .take_while(|c| c.is_ascii_digit())
                    .count();
                if b > 0 && lower[q + b..].starts_with('.') {
                    q += b + 1;
                    if lower[q..].starts_with(|c: char| !c.is_whitespace()) {
                        return true;
                    }
                }
            }
        }
        from = from + i + "muscle".len();
    }
    false
}

/// `re.search(r"FAMSA.*?version\s+\d+\.\d+\.\d+(?:-[A-Za-z0-9]+)?", p,
/// re.IGNORECASE | re.DOTALL)`
fn find_famsa(text: &str) -> bool {
    let lower = text.to_lowercase();
    match lower.find("famsa") {
        Some(i) => {
            find_pat(
                &lower[i..],
                "version",
                &[
                    Pat::Digits,
                    Pat::Lit("."),
                    Pat::Digits,
                    Pat::Lit("."),
                    Pat::Digits,
                ],
            ) || {
                // `version\s+` -- allow the whitespace the helper does not model
                let rest = &lower[i..];
                let mut from = 0;
                let mut hit = false;
                while let Some(j) = rest[from..].find("version") {
                    let mut pos = from + j + "version".len();
                    pos += rest[pos..]
                        .chars()
                        .take_while(|c| c.is_whitespace())
                        .count();
                    let a = rest[pos..]
                        .chars()
                        .take_while(|c| c.is_ascii_digit())
                        .count();
                    if a > 0 {
                        hit = true;
                        break;
                    }
                    from = from + j + "version".len();
                }
                hit
            }
        }
        None => false,
    }
}

/// `generate_alignments.py::check_aligner_sanity`
///
/// Warns (does not fail) when muscle is asked to align a large dataset.
pub fn check_aligner_sanity(aligner: &str, codons: bool, isolate_count: usize) -> bool {
    if aligner == "famsa" && !codons {
        panic!(
            "RuntimeError: FAMSA2 only supports amino-acid alignment.\
             Use --codons or --strict-codons to align."
        );
    } else if aligner == "muscle" && isolate_count > 300 {
        eprintln!(
            "UserWarning: MUSCLE is not optimised to run on more than a few\
             hundred isolates. Aligning the core genome may be very\
             slow or fail to complete. Use muscle-super5 for faster\
             alignment on larger datasets"
        );
    }
    true
}

// --- writing per-gene input files -------------------------------------------------------

/// `combined_DNA_CDS.fasta`, parsed once for a whole alignment stage.
///
/// **Not a Python object.** It exists to retire PORTING_PLAN.md §9 item 1: the Python
/// `output_sequence` calls `SeqIO.parse(outdir + "combined_DNA_CDS.fasta")` *inside* itself,
/// so the entire file is re-read and re-parsed once per gene cluster — O(n_genes × filesize),
/// with every joblib worker hammering the same file. On the `ci` parity dataset that is
/// 5,116 clusters × 19,507,986 bytes ≈ **99.8 GB parsed** and ~1.1 × 10⁹ heap allocations
/// per run, and it is quadratic in genome count. This type hoists the parse to the caller so
/// it happens exactly once.
///
/// The caller-side hoist is the one deviation from PORTING_PLAN.md rule 2 in this module,
/// and it is deliberate: the loop body itself is transcribed unchanged, only the *number of
/// times the file is read* differs.
///
/// # What has to hold
///
/// **File order is load-bearing.** The Python appends matches while iterating the FASTA, so
/// the records handed to the aligner are in **file** order. That is *not* `seqIDs` order:
/// `node.seq_ids` is a `BTreeSet<String>` ordered lexicographically, so it puts `"0_10_0"`
/// before `"0_2_0"` while the file does the opposite. Record order reaches mafft's input
/// file, and mafft's output depends on it — so [`output_sequence`] looks positions up and
/// then **sorts them ascending**, reproducing a file scan exactly. Never iterate `by_id`,
/// and never iterate `seq_ids` to build the output.
///
/// **Duplicate IDs.** `by_id` maps an ID to *every* position it occupies, so a FASTA with a
/// repeated ID yields the same repeated records a scan would. (The `ci` file has 0 duplicates
/// among its 20,516 records; the `Vec` makes the equivalence unconditional rather than
/// dataset-dependent.)
///
/// **Panics.** `isolate_num` is parsed and `isolate_list` indexed for *every* record, exactly
/// as the per-gene scan did, so a malformed or out-of-range record ID still aborts the run
/// even when no gene selects it.
pub struct CombinedDna {
    /// Records in file order.
    records: Vec<crate::support::seqio::SeqRecord>,
    /// `isolate_list[n].replace(";", "") + ";" + seq.id` per record, same order.
    ///
    /// Precomputed here because the Python recomputes it for every record on every gene:
    /// 5,116 × 20,516 = 105 M `format!` + `replace` pairs on `ci`.
    names: Vec<String>,
    /// Sequence ID -> its positions in `records`, ascending.
    by_id: std::collections::HashMap<String, Vec<usize>>,
}

impl CombinedDna {
    /// Parse `{outdir}combined_DNA_CDS.fasta` and precompute the per-record isolate names.
    ///
    /// Safe to hoist out of the per-gene loop because nothing in the alignment stage writes
    /// this file: it is produced by `prokka::output_files` and appended by `find_missing`,
    /// both long before `main` reaches the alignment stage, and the alignment workers only
    /// write into `temp_directory` and `outdir/aligned_gene_sequences/`.
    pub fn load(outdir: &str, isolate_list: &[String]) -> CombinedDna {
        // `isolate_list[i].replace(";", "")` depends only on `i`, so it is computed once per
        // isolate rather than once per record. Indexing is unchanged, so an out-of-range
        // isolate_num panics exactly as before.
        let clean_isolates: Vec<String> = isolate_list.iter().map(|s| s.replace(';', "")).collect();

        let records =
            crate::support::seqio::parse_fasta_file(&format!("{outdir}combined_DNA_CDS.fasta"));

        let mut names: Vec<String> = Vec::with_capacity(records.len());
        let mut by_id: std::collections::HashMap<String, Vec<usize>> =
            std::collections::HashMap::with_capacity(records.len());
        for (i, seq) in records.iter().enumerate() {
            let isolate_num: usize = seq.id.split('_').next().unwrap().parse().unwrap();
            names.push(format!("{};{}", clean_isolates[isolate_num], seq.id));
            by_id.entry(seq.id.clone()).or_default().push(i);
        }

        CombinedDna {
            records,
            names,
            by_id,
        }
    }
}

/// `generate_alignments.py::output_sequence`
///
/// Returns the written path, or `None` when the gene had a single sequence (in which case
/// the record is written straight to the aligned directory and no alignment is scheduled).
///
/// Takes a pre-parsed [`CombinedDna`] in place of the Python's `isolate_list`, which is the
/// fix for PORTING_PLAN.md §9 item 1 — see that type for why it is equivalent.
pub fn output_sequence(
    node: &NodeAttrs,
    combined_dna: &CombinedDna,
    temp_directory: &str,
    outdir: &str,
) -> Option<String> {
    use crate::support::seqio::{write_fasta_file, SeqRecord};

    // Get the name of the sequences for the gene of interest
    let sequence_ids = &node.seq_ids;
    let mut output_sequences: Vec<SeqRecord> = Vec::new();
    // Counter for the number of sequences for downstream check of >1
    let mut isolate_no = 0usize;

    // Look for gene sequences among all genes. The Python scans the whole combined FASTA
    // here and keeps the records whose id is in `seqIDs`, which yields them in FILE order;
    // sorting the looked-up positions reproduces that scan exactly. Iterating `sequence_ids`
    // instead would emit them in BTreeSet (lexicographic) order, which is a different order
    // and reaches the aligner's input file.
    let mut positions: Vec<usize> = Vec::new();
    for seq_id in sequence_ids.iter() {
        if let Some(idxs) = combined_dna.by_id.get(seq_id.as_str()) {
            positions.extend_from_slice(idxs);
        }
    }
    positions.sort_unstable();

    for i in positions {
        output_sequences.push(SeqRecord::new(
            combined_dna.records[i].seq.clone(),
            combined_dna.names[i].clone(),
            String::new(),
        ));
        isolate_no += 1;
    }

    // set filename to gene name, if more than one sequence to be aligned
    let outname = if isolate_no > 1 {
        get_temp_dna_input_path(node, temp_directory)
    } else {
        // If only one sequence, output it to the aligned directory and break
        let outname = get_expected_gene_alignment_path(node, outdir, false);
        write_fasta_file(&output_sequences, &outname);
        return None;
    };

    // Write them to disk
    write_fasta_file(&output_sequences, &outname);
    Some(outname)
}

/// `generate_alignments.py::output_dna_and_protein`
///
/// Returns `(protein_path, dna_path)`; the caller filters on `x[0]` being truthy, so both
/// slots are `None` when the gene needs no MSA.
pub fn output_dna_and_protein(
    node: &NodeAttrs,
    isolate_list: &[String],
    temp_directory: &str,
    outdir: &str,
    all_proteins: &std::collections::HashMap<String, SeqRecord>,
    all_dna: &std::collections::HashMap<String, SeqRecord>,
) -> (Option<String>, Option<String>) {
    use crate::support::seqio::{write_fasta_file, SeqRecord};

    let sequence_ids = &node.seq_ids;
    let mut output_dna: Vec<SeqRecord> = Vec::new();
    let mut output_protein: Vec<SeqRecord> = Vec::new();
    let mut isolate_no = 0usize;

    for seq_id in sequence_ids.iter() {
        let isolate_num: usize = seq_id.split('_').next().unwrap().parse().unwrap();
        let isolate_name = format!("{};{}", isolate_list[isolate_num].replace(';', ""), seq_id);
        output_dna.push(SeqRecord::new(
            all_dna[seq_id].seq.clone(),
            isolate_name.clone(),
            String::new(),
        ));
        output_protein.push(SeqRecord::new(
            all_proteins[seq_id].seq.clone(),
            isolate_name,
            String::new(),
        ));
        isolate_no += 1;
    }

    if isolate_no > 1 {
        let prot_name = get_expected_protein_input_path(node, temp_directory);
        let dna_name = get_expected_unaligned_dna_path(node, outdir);
        write_fasta_file(&output_protein, &prot_name);
        write_fasta_file(&output_dna, &dna_name);
        (Some(prot_name), Some(dna_name))
    } else {
        // a single-sequence gene needs no MSA -- write the final output directly
        let outname = get_expected_gene_alignment_path(node, outdir, true);
        write_fasta_file(&output_dna, &outname);
        (None, None)
    }
}

// --- aligner command construction -------------------------------------------------------

/// One entry of the command lists returned by the `get_*_commands` functions.
///
/// The Python builds a list of `[command_string, input_path, output_path]`-shaped tuples
/// whose exact arity varies by aligner; `align_sequences` indexes `command[0]`, `[1]`, `[2]`.
#[derive(Debug, Clone, Default)]
pub struct AlignCommand {
    /// `command[0]` — the shell command, run verbatim. `None` for a single-sequence gene,
    /// which `align_sequences` skips.
    pub command: Option<String>,
    /// `command[1]` — the input file, removed after the alignment runs.
    pub input_path: Option<String>,
    /// `get_align_dna_to_alignment_commands` returns `command[0]` as an **argv list**, not
    /// a string, and `realign_dna_sequences` runs `command[0][:-1]` while writing stdout to
    /// `command[0][-1]`. Kept as a separate field so the two shapes stay distinguishable.
    pub argv: Option<Vec<String>>,
}

/// `generate_alignments.py::get_alignment_commands`
pub fn get_alignment_commands(
    fastafile_name: &str,
    outdir: &str,
    aligner: &str,
    _threads: i64,
) -> AlignCommand {
    let gene_name = gene_name_of(fastafile_name);
    let command = match aligner {
        "prank" => format!(
            "prank -d={fastafile_name} -o={outdir}aligned_gene_sequences/{gene_name} -f=8"
        ),
        "mafft" => format!("mafft --auto --adjustdirection --thread 1 --nuc {fastafile_name}"),
        "clustal" => format!(
            "clustalo  -i {fastafile_name} -t DNA --threads 1 -o {outdir}aligned_gene_sequences/{gene_name}.aln.fas"
        ),
        "muscle" => format!(
            "muscle  -align {fastafile_name} -nt  -threads 1 -output {outdir}aligned_gene_sequences/{gene_name}.aln.fas"
        ),
        "muscle-super5" => format!(
            "muscle  -super5 {fastafile_name} -nt  -threads 1 -output {outdir}aligned_gene_sequences/{gene_name}.aln.fas"
        ),
        // FAMSA only supports amino acids; check_aligner_sanity should have caught this.
        "famsa" => panic!(
            "RuntimeError: FAMSA2 only supports amino-acid alignment.\
             Use --codons or --strict-codons to align."
        ),
        other => panic!("unknown aligner: {other}"),
    };
    AlignCommand {
        command: Some(command),
        input_path: Some(fastafile_name.to_string()),
        argv: None,
    }
}

/// `fastafile_name.split("/")[-1].split(".")[0]` — the basename up to the FIRST dot.
///
/// Note this is not `os.path.splitext`: a gene named `abc.def` truncates to `abc`.
fn gene_name_of(path: &str) -> String {
    path.rsplit('/')
        .next()
        .unwrap_or(path)
        .split('.')
        .next()
        .unwrap_or("")
        .to_string()
}

/// `generate_alignments.py::get_protein_commands`
pub fn get_protein_commands(
    fastafile_name: Option<&str>,
    outdir: &str,
    aligner: &str,
    _threads: i64,
) -> AlignCommand {
    let Some(fastafile_name) = fastafile_name else {
        return AlignCommand::default();
    };
    let gene_name = gene_name_of(fastafile_name);
    let command = match aligner {
        "prank" => format!("prank -d={fastafile_name} -o={gene_name} -f=8"),
        "mafft" => format!("mafft --auto --amino {fastafile_name}"),
        "clustal" => format!(
            "clustalo  -i {fastafile_name} -t Protein --threads 1 -o {outdir}aligned_protein_sequences/{gene_name}.aln.fas"
        ),
        "muscle" => format!(
            "muscle  -align {fastafile_name} -amino  -threads 1 -output {outdir}aligned_protein_sequences/{gene_name}.aln.fas"
        ),
        "muscle-super5" => format!(
            "muscle  -super5 {fastafile_name} -amino  -threads 1 -output {outdir}aligned_protein_sequences/{gene_name}.aln.fas"
        ),
        "famsa" => format!(
            "famsa  -t 1 {fastafile_name} {outdir}aligned_protein_sequences/{gene_name}.aln.fas"
        ),
        other => panic!("unknown aligner: {other}"),
    };
    AlignCommand {
        command: Some(command),
        input_path: Some(fastafile_name.to_string()),
        argv: None,
    }
}

/// `generate_alignments.py::get_align_dna_to_alignment_commands`
pub fn get_align_dna_to_alignment_commands(
    bad_dna_seqs_file: &str,
    codonalignment_file: &str,
    outdir: &str,
    aligner: &str,
) -> AlignCommand {
    let gene_name = gene_name_of(codonalignment_file);
    let argv: Vec<String> = match aligner {
        "prank" => {
            panic!("Exception: This is a bug! Panaroo does not supports codon alignment with PRANK")
        }
        // default to MAFFT for profile alignment (other aligners do not support it)
        "mafft" | "muscle" | "muscle-super5" | "famsa" => vec![
            "mafft".into(),
            "--add".into(),
            bad_dna_seqs_file.to_string(),
            codonalignment_file.to_string(),
            format!("{outdir}aligned_gene_sequences/{gene_name}.aln.fas"),
        ],
        "clustal" => vec![
            "clustalo".into(),
            "--in".into(),
            bad_dna_seqs_file.to_string(),
            "--profile1".into(),
            codonalignment_file.to_string(),
            "--out".into(),
            // NOTE: upstream omits the '.' before "aln.fas" on this branch only.
            format!("{outdir}aligned_gene_sequences/{gene_name}aln.fas"),
        ],
        other => panic!("unknown aligner: {other}"),
    };
    AlignCommand {
        command: None,
        input_path: Some(bad_dna_seqs_file.to_string()),
        argv: Some(argv),
    }
}

/// `generate_alignments.py::align_sequences`
///
/// Runs one aligner invocation. Aligners differ in whether they write their own output file
/// or print to stdout; the prank branch also renames `*.best.fas`.
pub fn align_sequences(command: &AlignCommand, outdir: &str, aligner: &str) -> bool {
    // Avoid running alignments on single-isolate genes
    let Some(cmd) = &command.command else {
        return false;
    };

    if aligner == "mafft" {
        // mafft writes to stdout; the output file name is derived from the last token of
        // the command, which is the input path.
        let name = gene_name_of(cmd.split_whitespace().last().unwrap_or(""));
        #[cfg(feature = "mafft-embedded")]
        let stdout = crate::mafft_embedded::run_mafft(cmd);
        #[cfg(not(feature = "mafft-embedded"))]
        let (stdout, _stderr) = crate::support::proc::popen_communicate(cmd);
        std::fs::write(format!("{outdir}{name}.aln.fas"), stdout).expect("write alignment");
    } else {
        let r = crate::support::proc::run_shell_capture(cmd);
        if r.returncode != 0 {
            panic!("RuntimeError: {}", String::from_utf8_lossy(&r.stderr));
        }
    }
    if let Some(p) = &command.input_path {
        let _ = std::fs::remove_file(p);
    }
    true
}

/// `generate_alignments.py::realign_dna_sequences`
pub fn realign_dna_sequences(command: &AlignCommand, _outdir: &str, aligner: &str) -> bool {
    let argv = command.argv.as_ref().expect("profile-alignment argv");
    match aligner {
        "prank" => panic!(
            "Exception: This is a bug! Please report it. Panaroo does not support codon \
             alignment with PRANK"
        ),
        "mafft" | "muscle" | "muscle-super5" | "famsa" => {
            // run argv[:-1] and write stdout to argv[-1]
            let out_path = argv.last().unwrap();
            let (stdout, _err) = crate::support::proc::popen_argv(&argv[..argv.len() - 1]);
            std::fs::write(out_path, stdout).expect("write realignment");
        }
        "clustal" => {
            let r = crate::support::proc::run_argv(argv);
            if r.returncode != 0 {
                panic!("RuntimeError: {}", String::from_utf8_lossy(&r.stderr));
            }
        }
        other => panic!("unknown aligner: {other}"),
    }
    // Delete the bad DNA seqs file
    if let Some(p) = &command.input_path {
        let _ = std::fs::remove_file(p);
    }
    true
}

/// `generate_alignments.py::multi_align_sequences`
pub fn multi_align_sequences(
    commands: &[AlignCommand],
    outdir: &str,
    threads: i64,
    aligner: &str,
) -> Vec<bool> {
    let outdir = outdir.to_string();
    let aligner = aligner.to_string();
    crate::support::parallel::parallel_map(threads, commands.to_vec(), move |c| {
        align_sequences(&c, &outdir, &aligner)
    })
}

/// `generate_alignments.py::multi_realign_sequences`
pub fn multi_realign_sequences(
    commands: &[AlignCommand],
    outdir: &str,
    threads: i64,
    aligner: &str,
) -> Vec<bool> {
    let outdir = outdir.to_string();
    let aligner = aligner.to_string();
    crate::support::parallel::parallel_map(threads, commands.to_vec(), move |c| {
        realign_dna_sequences(&c, &outdir, &aligner)
    })
}

// --- codon alignment --------------------------------------------------------------------

/// `generate_alignments.py::read_sequences`
pub fn read_sequences(handle: &str) -> Vec<SeqRecord> {
    crate::support::seqio::parse_fasta_file(handle)
}

/// `generate_alignments.py::read_alignment`
pub fn read_alignment(handle: &str) -> MultipleSeqAlignment {
    crate::support::align_io::read_fasta_alignment(handle)
}

/// `generate_alignments.py::reorder_protein_alignment_to_match_dna`
///
/// Raises on duplicate or mismatched IDs.
pub fn reorder_protein_alignment_to_match_dna(
    dna_records: &[SeqRecord],
    protein_alignment: &MultipleSeqAlignment,
    gene_name: &str,
) -> MultipleSeqAlignment {
    use std::collections::HashSet;
    let dna_ids: Vec<&str> = dna_records.iter().map(|r| r.id.as_str()).collect();
    let protein_ids: Vec<&str> = protein_alignment
        .records
        .iter()
        .map(|r| r.id.as_str())
        .collect();

    let dna_set: HashSet<&str> = dna_ids.iter().copied().collect();
    let prot_set: HashSet<&str> = protein_ids.iter().copied().collect();
    if dna_set.len() != dna_ids.len() {
        panic!("ValueError: Duplicate DNA sequence IDs found for gene: {gene_name}");
    }
    if prot_set.len() != protein_ids.len() {
        panic!("ValueError: Duplicate protein sequence IDs found for gene: {gene_name}");
    }
    if dna_set != prot_set {
        panic!("ValueError: DNA and protein sequence IDs do not match for gene: {gene_name}");
    }

    let by_id: std::collections::HashMap<&str, &SeqRecord> = protein_alignment
        .records
        .iter()
        .map(|r| (r.id.as_str(), r))
        .collect();
    MultipleSeqAlignment {
        records: dna_records
            .iter()
            .map(|r| by_id[r.id.as_str()].clone())
            .collect(),
    }
}

/// `generate_alignments.py::multithread_codonalign_build`
///
/// Wraps `Bio.codonalign.build`, which is what makes the codon path expensive and warning-
/// noisy. Returns `None` when the build fails.
pub fn multithread_codonalign_build(
    protein: &MultipleSeqAlignment,
    dna: &[SeqRecord],
    name: &str,
) -> (String, MultipleSeqAlignment) {
    // `Bio.codonalign.build(protein, dna, codon_table=generic_by_id[11])`.
    //
    // For the inputs that reach here the mapping is mechanical, and OBSERVED to be exactly
    // that: every caller has already filtered to sequences whose translation equals the
    // ungapped protein, so `build` walks the aligned protein emitting three DNA bases per
    // residue and `---` per gap. Biopython sets the result's description to
    // `<unknown description>`; Panaroo blanks it immediately afterwards, so this leaves it
    // empty.
    //
    // If a caller ever passes a pair that does *not* satisfy that invariant, Biopython
    // would attempt a frameshift-tolerant alignment instead. Rather than approximate that,
    // this panics -- an approximation would be silently wrong.
    let mut records = Vec::with_capacity(dna.len());
    for (i, d) in dna.iter().enumerate() {
        let p = &protein.records[i];
        let db = d.seq.as_bytes();
        let mut out = String::with_capacity(p.seq.len() * 3);
        let mut k = 0usize;
        for c in p.seq.chars() {
            if c == '-' {
                out.push_str("---");
            } else {
                if k + 3 > db.len() {
                    panic!(
                        "codonalign::build: DNA exhausted for {name}/{} -- the protein does \
                         not correspond to the DNA (Biopython would try a frameshift \
                         alignment here; not reproduced)",
                        d.id
                    );
                }
                out.push_str(&d.seq[k..k + 3]);
                k += 3;
            }
        }
        records.push(SeqRecord::new(out, d.id.clone(), String::new()));
    }
    (name.to_string(), MultipleSeqAlignment { records })
}

/// `generate_alignments.py::unambiguous_degenerate_codons`
///
/// Codons whose third position can be any base without changing the residue. Used as a
/// cheap pre-filter: `Bio.codonalign` cannot cope with a degenerate codon, so a DNA
/// sequence containing one is rejected and realigned at the protein level instead.
pub const UNAMBIGUOUS_DEGENERATE_CODONS: [&str; 8] =
    ["ACN", "TCN", "CTN", "CCN", "CGN", "GTN", "GCN", "GGN"];

/// `generate_alignments.py::reverse_translate_sequences`
///
/// The codon-alignment driver: QC-filters DNA against protein, builds codon alignments,
/// then realigns the sequences that failed QC back onto the codon alignment.
#[allow(clippy::too_many_arguments)]
pub fn reverse_translate_sequences(
    protein_sequence_files: &[String],
    dna_sequence_files: &[String],
    strict: bool,
    outdir: &str,
    temp_directory: &str,
    aligner: &str,
    threads: i64,
    completed_alignments_found: usize,
) -> Vec<String> {
    use crate::support::align_io::write_fasta_alignment;
    use crate::support::pydict::PyDict;
    use crate::support::seqio::write_fasta_file;

    // Check that the dna and protein files match up
    for index in 0..protein_sequence_files.len() {
        let gene_id = gene_name_of(&protein_sequence_files[index]);
        if gene_id != gene_name_of(&dna_sequence_files[index]) {
            println!("{}", protein_sequence_files[index]);
            println!("{}", dna_sequence_files[index]);
            panic!("ValueError: DNA and protien sequence IDs do not match!");
        }
    }

    // Read in files (multithreaded)
    let dna_sequences: Vec<Vec<SeqRecord>> =
        crate::support::parallel::parallel_map(threads, dna_sequence_files.to_vec(), |x| {
            read_sequences(&x)
        });
    let protein_alignments: Vec<MultipleSeqAlignment> =
        crate::support::parallel::parallel_map(threads, protein_sequence_files.to_vec(), |x| {
            read_alignment(&x)
        });

    let mut clean_dna: Vec<Vec<SeqRecord>> = Vec::new();
    let mut clean_proteins: Vec<MultipleSeqAlignment> = Vec::new();
    let mut reject_dna_files: PyDict<String, String> = PyDict::new();

    print_stage_progress(
        "Getting sequences",
        completed_alignments_found,
        dna_sequences.len(),
        None,
    );

    let trans_table = get_trans_table(11);

    for index in 0..dna_sequences.len() {
        let dna = &dna_sequences[index];
        let gene_name = gene_name_of(&dna_sequence_files[index]);
        let protein =
            reorder_protein_alignment_to_match_dna(dna, &protein_alignments[index], &gene_name);
        let mut seqids_to_remove: Vec<String> = Vec::new();

        for seq_index in 0..dna.len() {
            // sequentially checked QC failure variables
            let nogapped_protein_seq = protein.records[seq_index].seq.replace('-', "");
            let dna_seq = &dna[seq_index].seq;

            // fail if the sequence is not divisible by 3
            let fail0 = dna_seq.len() % 3 != 0;

            let fail1 = if !fail0 {
                // fail if the translated sequence isn't the same as the protein
                translate(dna_seq, &trans_table).trim_matches('*') != nogapped_protein_seq
            } else {
                false
            };

            // fail if there is a run of > 1 unknown nucleotides
            let fail2 = if !fail1 {
                dna_seq.contains("NN")
            } else {
                false
            };

            // Fail if the DNA contains a degenerate codon, codonalign cannot cope
            let mut fail3 = false;
            if !fail1 && !fail2 {
                // cheaper filter first
                if dna_seq.contains('N') {
                    for codon in UNAMBIGUOUS_DEGENERATE_CODONS {
                        if dna_seq.contains(codon) {
                            fail3 = true;
                        }
                    }
                }
            }

            if fail0 || fail1 || fail2 || fail3 {
                // `list(set([dna_id, protein_id]))` -- the two ids are equal after the
                // reorder, so this contributes one entry; it is only membership-tested.
                let mut pair = vec![
                    dna[seq_index].id.clone(),
                    protein.records[seq_index].id.clone(),
                ];
                pair.sort();
                pair.dedup();
                seqids_to_remove.extend(pair);
            }
        }

        if !seqids_to_remove.is_empty() {
            let remove: std::collections::HashSet<&String> = seqids_to_remove.iter().collect();
            let mut reject_dna = Vec::new();
            let mut clean_nucs = Vec::new();
            for sequence in dna {
                if remove.contains(&sequence.id) {
                    reject_dna.push(sequence.clone());
                } else {
                    clean_nucs.push(sequence.clone());
                }
            }
            let clean_prots: Vec<SeqRecord> = protein
                .records
                .iter()
                .filter(|s| !remove.contains(&s.id))
                .cloned()
                .collect();

            clean_dna.push(clean_nucs);
            clean_proteins.push(MultipleSeqAlignment {
                records: clean_prots,
            });

            let reject_outname = format!("{temp_directory}{gene_name}_untrans_dna.fasta");
            write_fasta_file(&reject_dna, &reject_outname);
            reject_dna_files.insert(gene_name.clone(), reject_outname);
        } else {
            clean_dna.push(dna.clone());
            clean_proteins.push(protein);
        }
    }

    // build codon alignments
    print_stage_progress(
        "Reverse translating DNA",
        completed_alignments_found,
        clean_proteins.len(),
        None,
    );

    let jobs: Vec<(MultipleSeqAlignment, Vec<SeqRecord>, String)> = (0..clean_proteins.len())
        .map(|i| {
            (
                clean_proteins[i].clone(),
                clean_dna[i].clone(),
                gene_name_of(&dna_sequence_files[i]),
            )
        })
        .collect();
    let all_codon_alignments: Vec<(String, MultipleSeqAlignment)> =
        crate::support::parallel::parallel_map(threads, jobs, |(p, d, n)| {
            multithread_codonalign_build(&p, &d, &n)
        });

    let mut completed: PyDict<String, MultipleSeqAlignment> = PyDict::new();
    let mut missing: PyDict<String, MultipleSeqAlignment> = PyDict::new();
    for (name, aln) in all_codon_alignments {
        if reject_dna_files.contains_key(&name) {
            missing.insert(name, aln);
        } else {
            completed.insert(name, aln);
        }
    }
    // `<unknown description>` is stripped upstream; our records already carry "".

    // output successful codon alignments
    for (x, aln) in completed.items() {
        write_fasta_alignment(aln, &format!("{outdir}aligned_gene_sequences/{x}.aln.fas"));
    }

    if strict {
        // output alignments missing DNA as complete
        for (x, aln) in missing.items() {
            write_fasta_alignment(aln, &format!("{outdir}aligned_gene_sequences/{x}.aln.fas"));
        }
        return list_aligned_dir(outdir);
    }

    // output alignments missing some DNA sequences to tmpdir
    for (x, aln) in missing.items() {
        write_fasta_alignment(aln, &format!("{temp_directory}{x}.aln.fas"));
    }

    println!("{} DNA realignments to perform...", missing.len());

    let mut dna2codons_commands = Vec::new();
    for (gene_name, reject) in reject_dna_files.items() {
        dna2codons_commands.push(get_align_dna_to_alignment_commands(
            reject,
            &format!("{temp_directory}{gene_name}.aln.fas"),
            outdir,
            aligner,
        ));
    }

    println!("Aligning untranslatable DNA...");
    multi_realign_sequences(
        &dna2codons_commands,
        &format!("{outdir}aligned_gene_sequences/"),
        threads,
        aligner,
    );

    list_aligned_dir(outdir)
}

/// `sorted(os.listdir(outdir + "aligned_gene_sequences/"))` — Tier D sorts this, since raw
/// `listdir` is filesystem order.
fn list_aligned_dir(outdir: &str) -> Vec<String> {
    let dir = format!("{outdir}aligned_gene_sequences/");
    let mut v: Vec<String> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read_dir {dir}: {e}"))
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
        .collect();
    v.sort();
    v
}

/// `generate_alignments.py::write_alignment_header`
///
/// Writes `core_alignment_header.embl` / `pan_genome_alignment_header.embl`.
pub fn write_alignment_header(alignment_list: &[(String, usize)], outdir: &str, filename: &str) {
    use std::io::Write;
    let mut out_entries: Vec<String> = Vec::new();
    // Set the tracking variables for gene positions
    let mut gene_start = 1usize;
    let mut gene_end = 0usize;
    for (gene_name, gene_len) in alignment_list {
        gene_end += gene_len;
        out_entries.push(format!(
            "FT   feature         {gene_start}..{gene_end}\n\
             FT                   /label={gene_name}\n\
             FT                   /locus_tag={gene_name}\n"
        ));
        gene_start += gene_len;
    }
    // The header and footer carry fixed placeholder numbers upstream -- they are not
    // computed from the alignment. Reproduced verbatim.
    // NOTE: written as one literal on purpose. A `\`-continuation here would swallow the
    // leading whitespace of the next line, and the column alignment is part of the format.
    let header = "ID   Genome standard; DNA; PRO; 1234 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\n";
    let footer = "XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T; 0 other;\n//\n";

    let f = std::fs::File::create(format!("{outdir}{filename}")).expect("create embl header");
    let mut w = std::io::BufWriter::new(f);
    write!(w, "{header}").unwrap();
    for e in &out_entries {
        write!(w, "{e}").unwrap();
    }
    write!(w, "{footer}").unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::support::intbitset::IntBitSet;
    use std::collections::BTreeSet;

    /// Minimal node carrying just what the path builders read.
    fn node(name: &str, n_seqs: usize) -> NodeAttrs {
        NodeAttrs {
            size: n_seqs,
            centroid: vec![],
            max_len_id: 0,
            members: IntBitSet::default(),
            seq_ids: (0..n_seqs)
                .map(|i| format!("{i}_0_1"))
                .collect::<BTreeSet<_>>(),
            has_end: false,
            protein: vec![],
            dna: vec![],
            annotation: String::new(),
            description: String::new(),
            lengths: vec![],
            long_centroid_id: (0, String::new()),
            paralog: false,
            merged_dna: false,
            prev_centroids: None,
            name: Some(name.to_string()),
            genome_ids: None,
            gene_ids: None,
            degrees: None,
            high_var: None,
            gml_late_attrs_before_name: false,
        }
    }

    // Every expectation below was produced by running the reference Python.

    #[test]
    fn path_builders_match_python_for_a_short_name() {
        let multi = node("grpA", 2);
        assert_eq!(get_alignment_basename(&multi), "grpA");
        assert_eq!(get_temp_dna_input_path(&multi, "/t/"), "/t/grpA.fasta");
        assert_eq!(
            get_expected_gene_alignment_path(&multi, "/o", false),
            "/o/aligned_gene_sequences/grpA.aln.fas"
        );
        assert_eq!(
            get_expected_gene_alignment_path(&multi, "/o", true),
            "/o/aligned_gene_sequences/grpA.aln.fas"
        );
        assert_eq!(
            get_expected_protein_input_path(&multi, "/t"),
            "/t/grpA.fasta"
        );
        assert_eq!(
            get_expected_protein_alignment_path(&multi, "/o"),
            "/o/aligned_protein_sequences/grpA.aln.fas"
        );
        assert_eq!(
            get_expected_unaligned_dna_path(&multi, "/o"),
            "/o/unaligned_dna_sequences/grpA.fasta"
        );
        assert_eq!(
            get_resume_manifest_path("/o"),
            "/o/alignment_resume_state.json"
        );
        assert!(node_requires_msa(&multi));
    }

    #[test]
    fn single_sequence_gene_skips_the_aln_extension() {
        // python: get_expected_gene_alignment_path(node1, '/o', False)
        //           -> '/o/aligned_gene_sequences/grpB.fasta'
        //         get_expected_gene_alignment_path(node1, '/o', True)
        //           -> '/o/aligned_gene_sequences/grpB.aln.fas'
        let single = node("grpB", 1);
        assert_eq!(
            get_expected_gene_alignment_path(&single, "/o", false),
            "/o/aligned_gene_sequences/grpB.fasta"
        );
        assert_eq!(
            get_expected_gene_alignment_path(&single, "/o", true),
            "/o/aligned_gene_sequences/grpB.aln.fas"
        );
        assert!(!node_requires_msa(&single));
    }

    #[test]
    fn long_names_truncate_exactly_as_python_does() {
        // A 300-character name. python:
        //   len(get_alignment_basename(n))            == 236
        //   len(get_temp_dna_input_path(n, '/t/'))    == 254   <- 248 + len('.fasta')
        let long = node(&"g".repeat(300), 2);
        assert_eq!(get_alignment_basename(&long).len(), 236);
        let tmp = get_temp_dna_input_path(&long, "/t/");
        assert_eq!(tmp.len(), 254);
        assert!(tmp.ends_with(".fasta"));
    }

    #[test]
    fn join_matches_os_path_join() {
        // os.path.join('/o', 'a', 'b')  -> '/o/a/b'
        // os.path.join('/o/', 'a')      -> '/o/a'   (no doubled separator)
        // os.path.join('/o', '/abs')    -> '/abs'   (absolute component resets)
        assert_eq!(join("/o", &["a", "b"]), "/o/a/b");
        assert_eq!(join("/o/", &["a"]), "/o/a");
        assert_eq!(join("/o", &["/abs"]), "/abs");
    }
}
