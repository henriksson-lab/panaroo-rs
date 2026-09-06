//! Translation of `panaroo/panaroo/generate_output.py` — every user-facing output file.

use crate::support::graph::Graph;
use std::collections::HashMap;
use std::io::Write;

/// `(length, has_internal_stop, is_valid)` — the `ids_len_stop` values built in
/// `__main__::main` and threaded through here.
pub type IdLenStop = (usize, bool, bool);

/// One entry of `gene_alignments` in `concatenate_core_genome_alignments`:
/// `(gene_name, genome_id -> (record_id, aligned_seq), gene_length, entropy)`.
///
/// The Python builds this as a bare tuple; naming it keeps the call sites readable without
/// changing the shape.
type GeneAlignment = (
    String,
    crate::support::pydict::PyDict<String, (String, String)>,
    usize,
    f64,
);

fn push_gaps(seq: &mut String, n: usize) {
    seq.extend(std::iter::repeat('-').take(n));
}

fn write_delimited_cells<'a, W, I>(writer: &mut W, cells: I, sep: char)
where
    W: Write,
    I: IntoIterator<Item = &'a str>,
{
    let mut first = true;
    for cell in cells {
        if first {
            first = false;
        } else {
            write!(writer, "{sep}").unwrap();
        }
        write!(writer, "{cell}").unwrap();
    }
    writeln!(writer).unwrap();
}

fn seq_sample_key(seq_id: &str) -> &str {
    let last = seq_id.rfind('_').unwrap_or(seq_id.len());
    let prefix = &seq_id[..last];
    match prefix.rfind('_') {
        Some(second_last) => &seq_id[..second_last],
        None => "",
    }
}

fn roary_gene_name(annotation: &str) -> String {
    let mut out = String::new();
    let mut first = true;
    for gene in annotation
        .trim()
        .trim_matches(';')
        .split(';')
        .filter(|gn| !gn.is_empty())
    {
        if first {
            first = false;
        } else {
            out.push_str("~~~");
        }
        out.extend(
            gene.chars()
                .filter(|e| e.is_alphanumeric() || *e == '_' || *e == '~'),
        );
    }
    out
}

fn length_stats_and_first_mode(
    lengths: &[usize],
    counts: &mut Vec<(usize, usize)>,
) -> (usize, usize, f64, usize) {
    counts.clear();
    let mut min = lengths[0];
    let mut max = lengths[0];
    let mut sum = 0.0;
    for &length in lengths {
        min = min.min(length);
        max = max.max(length);
        sum += length as f64;
        match counts.iter_mut().find(|(value, _)| *value == length) {
            Some((_, count)) => *count += 1,
            None => counts.push((length, 1)),
        }
    }
    let mut mode = counts[0].0;
    let mut best_count = counts[0].1;
    for &(length, count) in &counts[1..] {
        if count > best_count {
            mode = length;
            best_count = count;
        }
    }
    (min, max, sum / lengths.len() as f64, mode)
}

fn write_core_alignment_records(
    path: &str,
    isolates: &std::collections::BTreeSet<String>,
    gene_alignments: &[GeneAlignment],
    seq_capacity: usize,
    hc_threshold: Option<f64>,
) -> usize {
    use crate::support::seqio::write_fasta_record;

    let file =
        std::fs::File::create(path).unwrap_or_else(|e| panic!("could not create {path}: {e}"));
    let mut writer = std::io::BufWriter::new(file);
    let mut keep_count = 0usize;
    match hc_threshold {
        None => {
            for iso in isolates {
                let mut seq = String::with_capacity(seq_capacity);
                for gene in gene_alignments {
                    match gene.1.get(iso) {
                        Some((_, s)) => seq.push_str(s),
                        None => push_gaps(&mut seq, gene.2),
                    }
                }
                write_fasta_record(iso, &seq, &mut writer);
            }
        }
        Some(threshold) => {
            for iso in isolates {
                let mut seq = String::with_capacity(seq_capacity);
                for gene in gene_alignments {
                    if gene.3 <= threshold {
                        keep_count += 1;
                        match gene.1.get(iso) {
                            Some((_, s)) => seq.push_str(s),
                            None => push_gaps(&mut seq, gene.2),
                        }
                    }
                }
                write_fasta_record(iso, &seq, &mut writer);
            }
        }
    }
    keep_count
}

/// `generate_output.py::generate_roary_gene_presence_absence`
///
/// Writes `gene_presence_absence_roary.csv`, `gene_presence_absence.csv` and
/// `gene_presence_absence.Rtab`. Also assigns `G.nodes[node]['name']`, which every
/// downstream path uses to build filenames.
///
/// **The most parity-sensitive function in the codebase.** Three separate hazards:
///
///  - `for node in component` iterates a Python `set[usize]` from `connected_components`,
///    fixing row order and the `Genome Fragment` / `Order within Fragment` columns
///    (PORTING_PLAN.md §6.1).
///  - `for seq in G.nodes[node]["seqIDs"]` iterates a Python `set[String]`; the first
///    element seen claims the primary slot in each cell and the rest are appended after
///    `;` in iteration order (§6.1).
///  - `max(lengths, key=lengths.count)` returns the **first** modal value (§6.4);
///  - `np.mean(lengths)` needs Python-compatible float formatting (§6.2, §6.3).
pub fn generate_roary_gene_presence_absence(
    g: &mut Graph,
    mems_to_isolates: &crate::support::pydict::PyDict<usize, String>,
    orig_ids: &HashMap<String, String>,
    ids_len_stop: &HashMap<String, IdLenStop>,
    output_dir: &str,
) {
    use crate::support::pyfmt::py_str_f64;
    use std::collections::HashMap as Map;

    // arange isolates
    let mut isolates: Vec<String> = Vec::new();
    let mut mems_to_index: Map<String, usize> = Map::new();
    for (i, (mem, name)) in mems_to_isolates.items().enumerate() {
        isolates.push(name.clone());
        mems_to_index.insert(mem.to_string(), i);
    }

    let mut roary = std::io::BufWriter::new(
        std::fs::File::create(format!("{output_dir}gene_presence_absence_roary.csv")).unwrap(),
    );
    let mut csv = std::io::BufWriter::new(
        std::fs::File::create(format!("{output_dir}gene_presence_absence.csv")).unwrap(),
    );
    let mut rtab = std::io::BufWriter::new(
        std::fs::File::create(format!("{output_dir}gene_presence_absence.Rtab")).unwrap(),
    );

    let header = [
        "Gene",
        "Non-unique Gene name",
        "Annotation",
        "No. isolates",
        "No. sequences",
        "Avg sequences per isolate",
        "Genome Fragment",
        "Order within Fragment",
        "Accessory Fragment",
        "Accessory Order with Fragment",
        "QC",
        "Min group size nuc",
        "Max group size nuc",
        "Avg group size nuc",
    ];
    write_delimited_cells(
        &mut roary,
        header
            .iter()
            .copied()
            .chain(isolates.iter().map(String::as_str)),
        ',',
    );
    write_delimited_cells(
        &mut csv,
        header[..3]
            .iter()
            .copied()
            .chain(isolates.iter().map(String::as_str)),
        ',',
    );
    write_delimited_cells(
        &mut rtab,
        std::iter::once("Gene").chain(isolates.iter().map(String::as_str)),
        '\t',
    );

    // Iterate through components writing out to file
    let mut used_gene_names: std::collections::HashSet<String> =
        std::collections::HashSet::from([String::new()]);
    let mut unique_id_count = 0usize;
    let mut frag = 0usize;
    let mut entry_list: Vec<Vec<String>> = Vec::new();
    let mut entry_ext_list: Vec<Vec<String>> = Vec::new();
    let mut entry_sizes: Vec<(usize, usize)> = Vec::new();
    let mut entry_count = 0usize;
    let mut length_counts: Vec<(usize, usize)> = Vec::new();

    for component in crate::support::graph::connected_components(g) {
        frag += 1;
        let mut count = 0usize;
        // Tier D: sorted(component)
        for node in component {
            count += 1;
            let lengths = &g.node(node).lengths;
            let (length_min, length_max, length_mean, len_mode) =
                length_stats_and_first_mode(lengths, &mut length_counts);

            let name = roary_gene_name(&g.node(node).annotation);

            let mut entry: Vec<String> = Vec::new();
            if !used_gene_names.contains(&name.to_lowercase()) {
                entry.push(name.clone());
                used_gene_names.insert(name.to_lowercase());
                g.node_mut(node).name = Some(name.clone());
            } else {
                let n = format!("group_{unique_id_count}");
                g.node_mut(node).name = Some(n.clone());
                entry.push(n);
                unique_id_count += 1;
            }

            let nd = g.node(node);
            entry.push(nd.annotation.clone());
            entry.push(nd.description.clone());
            entry.push(nd.size.to_string());
            entry.push(nd.seq_ids.len().to_string());
            entry.push(py_str_f64(nd.seq_ids.len() as f64 / nd.size as f64));
            entry.push(frag.to_string());
            entry.push(count.to_string());
            entry.push(String::new());
            entry.push(String::new());
            entry.push(String::new());
            entry.push(length_min.to_string());
            entry.push(length_max.to_string());
            entry.push(py_str_f64(length_mean));

            let mut pres_abs = vec![String::new(); isolates.len()];
            let mut pres_abs_ext = vec![String::new(); isolates.len()];
            let mut entry_size = 0usize;
            // Tier D: sorted(seqIDs) -- a BTreeSet is already sorted
            for seq in nd.seq_ids.iter() {
                let sample_key = seq_sample_key(seq);
                let sample_id = *mems_to_index.get(sample_key).unwrap();
                let val = orig_ids.get(seq).cloned().unwrap_or_else(|| seq.clone());
                if pres_abs[sample_id].is_empty() {
                    pres_abs[sample_id] = val.clone();
                    pres_abs_ext[sample_id] = val;
                    entry_size += 1;
                } else {
                    // this is similar to PIRATE output
                    pres_abs[sample_id].push(';');
                    pres_abs[sample_id].push_str(&val);
                    pres_abs_ext[sample_id].push(';');
                    pres_abs_ext[sample_id].push_str(&val);
                }
                let ils = &ids_len_stop[seq];
                if ((ils.0 as f64 - len_mode as f64).abs() / len_mode as f64)
                    > (0.05 * len_mode as f64)
                {
                    pres_abs_ext[sample_id].push_str("_len");
                }
                if !ils.2 {
                    pres_abs_ext[sample_id].push_str("_pseudo");
                }
            }

            let mut ext = entry[..3].to_vec();
            ext.extend(pres_abs_ext);
            entry.extend(pres_abs.iter().cloned());
            entry_list.push(entry);
            entry_ext_list.push(ext);
            entry_sizes.push((entry_size, entry_count));
            entry_count += 1;
        }
    }

    // sort so that the most common genes are first (as in roary)
    entry_sizes.sort_by(|a, b| b.cmp(a));
    for (_s, i) in entry_sizes {
        write_delimited_cells(&mut roary, entry_list[i].iter().map(String::as_str), ',');
        write_delimited_cells(&mut csv, entry_ext_list[i].iter().map(String::as_str), ',');
        write!(rtab, "{}\t", entry_list[i][0]).unwrap();
        for (j, entry) in entry_list[i][14..].iter().enumerate() {
            if j != 0 {
                write!(rtab, "\t").unwrap();
            }
            let call = if entry.is_empty() { "0" } else { "1" };
            write!(rtab, "{call}").unwrap();
        }
        writeln!(rtab).unwrap();
    }
}

/// `generate_output.py::generate_pan_genome_reference`
///
/// Writes `pan_genome_reference.fa`.
pub fn generate_pan_genome_reference(
    g: &Graph,
    output_dir: &str,
    ids_len_stop: &HashMap<String, IdLenStop>,
    split_paralogs: bool,
) {
    use crate::support::seqio::{write_fasta_file, SeqRecord};
    use std::collections::HashSet;

    // need to treat paralogs differently?
    let mut centroids: HashSet<String> = HashSet::new();
    let mut records: Vec<SeqRecord> = Vec::new();
    let mut representatives: HashMap<String, String> = HashMap::new();

    for node in g.nodes() {
        let nd = g.node(node);
        if !split_paralogs && nd.centroid.iter().any(|c| centroids.contains(c)) {
            continue;
        }

        let mut best = nd.centroid[0].clone();
        for centroid in &nd.centroid {
            // skip sequences that are not valid genes
            if !ids_len_stop[centroid].2 {
                continue;
            }
            if ids_len_stop[centroid].0 > ids_len_stop[&best].0 {
                best = centroid.clone();
            }
            centroids.insert(centroid.clone());
        }
        representatives.insert(best, nd.name.clone().expect("node name set"));
    }

    let text =
        std::fs::read_to_string(format!("{output_dir}gene_data.csv")).expect("read gene_data.csv");
    for line in text.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() < 6 {
            continue;
        }
        if let Some(name) = representatives.get(f[2]) {
            records.push(SeqRecord::new(
                f[5].to_string(),
                name.clone(),
                String::new(),
            ));
        }
    }

    write_fasta_file(&records, &format!("{output_dir}pan_genome_reference.fa"));
}

/// `generate_output.py::generate_common_struct_presence_absence`
///
/// Writes `struct_presence_absence.Rtab`.
pub fn generate_common_struct_presence_absence(
    g: &Graph,
    output_dir: &str,
    mems_to_isolates: &crate::support::pydict::PyDict<usize, String>,
    min_variant_support: usize,
) {
    use crate::support::pydict::PyDict;
    use std::io::Write;

    // arange isolates
    let mut isolates: Vec<String> = Vec::new();
    let mut members: Vec<usize> = Vec::new();
    for (mem, name) in mems_to_isolates.items() {
        isolates.push(name.clone());
        members.push(*mem);
    }

    // `struct_variants` is a plain dict keyed by a 3-tuple -- insertion ordered, so the
    // column order in the Rtab follows node order then edge-pair order.
    let mut struct_variants: PyDict<(usize, usize, usize), crate::support::intbitset::IntBitSet> =
        PyDict::new();
    for node in g.nodes() {
        if g.degree(node) < 3 {
            continue; // skip as linear
        }
        let edges = g.edges_of(&[node]);
        for i in 0..edges.len() {
            for j in (i + 1)..edges.len() {
                let (p0, p1) = (edges[i], edges[j]);
                let in_both = g
                    .edge(p0.0, p0.1)
                    .members
                    .intersection(&g.edge(p1.0, p1.1).members);
                if in_both.len() >= min_variant_support {
                    struct_variants.insert((p0.0, p0.1, p1.1), in_both);
                }
            }
        }
    }

    let mut out = std::io::BufWriter::new(
        std::fs::File::create(format!("{output_dir}struct_presence_absence.Rtab")).unwrap(),
    );
    write!(out, "Gene").unwrap();
    for iso in &isolates {
        write!(out, "\t{iso}").unwrap();
    }
    writeln!(out).unwrap();
    for (variant, in_both) in struct_variants.items() {
        write!(
            out,
            "{}-{}-{}",
            g.node(variant.1).name.as_deref().unwrap_or(""),
            g.node(variant.0).name.as_deref().unwrap_or(""),
            g.node(variant.2).name.as_deref().unwrap_or("")
        )
        .unwrap();
        for &member in &members {
            let call = if in_both.contains(member) { "1" } else { "0" };
            write!(out, "\t{call}").unwrap();
        }
        writeln!(out).unwrap();
    }
}

/// `generate_output.py::generate_pan_genome_alignment`
#[allow(clippy::too_many_arguments)]
pub fn generate_pan_genome_alignment(
    g: &Graph,
    temp_dir: &str,
    output_dir: &str,
    threads: i64,
    aligner: &str,
    codons: bool,
    strict: bool,
    isolates: &[String],
    resume: bool,
) {
    use crate::generate_alignments::*;

    // Make a folder for the output alignments
    let _ = std::fs::create_dir(format!("{output_dir}aligned_gene_sequences"));

    let gene_ids = g.nodes();
    let nodes: Vec<(usize, &crate::support::graph::NodeAttrs)> =
        gene_ids.iter().map(|&id| (id, g.node(id))).collect();

    let pending_gene_ids = get_pending_gene_ids(&nodes, output_dir, codons, resume);
    let total_gene_count = gene_ids.len();

    if codons || strict {
        let (protein_pending, reverse_translate_pending) =
            get_pending_codon_gene_ids(&nodes, output_dir, resume);
        println!("Codon alignment is experimental in Panaroo...");
        let _ = std::fs::create_dir(format!("{output_dir}aligned_protein_sequences"));
        let _ = std::fs::create_dir(format!("{output_dir}unaligned_dna_sequences"));

        // transform to dicts for fast lookup
        let proteins_dic = fasta_by_id(&format!("{output_dir}combined_protein_CDS.fasta"));
        let nucleotides_dic = fasta_by_id(&format!("{output_dir}combined_DNA_CDS.fasta"));
        let clean_isolates: Vec<String> = isolates.iter().map(|s| s.replace(';', "")).collect();

        // File output must stay single threaded (see the upstream comment)
        let mut output_files = Vec::new();
        for gene in &protein_pending {
            output_files.push(output_dna_and_protein(
                g.node(*gene),
                &clean_isolates,
                temp_dir,
                output_dir,
                &proteins_dic,
                &nucleotides_dic,
            ));
        }
        let filtered: Vec<(String, String)> = output_files
            .into_iter()
            .filter_map(|(p, d)| match (p, d) {
                (Some(p), Some(d)) => Some((p, d)),
                _ => None,
            })
            .collect();
        let unaligned_protein_files: Vec<String> = filtered.iter().map(|x| x.0.clone()).collect();

        let commands: Vec<AlignCommand> = unaligned_protein_files
            .iter()
            .map(|f| get_protein_commands(Some(f), output_dir, aligner, threads))
            .collect();
        print_stage_progress(
            "Protein alignments",
            total_gene_count - protein_pending.len(),
            commands.len(),
            Some(total_gene_count),
        );
        multi_align_sequences(
            commands,
            &format!("{output_dir}aligned_protein_sequences/"),
            threads,
            aligner,
        );

        let (protein_sequences, unaligned_dna_files) =
            get_codon_pending_files(&nodes, output_dir, &reverse_translate_pending);

        // Check all alignments completed
        for file in &protein_sequences {
            if !std::path::Path::new(file).is_file() {
                println!("{file}");
                panic!("RuntimeError: Some alignments failed to complete!");
            }
        }

        if !reverse_translate_pending.is_empty() {
            let completed = total_gene_count - reverse_translate_pending.len();
            reverse_translate_sequences(
                &protein_sequences,
                &unaligned_dna_files,
                strict,
                output_dir,
                temp_dir,
                aligner,
                threads,
                completed,
            );
        }
    } else {
        let mut temp_dir = temp_dir.to_string();
        if aligner == "none" {
            temp_dir = format!("{output_dir}unaligned_gene_sequences/");
            let _ = std::fs::create_dir(&temp_dir);
        }

        print_stage_progress(
            "Gene alignments",
            total_gene_count - pending_gene_ids.len(),
            pending_gene_ids.len(),
            Some(total_gene_count),
        );

        let td = temp_dir.clone();
        let od = output_dir.to_string();
        // A borrow, not a clone: `output_sequence` only reads the node, and `NodeAttrs`
        // owns the `dna`/`protein` sequence lists -- cloning every pending node duplicated
        // the whole gene set for the duration of the stage.
        let jobs: Vec<&crate::support::graph::NodeAttrs> =
            pending_gene_ids.iter().map(|&x| g.node(x)).collect();
        // PORTING_PLAN.md §9 item 1: parse combined_DNA_CDS.fasta ONCE for the stage rather
        // than once per gene inside `output_sequence`. Guarded on non-empty so that a run
        // with nothing pending still never opens the file, as before.
        let unaligned: Vec<Option<String>> = if jobs.is_empty() {
            Vec::new()
        } else {
            let combined_dna = CombinedDna::load(output_dir, isolates);
            crate::support::parallel::parallel_map(threads, jobs, move |node| {
                output_sequence(node, &combined_dna, &td, &od)
            })
        };
        let unaligned: Vec<String> = unaligned.into_iter().flatten().collect();

        if aligner == "none" {
            println!("No aligner specified. Returning unaligned gene fasta files.");
            return;
        }

        let commands: Vec<AlignCommand> = unaligned
            .iter()
            .map(|f| get_alignment_commands(f, output_dir, aligner, threads))
            .collect();
        multi_align_sequences(
            commands,
            &format!("{output_dir}aligned_gene_sequences/"),
            threads,
            aligner,
        );
    }
}

/// `dict(zip([x.id for x in recs], recs))` over a FASTA file.
///
/// Not a Python function -- the Python inlines this twice in the codon branch. Kept as one
/// helper because it is a lookup table, not logic.
fn fasta_by_id(path: &str) -> HashMap<String, crate::support::seqio::SeqRecord> {
    crate::support::seqio::parse_fasta_file(path)
        .into_iter()
        .map(|r| (r.id.clone(), r))
        .collect()
}

/// `generate_output.py::get_core_gene_nodes`
///
/// Nodes present in at least `threshold * num_isolates` genomes. With `subset`, takes a
/// `random.sample` — which needs a Mersenne Twister clone for parity on that path
/// (PORTING_PLAN.md §6.7).
pub fn get_core_gene_nodes(
    g: &Graph,
    threshold: f64,
    num_isolates: usize,
    subset: Option<usize>,
) -> Vec<usize> {
    // Get the core genes based on percent threshold.
    // Note the Python compares `float(size) / float(num_isolates) >= threshold`, a division
    // -- not `size >= threshold * num_isolates`. With float thresholds like 0.95 the two
    // disagree at the boundary, so the division is transcribed literally.
    let mut core_nodes: Vec<usize> = Vec::new();
    for node in g.nodes() {
        if g.node(node).size as f64 / num_isolates as f64 >= threshold {
            core_nodes.push(node);
        }
    }
    if let Some(n) = subset {
        if n > core_nodes.len() {
            panic!(
                "RuntimeError: Cannot subset core genes to {n}, only {} are available!",
                core_nodes.len()
            );
        }
        // `random.shuffle(core_nodes)` then `core_nodes[:subset]`. Reproducing which genes
        // survive needs CPython's Mersenne Twister and its exact Fisher-Yates order; see
        // PORTING_PLAN.md §6.7. Only reachable via --core_subset.
        panic!(
            "noimpl: --core_subset needs a CPython Mersenne Twister clone for \
             random.shuffle parity (PORTING_PLAN.md §6.7)"
        );
    }
    core_nodes
}

/// `generate_output.py::update_col_counts`
///
/// ```python
/// s = np.array(bytearray(s.lower().encode()), dtype=np.int8)
/// s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
/// col_counts[0,s==97] += 1   # a
/// col_counts[1,s==99] += 1   # c
/// col_counts[2,s==103] += 1  # g
/// col_counts[3,s==116] += 1  # t
/// col_counts[4,s==110] += 1  # n -- everything else, including gaps
/// ```
///
/// **Layout note.** The Python array is `5 x L` (row per base, column per alignment
/// position); here it is `L x 5` (one `[a,c,g,t,n]` count per position). That is a pure
/// transposition with no behavioural effect, and it keeps each column's counts contiguous
/// for [`calc_hc`]. Both this function and `calc_hc` use the transposed form consistently.
pub fn update_col_counts(col_counts: &mut [[i64; 5]], s: &str) {
    for (j, b) in s.bytes().enumerate() {
        let idx = match b.to_ascii_lowercase() {
            b'a' => 0,
            b'c' => 1,
            b'g' => 2,
            b't' => 3,
            _ => 4, // n, gaps, ambiguity codes -- all folded to 'n' by the mask above
        };
        col_counts[j][idx] += 1;
    }
}

/// `generate_output.py::calc_hc`
///
/// ```python
/// with np.errstate(divide='ignore', invalid='ignore'):
///     col_counts = col_counts/np.sum(col_counts,0)
///     hc = -np.nansum(col_counts[0:4,:]*np.log(col_counts[0:4,:]), 0)
/// return(np.sum((1-col_counts[4,:]) * hc)/np.sum(1-col_counts[4,:]))
/// ```
///
/// Mean per-column Shannon entropy, weighted by the non-`N` fraction. Used to drop
/// high-entropy columns from the core alignment.
///
/// Three details carry the float result:
///   - `0 * log(0)` is `NaN` in both numpy and Rust, and `nansum` drops it — that is how
///     absent bases contribute nothing rather than `-inf`.
///   - an all-`N` column has sum > 0, so no division by zero; an **empty** column divides
///     `0/0` to `NaN`, which `np.sum` (not `nansum`) then propagates to the result. That
///     poisoning is preserved deliberately.
///   - the outer sums run over every alignment column, so they go through
///     [`np_sum_f64`](crate::support::npmath::np_sum_f64) for numpy's pairwise order.
///     The inner `nansum` is over 4 elements, which numpy accumulates sequentially.
pub fn calc_hc(col_counts: &[[i64; 5]]) -> f64 {
    let n = col_counts.len();
    let mut numerator = Vec::with_capacity(n);
    let mut denominator = Vec::with_capacity(n);

    for col in col_counts {
        let total: f64 = col.iter().map(|&c| c as f64).sum::<f64>();
        // col_counts / np.sum(col_counts, 0) -- 0/0 yields NaN, as in numpy under errstate
        let p: [f64; 5] = std::array::from_fn(|i| col[i] as f64 / total);

        // hc = -nansum(p[0:4] * log(p[0:4]))
        let terms: [f64; 4] = std::array::from_fn(|i| p[i] * p[i].ln());
        let hc = -crate::support::npmath::np_nansum_f64(&terms);

        let w = 1.0 - p[4];
        numerator.push(w * hc);
        denominator.push(w);
    }

    crate::support::npmath::np_sum_f64(&numerator)
        / crate::support::npmath::np_sum_f64(&denominator)
}

/// `generate_output.py::concatenate_core_genome_alignments`
///
/// Writes `core_gene_alignment.aln` and `core_alignment_header.embl`; with `hc_threshold`
/// set, also the filtered variants.
pub fn concatenate_core_genome_alignments(
    core_names: &[String],
    output_dir: &str,
    hc_threshold: Option<f64>,
) {
    use crate::support::pydict::PyDict;
    use std::collections::BTreeSet;
    use std::io::Write;

    let alignments_dir = format!("{output_dir}/aligned_gene_sequences/");
    // Tier D: sorted(os.listdir(...)) -- raw listdir is filesystem order, which is not even
    // stable across machines, and it fixes the column order of core_gene_alignment.aln.
    let mut alignment_filenames: Vec<String> = std::fs::read_dir(&alignments_dir)
        .unwrap_or_else(|e| panic!("read_dir {alignments_dir}: {e}"))
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
        .collect();
    alignment_filenames.sort();

    let core_set: std::collections::HashSet<&str> = core_names.iter().map(|s| s.as_str()).collect();
    let core_filenames: Vec<String> = alignment_filenames
        .into_iter()
        .filter(|x| core_set.contains(x.split('.').next().unwrap_or("")))
        .collect();

    // Read in all these alignments.
    // Each entry: (gene_name, per-genome sequence, gene_length, entropy)
    let mut gene_alignments: Vec<GeneAlignment> = Vec::new();
    let mut isolates: BTreeSet<String> = BTreeSet::new(); // Tier D: sorted(isolates)

    for filename in &core_filenames {
        let gene_name = {
            let base = filename.rsplit('/').next().unwrap_or(filename);
            match base.rfind('.') {
                Some(i) if i > 0 => base[..i].to_string(),
                _ => base.to_string(),
            }
        };
        let alignment =
            crate::support::align_io::read_fasta_alignment(&format!("{alignments_dir}{filename}"));
        let mut gene_dict: PyDict<String, (String, String)> = PyDict::new();
        let mut gene_length = 0usize;
        let mut col_counts: Vec<[i64; 5]> = Vec::new();

        for record in &alignment.records {
            if gene_dict.is_empty() && col_counts.is_empty() {
                gene_length = record.seq.len();
                col_counts = vec![[0i64; 5]; gene_length];
            }
            update_col_counts(&mut col_counts, &record.seq);

            let mut rid = record.id.clone();
            if rid.starts_with("_R_") {
                rid = rid[3..].to_string();
            }
            let genome_id = rid.split(';').next().unwrap().to_string();

            if let Some(existing) = gene_dict.get(&genome_id) {
                if record.seq.matches('-').count() < existing.1.matches('-').count() {
                    gene_dict.insert(genome_id.clone(), (rid.clone(), record.seq.clone()));
                }
            } else {
                gene_dict.insert(genome_id.clone(), (rid.clone(), record.seq.clone()));
            }
            isolates.insert(genome_id);
        }
        let hc = calc_hc(&col_counts);
        gene_alignments.push((gene_name, gene_dict, gene_length, hc));
    }

    let total_alignment_len: usize = gene_alignments.iter().map(|g| g.2).sum();
    write_core_alignment_records(
        &format!("{output_dir}core_gene_alignment.aln"),
        &isolates,
        &gene_alignments,
        total_alignment_len,
        None,
    );
    let header_list: Vec<(String, usize)> =
        gene_alignments.iter().map(|g| (g.0.clone(), g.2)).collect();
    crate::generate_alignments::write_alignment_header(
        &header_list,
        output_dir,
        "core_alignment_header.embl",
    );

    // Calculate threshold for h.
    let hc_threshold = match hc_threshold {
        Some(t) => t,
        None => {
            let allh: Vec<f64> = gene_alignments.iter().map(|g| g.3).collect();
            let q = np_quantile(&allh, &[0.25, 0.75]);
            let t = 0.01f64.max(q[1] + 1.5 * (q[1] - q[0]));
            println!(
                "Entropy threshold automatically set to {}.",
                crate::support::pyfmt::py_str_f64(t)
            );
            t
        }
    };

    {
        let f = std::fs::File::create(format!("{output_dir}alignment_entropy.csv")).unwrap();
        let mut w = std::io::BufWriter::new(f);
        for g in &gene_alignments {
            writeln!(w, "{},{}", g.0, crate::support::pyfmt::py_str_f64(g.3)).unwrap();
        }
    }

    let filtered: Vec<(String, usize)> = gene_alignments
        .iter()
        .filter(|g| g.3 <= hc_threshold)
        .map(|g| (g.0.clone(), g.2))
        .collect();
    let filtered_alignment_len: usize = filtered.iter().map(|g| g.1).sum();
    let keep_count = write_core_alignment_records(
        &format!("{output_dir}core_gene_alignment_filtered.aln"),
        &isolates,
        &gene_alignments,
        filtered_alignment_len,
        Some(hc_threshold),
    );
    crate::generate_alignments::write_alignment_header(
        &filtered,
        output_dir,
        "core_alignment_filtered_header.embl",
    );

    println!(
        "{} out of {} genes kept in filtered core genome",
        crate::support::pyfmt::py_str_f64(keep_count as f64 / isolates.len() as f64),
        gene_alignments.len()
    );
}

/// `np.quantile(a, qs)` with numpy's default `linear` interpolation.
///
/// Reproduces NumPy 1.26's behaviour: sort, then for each `q` take
/// `virtual_index = q * (n - 1)` and linearly interpolate between the two neighbouring
/// order statistics. See `NOTICE.md` for NumPy attribution.
fn np_quantile(a: &[f64], qs: &[f64]) -> Vec<f64> {
    let mut v = a.to_vec();
    v.sort_by(|x, y| x.partial_cmp(y).unwrap());
    qs.iter()
        .map(|&q| {
            if v.is_empty() {
                return f64::NAN;
            }
            let idx = q * (v.len() - 1) as f64;
            let lo = idx.floor() as usize;
            let hi = idx.ceil() as usize;
            if lo == hi {
                v[lo]
            } else {
                v[lo] + (idx - lo as f64) * (v[hi] - v[lo])
            }
        })
        .collect()
}

/// `generate_output.py::generate_core_genome_alignment`
#[allow(clippy::too_many_arguments)]
pub fn generate_core_genome_alignment(
    g: &Graph,
    temp_dir: &str,
    output_dir: &str,
    threads: i64,
    aligner: &str,
    isolates: &[String],
    threshold: f64,
    codons: bool,
    strict: bool,
    num_isolates: usize,
    hc_threshold: Option<f64>,
    subset: Option<usize>,
    resume: bool,
) {
    use crate::generate_alignments::*;

    let _ = std::fs::create_dir(format!("{output_dir}aligned_gene_sequences"));

    let core_genes = get_core_gene_nodes(g, threshold, num_isolates, subset);
    let core_gene_names: Vec<String> = core_genes
        .iter()
        .map(|&x| g.node(x).name.clone().expect("node name set"))
        .collect();

    if core_genes.is_empty() {
        println!(
            "No gene clusters were present above the core frequency threshold! \
             Try adjusting the '--core_threshold' parameter"
        );
        return;
    }

    let nodes: Vec<(usize, &crate::support::graph::NodeAttrs)> =
        core_genes.iter().map(|&id| (id, g.node(id))).collect();
    let pending_gene_ids = get_pending_gene_ids(&nodes, output_dir, codons, resume);
    let total_gene_count = core_genes.len();

    if codons || strict {
        let (protein_pending, reverse_translate_pending) =
            get_pending_codon_gene_ids(&nodes, output_dir, resume);
        println!("Codon alignment is experimental in Panaroo...");
        let _ = std::fs::create_dir(format!("{output_dir}aligned_protein_sequences"));
        let _ = std::fs::create_dir(format!("{output_dir}unaligned_dna_sequences"));

        let proteins_dic = fasta_by_id(&format!("{output_dir}combined_protein_CDS.fasta"));
        let nucleotides_dic = fasta_by_id(&format!("{output_dir}combined_DNA_CDS.fasta"));
        let clean_isolates: Vec<String> = isolates.iter().map(|s| s.replace(';', "")).collect();

        let mut output_files = Vec::new();
        for gene in &protein_pending {
            output_files.push(output_dna_and_protein(
                g.node(*gene),
                &clean_isolates,
                temp_dir,
                output_dir,
                &proteins_dic,
                &nucleotides_dic,
            ));
        }
        let unaligned_protein_files: Vec<String> =
            output_files.into_iter().filter_map(|(p, _)| p).collect();

        let commands: Vec<AlignCommand> = unaligned_protein_files
            .iter()
            .map(|f| get_protein_commands(Some(f), output_dir, aligner, threads))
            .collect();
        print_stage_progress(
            "Protein alignments",
            total_gene_count - protein_pending.len(),
            commands.len(),
            Some(total_gene_count),
        );
        multi_align_sequences(
            commands,
            &format!("{output_dir}aligned_protein_sequences/"),
            threads,
            aligner,
        );

        let (protein_sequences, unaligned_dna_files) =
            get_codon_pending_files(&nodes, output_dir, &reverse_translate_pending);
        for file in &protein_sequences {
            if !std::path::Path::new(file).is_file() {
                println!("{file}");
                panic!("RuntimeError: Some alignments failed to complete!");
            }
        }
        if !reverse_translate_pending.is_empty() {
            let completed = total_gene_count - reverse_translate_pending.len();
            reverse_translate_sequences(
                &protein_sequences,
                &unaligned_dna_files,
                strict,
                output_dir,
                temp_dir,
                aligner,
                threads,
                completed,
            );
        }
    } else {
        let mut temp_dir = temp_dir.to_string();
        if aligner == "none" {
            temp_dir = format!("{output_dir}unaligned_gene_sequences/");
            let _ = std::fs::create_dir(&temp_dir);
        }

        print_stage_progress(
            "Gene alignments",
            total_gene_count - pending_gene_ids.len(),
            pending_gene_ids.len(),
            Some(total_gene_count),
        );
        let td = temp_dir.clone();
        let od = output_dir.to_string();
        // A borrow, not a clone -- see the matching comment in generate_pan_genome_alignment.
        let jobs: Vec<&crate::support::graph::NodeAttrs> =
            pending_gene_ids.iter().map(|&x| g.node(x)).collect();
        // PORTING_PLAN.md §9 item 1: one parse of combined_DNA_CDS.fasta for the stage.
        let unaligned: Vec<Option<String>> = if jobs.is_empty() {
            Vec::new()
        } else {
            let combined_dna = CombinedDna::load(output_dir, isolates);
            crate::support::parallel::parallel_map(threads, jobs, move |node| {
                output_sequence(node, &combined_dna, &td, &od)
            })
        };

        if aligner == "none" {
            println!("No aligner specified. Returning unaligned gene fasta files.");
            return;
        }

        let unaligned: Vec<String> = unaligned.into_iter().flatten().collect();
        let commands: Vec<AlignCommand> = unaligned
            .iter()
            .map(|f| get_alignment_commands(f, output_dir, aligner, threads))
            .collect();
        multi_align_sequences(
            commands,
            &format!("{output_dir}aligned_gene_sequences/"),
            threads,
            aligner,
        );
    }

    // Concatenate them together to produce the two output files
    concatenate_core_genome_alignments(&core_gene_names, output_dir, hc_threshold);
}

/// `generate_output.py::generate_summary_stats`
///
/// Reads `gene_presence_absence.Rtab` back off disk and writes `summary_statistics.txt`.
pub fn generate_summary_stats(output_dir: &str) {
    let text = std::fs::read_to_string(format!("{output_dir}gene_presence_absence_roary.csv"))
        .expect("read gene_presence_absence_roary.csv");
    let mut lines = text.lines();
    let header = lines
        .next()
        .expect("gene_presence_absence_roary.csv header");
    let no_samples = header.split(',').count() - 14;

    let (mut no_core, mut no_soft_core, mut no_shell, mut no_cloud, mut total) = (0, 0, 0, 0, 0);
    for gene in lines {
        let no_isolates = gene
            .split(',')
            .nth(3)
            .expect("gene_presence_absence_roary.csv No. isolates column");
        let proportion_present = no_isolates.parse::<f64>().unwrap() / no_samples as f64 * 100.0;
        if proportion_present >= 99.0 {
            no_core += 1;
        } else if proportion_present >= 95.0 {
            no_soft_core += 1;
        } else if proportion_present >= 15.0 {
            no_shell += 1;
        } else {
            no_cloud += 1;
        }
        total += 1;
    }

    // Note: no trailing newline -- the Python writes the joined string as-is.
    let output = format!(
        "Core genes\t(99% <= strains <= 100%)\t{no_core}\n\
         Soft core genes\t(95% <= strains < 99%)\t{no_soft_core}\n\
         Shell genes\t(15% <= strains < 95%)\t{no_shell}\n\
         Cloud genes\t(0% <= strains < 15%)\t{no_cloud}\n\
         Total genes\t(0% <= strains <= 100%)\t{total}"
    );
    std::fs::write(format!("{output_dir}summary_statistics.txt"), output).expect("write summary");
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values produced by the reference Python.

    #[test]
    fn seq_sample_key_matches_python_split_join_prefix() {
        assert_eq!(seq_sample_key("sample_gene_1"), "sample");
        assert_eq!(seq_sample_key("sample_name_gene_1"), "sample_name");
        assert_eq!(seq_sample_key("_gene_1"), "");
        assert_eq!(seq_sample_key("sample__1"), "sample");
        assert_eq!(seq_sample_key("sample_1"), "");
    }

    #[test]
    fn roary_gene_name_matches_join_then_filter() {
        assert_eq!(
            roary_gene_name(";;geneA;gene-B; ;x_y~z;;"),
            "geneA~~~geneB~~~~~~x_y~z"
        );
        assert_eq!(
            roary_gene_name(" ; weird/name ;two.words; "),
            "weirdname~~~twowords"
        );
        assert_eq!(roary_gene_name(";;;"), "");
    }

    #[test]
    fn length_stats_match_existing_numpy_helpers() {
        use crate::support::npmath::{np_max_i64, np_mean_f64, np_min_i64};
        use crate::support::pyfmt::py_str_f64;

        let lengths = [3usize, 10, 4, 10, 5];
        let mut counts = Vec::new();
        let (min, max, mean, mode) = length_stats_and_first_mode(&lengths, &mut counts);
        let li: Vec<i64> = lengths.iter().map(|&x| x as i64).collect();
        let lf: Vec<f64> = lengths.iter().map(|&x| x as f64).collect();
        assert_eq!(min.to_string(), np_min_i64(&li).to_string());
        assert_eq!(max.to_string(), np_max_i64(&li).to_string());
        assert_eq!(py_str_f64(mean), py_str_f64(np_mean_f64(&lf)));
        assert_eq!(mode, 10);
    }

    #[test]
    fn first_modal_length_matches_python_first_max() {
        let mut counts = Vec::new();
        assert_eq!(length_stats_and_first_mode(&[9, 5, 9, 5], &mut counts).3, 9);
        assert_eq!(
            length_stats_and_first_mode(&[4, 7, 7, 4, 7, 4], &mut counts).3,
            4
        );
        assert_eq!(length_stats_and_first_mode(&[8], &mut counts).3, 8);
    }

    #[test]
    fn summary_stats_uses_roary_counts() {
        let base =
            std::env::temp_dir().join(format!("panaroo-summary-test-{}", std::process::id()));
        std::fs::create_dir_all(&base).unwrap();
        let dir_path = crate::support::pytempfile::mkdtemp(&base).unwrap();
        let dir = format!("{}/", dir_path.display());
        std::fs::write(
            format!("{dir}gene_presence_absence_roary.csv"),
            concat!(
                "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc,s1,s2,s3,s4\n",
                "g1,,,4,,,,,,,,,,,a,b,c,d\n",
                "g2,,,3,,,,,,,,,,,a,b,c,\n",
                "g3,,,2,,,,,,,,,,,a,b,,\n",
                "g4,,,0,,,,,,,,,,,,,,\n",
            ),
        )
        .unwrap();

        generate_summary_stats(&dir);
        let summary = std::fs::read_to_string(format!("{dir}summary_statistics.txt")).unwrap();
        assert_eq!(
            summary,
            "Core genes\t(99% <= strains <= 100%)\t1\n\
             Soft core genes\t(95% <= strains < 99%)\t0\n\
             Shell genes\t(15% <= strains < 95%)\t2\n\
             Cloud genes\t(0% <= strains < 15%)\t1\n\
             Total genes\t(0% <= strains <= 100%)\t4"
        );
        std::fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn update_col_counts_matches_numpy_masking() {
        // python, after feeding "ACGTNN", "ACGTAA", "aCgTNa" into a zeros((5,6)) array:
        //   [[3 0 0 0 1 2]
        //    [0 3 0 0 0 0]
        //    [0 0 3 0 0 0]
        //    [0 0 0 3 0 0]
        //    [0 0 0 0 2 1]]
        // transposed to our L x 5 layout:
        let mut cc = [[0i64; 5]; 6];
        for s in ["ACGTNN", "ACGTAA", "aCgTNa"] {
            update_col_counts(&mut cc, s);
        }
        assert_eq!(cc[0], [3, 0, 0, 0, 0]);
        assert_eq!(cc[1], [0, 3, 0, 0, 0]);
        assert_eq!(cc[2], [0, 0, 3, 0, 0]);
        assert_eq!(cc[3], [0, 0, 0, 3, 0]);
        assert_eq!(cc[4], [1, 0, 0, 0, 2]);
        assert_eq!(cc[5], [2, 0, 0, 0, 1]);
    }

    #[test]
    fn calc_hc_matches_python() {
        // python: calc_hc(cc.astype(float)) -> 0.06045494935779483
        let mut cc = [[0i64; 5]; 6];
        for s in ["ACGTNN", "ACGTAA", "aCgTNa"] {
            update_col_counts(&mut cc, s);
        }
        assert_eq!(calc_hc(&cc), 0.06045494935779483);

        // python: cc2 from "AAA" and "ACA" -> calc_hc -> 0.23104906018664842
        let mut cc2 = [[0i64; 5]; 3];
        update_col_counts(&mut cc2, "AAA");
        update_col_counts(&mut cc2, "ACA");
        assert_eq!(calc_hc(&cc2), 0.23104906018664842);
    }
}
