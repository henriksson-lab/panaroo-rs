//! Translation of `panaroo/panaroo/prokka.py` — GFF3 input parsing and the first pass over
//! every genome.

use crate::support::gff::GffEntry;
use crate::support::pydict::PyDict;
use crate::support::seqio::SeqRecord;
use std::io::Write;

/// `prokka.py::bact_translation_table`
///
/// Panaroo's own hard-coded bacterial code (NCBI table 11), as a `5 x 5 x 5` array of
/// amino-acid bytes indexed `[base1][base2][base3]` with `A=0 C=1 G=2 T=3, anything else=4`.
/// Index 4 in any position yields `X`.
///
/// Transcribed from the literal in `prokka.py`. This one *is* Panaroo's, so it stays here;
/// the NCBI tables it is patched against live in [`crate::support::codon_table`].
#[rustfmt::skip]
pub const BACT_TRANSLATION_TABLE: [[[u8; 5]; 5]; 5] = [
    [ *b"KNKNX", *b"TTTTT", *b"RSRSX", *b"IIMIX", *b"XXXXX" ],
    [ *b"QHQHX", *b"PPPPP", *b"RRRRR", *b"LLLLL", *b"XXXXX" ],
    [ *b"EDEDX", *b"AAAAA", *b"GGGGG", *b"VVVVV", *b"XXXXX" ],
    [ *b"*Y*YX", *b"SSSSS", *b"*CWCX", *b"LFLFX", *b"XXXXX" ],
    [ *b"XXXXX", *b"XXXXX", *b"XXXXX", *b"XXXXX", *b"XXXXX" ],
];

/// `prokka.py::reduce_array`
///
/// ```python
/// reduce_array = np.full(200, 4)
/// reduce_array[[65, 97]] = 0   # A a
/// reduce_array[[67, 99]] = 1   # C c
/// reduce_array[[71, 103]] = 2  # G g
/// reduce_array[[84, 116]] = 3  # T t
/// ```
///
/// The Python indexes this with `np.int8`, so a byte above 127 becomes negative and wraps
/// to the *end* of the 200-element array — where every entry is still 4. So every
/// non-`ACGTacgt` byte reduces to 4 regardless, and [`reduce`] can just match.
pub fn reduce(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// The `[translation_table, set(start_codons)]` pair returned by `get_trans_table`.
///
/// The Python returns a two-element list of mixed type; a struct is the honest Rust
/// equivalent and does not change behaviour.
#[derive(Debug, Clone)]
pub struct TransTable {
    /// `translation_table[0]` — the `5 x 5 x 5` codon lookup.
    pub table: [[[u8; 5]; 5]; 5],
    /// `translation_table[1]` — `set(tb.start_codons)`, membership-tested only.
    pub start_codons: Vec<String>,
}

/// `prokka.py::get_trans_table`
///
/// Duplicated verbatim as [`crate::generate_alignments::get_trans_table`]. Both are kept;
/// see PORTING_PLAN.md §5.
pub fn get_trans_table(table: i64) -> TransTable {
    // swap to different codon table
    let mut translation_table = BACT_TRANSLATION_TABLE;
    let tb = crate::support::codon_table::generic_by_id(table).unwrap_or_else(|| {
        panic!("Invalid codon table! Must be available as a generic table in BioPython")
    });

    if table != 11 {
        for &(codon, aa) in tb.forward_table {
            if codon.contains('U') {
                continue;
            }
            let b = codon.as_bytes();
            translation_table[reduce(b[0])][reduce(b[1])][reduce(b[2])] = aa as u8;
        }
        for &codon in tb.stop_codons {
            if codon.contains('U') {
                continue;
            }
            let b = codon.as_bytes();
            translation_table[reduce(b[0])][reduce(b[1])][reduce(b[2])] = b'*';
        }
    }

    TransTable {
        table: translation_table,
        start_codons: tb.start_codons.iter().map(|s| s.to_string()).collect(),
    }
}

/// `prokka.py::translate`
///
/// Panaroo's own numpy lookup-table translator — **not** `Bio.Seq.translate`. Substitutes
/// `M` for the first residue when the sequence starts with a recognised start codon.
///
/// Duplicated verbatim as [`crate::generate_alignments::translate`]. Both are kept.
pub fn translate(seq: &str, translation_table: &TransTable) -> String {
    let b = seq.as_bytes();
    // The Python builds three index arrays with np.arange(0|1|2, len(seq), 3) and fancy-
    // indexes with all three at once, which requires len(seq) % 3 == 0 -- otherwise numpy
    // raises on the shape mismatch. Callers guarantee it.
    assert!(
        b.len() % 3 == 0,
        "translate: sequence length {} is not divisible by 3",
        b.len()
    );

    let mut pseq = String::with_capacity(b.len() / 3);
    for codon in b.chunks_exact(3) {
        let aa = translation_table.table[reduce(codon[0])][reduce(codon[1])][reduce(codon[2])];
        pseq.push(aa as char);
    }

    // Check for a different start codon.
    if b.len() >= 3
        && translation_table
            .start_codons
            .iter()
            .any(|c| c.as_bytes() == &b[0..3])
    {
        let mut out = String::with_capacity(pseq.len());
        out.push('M');
        out.push_str(&pseq[1..]);
        return out;
    }
    pseq
}

/// `prokka.py::create_temp_gff3`
///
/// Either converts GenBank via [`crate::biocode_convert::convert_gbk_gff3`], or splices a
/// separate FASTA onto a GFF3 under a `##FASTA` line. Also strips the PATRIC `accn|`
/// prefix.
pub fn create_temp_gff3(gff_file: &str, fasta_file: Option<&str>, temp_dir: &str) -> String {
    // create directory if it isn't present already
    let dir = format!("{temp_dir}temp_gffs");
    if !std::path::Path::new(&dir).exists() {
        std::fs::create_dir(&dir).unwrap_or_else(|e| panic!("mkdir {dir}: {e}"));
    }

    let prefix = file_stem(gff_file);
    let ext = extension(gff_file);
    let out_path = format!("{temp_dir}temp_gffs/{prefix}.gff");

    match fasta_file {
        None => {
            // GenBank input. Upstream wraps the conversion in a bare `try/except` that
            // swallows the ORIGINAL exception and re-raises a generic RuntimeError naming
            // the offending file, after printing a hint to stdout. Reproduced because the
            // printed hint is the only guidance a user gets on a non-compliant GenBank, and
            // because the wrapper is what turns any converter failure into a message that
            // names the file.
            let converted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                crate::biocode_convert::convert_gbk_gff3(gff_file, &out_path, true)
            }));
            if converted.is_err() {
                println!(
                    "Error reading Genbank input! These must compliant with \
                     Genbank/ENA/DDJB. This can be forced in Prokka by specifying the \
                     --compliance parameter."
                );
                panic!(
                    "RuntimeError: Error reading Genbank input: {gff_file}\n\
                     These must compliant with Genbank/ENA/DDJB. This can be forced in \
                     Prokka by specifying the --compliance parameter."
                );
            }
        }
        Some(fasta) => {
            if ext != ".gff" && ext != ".gff3" {
                panic!("RuntimeError: Invalid file extension! ({ext})");
            }
            let fasta_ext = extension(fasta);
            if !matches!(fasta_ext.as_str(), ".fasta" | ".fa" | ".fas" | ".fna") {
                panic!("RuntimeError: Invalid file extension! ({fasta_ext})");
            }

            // merge files into temporary gff3
            let gff_text = std::fs::read_to_string(gff_file)
                .unwrap_or_else(|e| panic!("could not read {gff_file}: {e}"));
            let mut gff_string = gff_text.trim().to_string();
            if gff_string.contains("\naccn") {
                // deal with PATRIC input format
                gff_string = gff_string.replace("accn|", "");
            }
            let fasta_text = std::fs::read_to_string(fasta)
                .unwrap_or_else(|e| panic!("could not read {fasta}: {e}"));

            let mut out = String::with_capacity(gff_string.len() + fasta_text.len() + 16);
            out.push_str(&gff_string);
            out.push_str("\n##FASTA\n");
            out.push_str(fasta_text.trim());
            std::fs::write(&out_path, out).unwrap_or_else(|e| panic!("write {out_path}: {e}"));
        }
    }

    out_path
}

/// `os.path.splitext(p)[1]` — the extension including its leading dot, or `""`.
fn extension(p: &str) -> String {
    let base = p.rsplit('/').next().unwrap_or(p);
    match base.rfind('.') {
        Some(i) if i > 0 => base[i..].to_string(),
        _ => String::new(),
    }
}

/// `prokka.py::clean_gff_string`
///
/// ```python
/// splitlines = gff_string.splitlines()
/// lines_to_delete = []
/// for index in range(len(splitlines)):
///     if '##sequence-region' in splitlines[index]:
///         lines_to_delete.append(index)
/// for index in sorted(lines_to_delete, reverse=True):
///     del splitlines[index]
/// cleaned_gff = "\n".join(splitlines)
/// return cleaned_gff
/// ```
///
/// The delete-by-descending-index dance is just a filter; the observable result is the
/// surviving lines joined by `\n`. Note `splitlines()` discards the final line terminator,
/// so a trailing newline in the input is **not** present in the output — verified against
/// Python.
///
/// The test is `in`, not `startswith`, so a line merely *containing* the directive
/// anywhere is dropped too. Preserved.
pub fn clean_gff_string(gff_string: &str) -> String {
    gff_string
        .lines()
        .filter(|l| !l.contains("##sequence-region"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// `prokka.py::get_gene_sequences`
///
/// Returns `(sequence_dictionary, protein_list)`. The keys of `sequence_dictionary` are the
/// clustering IDs `{file_number}_{scaffold_index}_{gene_index}` that every downstream file
/// is built on — see PORTING_PLAN.md §6.6 for why the feature and scaffold ordering has to
/// match gffutils exactly.
///
/// PORTING_PLAN.md §9 item 3: the inner loop scans every contig for every CDS and does not
/// `break` on a match. Preserve it — including the missing `break`, which is what makes the
/// `gene_sequence is None` check below the loop reachable.
///
/// // UPSTREAM BUG: without the `break`, a duplicated contig ID would emit the gene twice.
pub fn get_gene_sequences(
    gff_file_name: &str,
    file_number: usize,
    filter_seqs: bool,
    table: &TransTable,
) -> (PyDict<String, SeqRecord>, Vec<SeqRecord>) {
    // Get name and separate the prokka GFF into separate GFF and FASTA files
    if gff_file_name.contains(',') {
        println!("Problem reading GFF3 file:  {gff_file_name}");
        panic!("RuntimeError: Error reading prokka input!");
    }

    let raw = std::fs::read_to_string(gff_file_name)
        .unwrap_or_else(|e| panic!("could not read {gff_file_name}: {e}"));
    let mut sequence_dictionary: PyDict<String, SeqRecord> = PyDict::new();

    // Split file and parse. Note every comma is stripped from the *whole* file first --
    // that is what makes the unquoted CSV in `output_files` safe, and it happens before
    // gffutils percent-decodes, so a `%2C` still becomes a comma later.
    let lines = raw.replace(',', "");
    let split: Vec<&str> = lines.split("##FASTA").collect();
    if split.len() != 2 {
        println!("Problem reading GFF3 file:  {gff_file_name}");
        panic!("RuntimeError: Error reading prokka input!");
    }

    let sequences = crate::support::seqio::parse_fasta(split[1]);
    let parsed_gff = crate::support::gff::GffDb::create_db(&clean_gff_string(split[0]))
        .unwrap_or_else(|e| panic!("RuntimeError: Error reading prokka input! {e}"));

    // Get genes per scaffold
    let mut scaffold_genes: PyDict<String, Vec<(i64, SeqRecord)>> = PyDict::new();
    for entry in parsed_gff.all_features() {
        if !entry.featuretype.contains("CDS") {
            continue;
        }

        let mut gene_sequence_found = false;
        // PORTING_PLAN.md §9 item 3: linear scan of every contig, with no `break` after a
        // match. Preserved -- the missing `break` is what makes the `gene_sequence is None`
        // check below reachable at all. The scan itself is cheap (a few short string
        // compares per contig); only clone the contig id once the ids actually match.
        for sequence_index in 0..sequences.len() {
            if sequences[sequence_index].id == entry.seqid {
                let scaffold_id = sequences[sequence_index].id.clone();
                let full = &sequences[sequence_index].seq;
                let lo = (entry.start - 1).max(0) as usize;
                let hi = (entry.stop as usize).min(full.len());
                let mut seq = if lo < hi {
                    full[lo..hi].to_string()
                } else {
                    String::new()
                };
                if entry.strand == "-" {
                    seq = crate::support::seq::reverse_complement(&seq);
                }
                gene_sequence_found = true;

                let mut gene_name = attr_first_or_empty(entry, "gene");
                if gene_name.is_empty() {
                    gene_name = attr_first_or_empty(entry, "name");
                }

                let gene_description = match entry.attribute("product") {
                    Some(v) => v.join(";").replace(',', ""),
                    None => String::new(),
                };

                // clean entries if requested
                if entry.frame != "0" {
                    println!("Invalid gene! Panaroo currently does not support frame shifts.");
                    if filter_seqs {
                        continue;
                    } else {
                        panic!("ValueError: Invalid gene sequence!");
                    }
                }

                // The Python's `or` short-circuits, so the stop-codon test is only reached
                // when `len(gene_sequence) % 3 == 0` -- which is exactly the condition under
                // which `pad3` was a no-op. Kept lazy for the same reason: `translate`
                // asserts `len % 3 == 0`, and the first operand is what guarantees it.
                let mut bad = (seq.len() % 3 > 0) || (seq.len() < 34);
                if !bad {
                    let prot = translate(&seq, table);
                    bad = prot
                        .get(..prot.len().saturating_sub(1))
                        .unwrap_or("")
                        .contains('*');
                }
                if bad {
                    println!("invalid gene! file - id:  {gff_file_name}  -  {}", entry.id);
                    println!("Length: {} , Has stop: ...", seq.len());
                    if filter_seqs {
                        continue;
                    } else {
                        panic!("ValueError: Invalid gene sequence!");
                    }
                }

                let mut rec = SeqRecord {
                    id: entry.id.clone(),
                    name: gene_name,
                    description: gene_description,
                    seq,
                    annotations: Vec::new(),
                };
                rec.set_annotation("scaffold", &scaffold_id);

                if !scaffold_genes.contains_key(&scaffold_id) {
                    scaffold_genes.insert(scaffold_id.clone(), Vec::new());
                }
                scaffold_genes
                    .get_mut(&scaffold_id)
                    .unwrap()
                    .push((entry.start, rec));
            }
        }
        if !gene_sequence_found {
            println!("Sequence ID not found in Fasta! {}", entry.seqid);
            if filter_seqs {
                continue;
            } else {
                panic!("ValueError: Invalid gene sequence!");
            }
        }
    }

    if scaffold_genes.is_empty() {
        println!("No valid sequences found in GFF! {gff_file_name}");
        panic!("ValueError: Invalid GFF!");
    }

    // `sorted(..., key=lambda x: x[0])` is stable, so genes starting at the same coordinate
    // keep GFF order.
    for v in scaffold_genes.values_mut() {
        v.sort_by_key(|(start, _)| *start);
    }

    // scaff_count follows insertion order into scaffold_genes, i.e. the order contigs are
    // first seen while walking features -- which is why gffutils' feature order matters.
    // See PORTING_PLAN.md §6.6.
    let scaffolds: Vec<String> = scaffold_genes.keys().cloned().collect();
    for (scaff_count, scaffold) in scaffolds.iter().enumerate() {
        // `scaffold_genes` is dead after this loop and each key is visited exactly once
        // (the keys came from `PyDict::keys()`), so take the vector rather than deep-cloning
        // every SeqRecord in it. `mem::take` leaves an empty Vec behind, so no key is
        // removed and IndexMap order is untouched.
        let genes = std::mem::take(scaffold_genes.get_mut(scaffold).unwrap());
        for (gene_index, (_, rec)) in genes.into_iter().enumerate() {
            let clustering_id = format!("{file_number}_{scaff_count}_{gene_index}");
            sequence_dictionary.insert(clustering_id, rec);
        }
    }

    let proteins = translate_sequences(&sequence_dictionary, table);
    (sequence_dictionary, proteins)
}

/// `prokka.py::translate_sequences`
///
/// Raises on a premature stop codon or a length not divisible by three.
pub fn translate_sequences(
    sequence_dic: &PyDict<String, SeqRecord>,
    table: &TransTable,
) -> Vec<SeqRecord> {
    let mut protein_list = Vec::new();
    for (strain_id, sequence_record) in sequence_dic.items() {
        if sequence_record.seq.len() % 3 != 0 {
            panic!("ValueError: Coding sequence not divisible by 3, is it complete?!");
        }
        let mut protien_sequence = translate(&sequence_record.seq, table);
        if protien_sequence.ends_with('*') {
            protien_sequence.pop();
        }
        if protien_sequence.contains('*') {
            panic!("ValueError: Premature stop codon in a gene!");
        }
        protein_list.push(SeqRecord::new(
            protien_sequence,
            strain_id.clone(),
            strain_id.clone(),
        ));
    }
    protein_list
}

/// `prokka.py::output_files`
///
/// Appends to `combined_protein_CDS.fasta`, `combined_DNA_CDS.fasta` and `gene_data.csv`.
/// The CSV is written by hand with `",".join(...)` and no quoting — which is safe only
/// because `get_gene_sequences` strips every comma from the input up front.
pub fn output_files<P: Write, D: Write, C: Write>(
    dna_dictionary: &PyDict<String, SeqRecord>,
    protien_list: &[SeqRecord],
    prot_handle: &mut P,
    dna_handle: &mut D,
    csv_handle: &mut C,
    gff_filename: &str,
) {
    use crate::support::seqio::write_fasta;

    // Simple output for protien list
    write_fasta(protien_list, prot_handle);

    // Correct DNA ids to CD-Hit acceptable ids, and output
    let clustering_id_records: Vec<SeqRecord> = dna_dictionary
        .items()
        .map(|(clusteringid, rec)| {
            SeqRecord::new(rec.seq.clone(), clusteringid.clone(), clusteringid.clone())
        })
        .collect();
    write_fasta(&clustering_id_records, dna_handle);

    // `os.path.splitext(os.path.basename(gff_filename))[0]`
    let gff_name = file_stem(gff_filename);

    // Combine everything to a csv and output it.
    //
    // Written by hand with `",".join(...)` and no quoting. That is only safe because
    // `get_gene_sequences` strips every comma from the input file up front -- but note the
    // *description* can still contain a comma, because gffutils percent-decodes `%2C`
    // after the strip. Panaroo re-strips commas from the description for that reason; the
    // gene name is not re-stripped, which is upstream behaviour.
    for protien in protien_list {
        let clustering_id = &protien.id;
        let relevant_seqrecord = dna_dictionary
            .get(clustering_id)
            .unwrap_or_else(|| panic!("KeyError: {clustering_id}"));
        // `",".join(out_list)` -- written field by field so no ~1.4 KB temporary is built
        // per row. Same bytes, same order.
        writeln!(
            csv_handle,
            "{},{},{},{},{},{},{},{}",
            gff_name,
            relevant_seqrecord.scaffold(),
            clustering_id,
            relevant_seqrecord.id,
            protien.seq,
            relevant_seqrecord.seq,
            relevant_seqrecord.name,
            relevant_seqrecord.description,
        )
        .expect("write gene_data.csv");
    }
}

/// `os.path.splitext(os.path.basename(p))[0]`
///
/// Not a Python function — Rust's `Path::file_stem` differs on names like `a.b.c` (it strips
/// only the last extension, which matches) and on dotfiles (`.bashrc`: Python's splitext
/// gives `.bashrc`, Rust's file_stem gives `.bashrc` too). Spelled out here so the
/// behaviour is pinned rather than inherited.
fn file_stem(p: &str) -> String {
    let base = p.rsplit('/').next().unwrap_or(p);
    match base.rfind('.') {
        // splitext does not treat a leading dot as an extension separator
        Some(i) if i > 0 => base[..i].to_string(),
        _ => base.to_string(),
    }
}

/// `prokka.py::process_prokka_input`
///
/// PORTING_PLAN.md §9 item 9: the Python chunks the file list into batches of `n_cpu` and
/// spawns a fresh joblib pool per batch. Results are consumed in batch order, so a single
/// pool produces identical output.
///
/// The batching is not *only* a pool-lifetime detail, though: it also caps how many results
/// are alive at once. Each result here is one genome's `SeqRecord`s (~9-10 MB), so collecting
/// all `N` before writing any -- which is what an earlier version of this function did --
/// holds ~5 GB at `N = 500` where CPython holds ~`n_cpu` genomes' worth. See
/// [`crate::support::parallel::parallel_map_ordered_consume`], which keeps the single pool
/// while restoring CPython's memory bound and its in-order consumption.
pub fn process_prokka_input(
    gff_list: &[String],
    output_dir: &str,
    filter_seqs: bool,
    _quiet: bool,
    n_cpu: i64,
    table: i64,
) -> bool {
    use std::io::BufWriter;

    let trans_table = get_trans_table(table);

    // 1 MiB rather than BufWriter's 8 KiB default: phase 1 writes ~55 MB across these three
    // handles on the `ci` dataset, which is ~6700 write(2) calls at the default capacity.
    // Buffer capacity cannot change the byte stream or its order.
    const WRITE_BUF: usize = 1 << 20;
    let mut prot_handle = BufWriter::with_capacity(
        WRITE_BUF,
        std::fs::File::create(format!("{output_dir}combined_protein_CDS.fasta"))
            .expect("create combined_protein_CDS.fasta"),
    );
    let mut dna_handle = BufWriter::with_capacity(
        WRITE_BUF,
        std::fs::File::create(format!("{output_dir}combined_DNA_CDS.fasta"))
            .expect("create combined_DNA_CDS.fasta"),
    );
    let mut csv_handle = BufWriter::with_capacity(
        WRITE_BUF,
        std::fs::File::create(format!("{output_dir}gene_data.csv")).expect("create gene_data.csv"),
    );
    writeln!(
        csv_handle,
        "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description"
    )
    .expect("write header");

    // PORTING_PLAN.md §9 item 9: the Python chunks into batches of n_cpu and spawns a fresh
    // joblib pool per batch, then consumes each batch's results in order. One pool is
    // observationally identical because results are consumed in input order either way --
    // but the batching also bounds live results to n_cpu, which `parallel_map` would not.
    // `parallel_map_ordered_consume` keeps both properties: same sink order, same memory
    // bound, no per-batch barrier. See support::parallel.
    let jobs: Vec<(usize, String)> = gff_list.iter().cloned().enumerate().collect();
    crate::support::parallel::parallel_map_ordered_consume(
        n_cpu,
        jobs,
        |(gff_no, gff)| get_gene_sequences(&gff, gff_no, filter_seqs, &trans_table),
        |i, (dna, prot)| {
            output_files(
                &dna,
                &prot,
                &mut prot_handle,
                &mut dna_handle,
                &mut csv_handle,
                &gff_list[i],
            );
        },
    );
    true
}

/// Helper for the `entry.attributes["gene"][0]` / `["name"][0]` / `["product"]` reads in
/// `get_gene_sequences`, which fall back to `""` on `KeyError`.
///
/// Not a Python function — a `KeyError` idiom that has no direct Rust spelling.
pub fn attr_first_or_empty(entry: &GffEntry, key: &str) -> String {
    match entry.attribute(key) {
        Some(v) if !v.is_empty() => v[0].clone(),
        _ => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values from the reference Python.

    #[test]
    fn table_11_is_the_bacterial_table_unchanged() {
        // python: np.array_equal(get_trans_table(11)[0], bact_translation_table) -> True
        let tt = get_trans_table(11);
        assert_eq!(tt.table, BACT_TRANSLATION_TABLE);
        // python: sorted(get_trans_table(11)[1])
        let mut starts = tt.start_codons.clone();
        starts.sort();
        assert_eq!(
            starts,
            [
                "ATA", "ATC", "ATG", "ATT", "AUA", "AUC", "AUG", "AUU", "CTG", "CUG", "GTG", "GUG",
                "TTG", "UUG"
            ]
        );
    }

    #[test]
    fn table_4_patches_exactly_one_cell() {
        // python: the only cell differing from the bacterial table is
        //   (3, 2, 0)  b'*' -> b'W'      i.e. TGA, which is Trp in the mycoplasma code
        let tt = get_trans_table(4);
        let mut diffs = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    if tt.table[i][j][k] != BACT_TRANSLATION_TABLE[i][j][k] {
                        diffs.push((i, j, k, BACT_TRANSLATION_TABLE[i][j][k], tt.table[i][j][k]));
                    }
                }
            }
        }
        assert_eq!(diffs, [(3, 2, 0, b'*', b'W')]);
    }

    #[test]
    fn translate_matches_python_on_table_11() {
        let tt = get_trans_table(11);
        // Each expectation run through panaroo.prokka.translate.
        assert_eq!(translate("ATGAAATTTTAA", &tt), "MKF*");
        assert_eq!(translate("TTGAAATTT", &tt), "MKF"); // TTG is a start codon -> M
        assert_eq!(translate("GTGAAA", &tt), "MK"); // GTG likewise
        assert_eq!(translate("ATGNNNAAA", &tt), "MXK"); // N reduces to index 4 -> X
        assert_eq!(translate("atgaaa", &tt), "MK"); // lowercase accepted
        assert_eq!(translate("AAATTT", &tt), "KF"); // no start codon -> no M swap
        assert_eq!(translate("ATGAAATAGAAA", &tt), "MK*K"); // internal stop kept
    }

    #[test]
    fn translate_honours_the_selected_table() {
        // python: translate("ATGTGAAAA", tt11) == 'M*K';  with table 4 == 'MWK'
        assert_eq!(translate("ATGTGAAAA", &get_trans_table(11)), "M*K");
        assert_eq!(translate("ATGTGAAAA", &get_trans_table(4)), "MWK");
    }

    #[test]
    fn reduce_matches_the_numpy_lookup() {
        // python: [int(reduce_array[c]) for c in b"ACGTacgt"] -> [0,1,2,3,0,1,2,3]
        assert_eq!(
            b"ACGTacgt".iter().map(|&c| reduce(c)).collect::<Vec<_>>(),
            [0, 1, 2, 3, 0, 1, 2, 3]
        );
        // python: [int(reduce_array[c]) for c in b"N-*X"] -> [4,4,4,4]
        assert_eq!(
            b"N-*X".iter().map(|&c| reduce(c)).collect::<Vec<_>>(),
            [4, 4, 4, 4]
        );
    }

    #[test]
    fn clean_gff_string_drops_sequence_region_lines() {
        // python: clean_gff_string("##gff-version 3\n##sequence-region c1 1 99\n<CDS line>\n##sequence-region c2 1 5\n")
        //   -> '##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=a'
        let inp = "##gff-version 3\n##sequence-region c1 1 99\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=a\n##sequence-region c2 1 5\n";
        assert_eq!(
            clean_gff_string(inp),
            "##gff-version 3\nc1\tp\tCDS\t1\t9\t.\t+\t0\tID=a"
        );
    }
}
