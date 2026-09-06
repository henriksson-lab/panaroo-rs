//! Translation of `panaroo/panaroo/__main__.py`.
//!
//! `__main__.py::main` becomes [`panaroo_main`]; Rust's real `fn main` is a thin wrapper,
//! since the name is taken.

use panaroo::Args;

/// `__main__.py::SmartFormatter`
///
/// An `argparse.HelpFormatter` subclass that rewraps help text prefixed with `R|` at 55
/// columns. Only affects `--help` output, which is not part of the parity contract, but it
/// is kept so the module inventory matches.
pub struct SmartFormatter;

impl SmartFormatter {
    /// `SmartFormatter._split_lines`
    pub fn split_lines(&self, _text: &str, _width: usize) -> Vec<String> {
        panic!("noimpl: main::SmartFormatter::split_lines")
    }
}

/// `__main__.py::get_options`
///
/// Parses argv into [`Args`]. Options left unset stay `None` here and are filled in by
/// [`panaroo::set_default_args::set_default_args`] according to `--mode`.
pub fn get_options(argv: &[String]) -> Args {
    use clap::{Arg, ArgAction, Command};

    let m = Command::new("panaroo")
        .about("panaroo: an updated pipeline for pangenome investigation")
        .version(panaroo::VERSION) // argparse: --version, action="version"
        .arg(
            Arg::new("input_files")
                .short('i')
                .long("input")
                .required(true)
                .num_args(1..),
        )
        .arg(
            Arg::new("output_dir")
                .short('o')
                .long("out_dir")
                .required(true),
        )
        .arg(
            Arg::new("mode")
                .long("clean-mode")
                .required(true)
                .value_parser(["strict", "moderate", "sensitive"]),
        )
        .arg(
            Arg::new("filter_invalid")
                .long("remove-invalid-genes")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("id").short('c').long("threshold"))
        .arg(
            Arg::new("family_threshold")
                .short('f')
                .long("family_threshold"),
        )
        .arg(Arg::new("len_dif_percent").long("len_dif_percent"))
        .arg(Arg::new("family_len_dif_percent").long("family_len_dif_percent"))
        .arg(
            Arg::new("merge_paralogs")
                .long("merge_paralogs")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("search_radius").long("search_radius"))
        .arg(Arg::new("refind_prop_match").long("refind_prop_match"))
        .arg(
            Arg::new("refind_mode")
                .long("refind-mode")
                .value_parser(["default", "strict", "off"]),
        )
        .arg(Arg::new("min_trailing_support").long("min_trailing_support"))
        .arg(Arg::new("trailing_recursive").long("trailing_recursive"))
        .arg(Arg::new("edge_support_threshold").long("edge_support_threshold"))
        .arg(
            Arg::new("length_outlier_support_proportion").long("length_outlier_support_proportion"),
        )
        // `type=ast.literal_eval, choices=[True, False]`. Because `1 == True` and
        // `0 == False` in Python, the choices check also admits `1` and `0`.
        .arg(
            Arg::new("remove_by_consensus")
                .long("remove_by_consensus")
                .value_parser(["True", "False", "1", "0"]),
        )
        // NOTE: the flag is `--high_var_flag`; `cycle_threshold_min` is only the argparse
        // `dest`. Using the dest as the flag name breaks a working command in both
        // directions -- found by the CLI parity audit.
        .arg(Arg::new("cycle_threshold_min").long("high_var_flag"))
        .arg(Arg::new("min_edge_support_sv").long("min_edge_support_sv"))
        .arg(
            Arg::new("all_seq_in_graph")
                .long("all_seq_in_graph")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("clean_edges")
                .long("no_clean_edges")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("aln")
                .short('a')
                .long("alignment")
                .value_parser(["core", "pan"]),
        )
        .arg(Arg::new("alr").long("aligner").value_parser([
            "muscle",
            "muscle-super5",
            "famsa",
            "prank",
            "clustal",
            "mafft",
            "none",
        ]))
        .arg(Arg::new("codons").long("codons").action(ArgAction::SetTrue))
        // hyphen, not underscore -- upstream mixes the two and this one is hyphenated
        .arg(
            Arg::new("strict_codons")
                .long("strict-codons")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("core").long("core_threshold"))
        .arg(Arg::new("subset").long("core_subset"))
        .arg(Arg::new("hc_threshold").long("core_entropy_filter"))
        .arg(Arg::new("n_cpu").short('t').long("threads"))
        .arg(Arg::new("table").long("codon-table"))
        .arg(Arg::new("verbose").long("quiet").action(ArgAction::SetTrue))
        .get_matches_from(argv);

    let f = |k: &str| m.get_one::<String>(k).map(|s| s.to_string());
    let pf = |k: &str| f(k).map(|v| v.parse::<f64>().expect("float option"));
    let pi = |k: &str| f(k).map(|v| v.parse::<i64>().expect("int option"));

    Args {
        input_files: m
            .get_many::<String>("input_files")
            .unwrap()
            .cloned()
            .collect(),
        output_dir: f("output_dir").unwrap(),
        mode: f("mode").unwrap(),
        filter_invalid: m.get_flag("filter_invalid"),
        id: pf("id"),
        family_threshold: pf("family_threshold"),
        len_dif_percent: pf("len_dif_percent"),
        // `type=float` with NO `default=`, so this is None unless given -- the help text
        // claiming "default=0.0" is wrong. `set_default_args` never fills it either, so
        // cd-hit really is invoked with `-s None`. See PyNum::None.
        family_len_dif_percent: pf("family_len_dif_percent"),
        merge_paralogs: m.get_flag("merge_paralogs"),
        search_radius: pi("search_radius").unwrap_or(5000),
        refind_prop_match: pf("refind_prop_match").unwrap_or(0.2),
        refind_mode: f("refind_mode").unwrap_or_else(|| "default".into()),
        min_trailing_support: pi("min_trailing_support"),
        trailing_recursive: pi("trailing_recursive"),
        edge_support_threshold: pf("edge_support_threshold"),
        // UPSTREAM: the help text says default=0.01; the actual default is 0.1.
        // See ORIGINAL_CODE_BUG.md B5.
        length_outlier_support_proportion: pf("length_outlier_support_proportion").unwrap_or(0.1),
        remove_by_consensus: f("remove_by_consensus").map(|v| v == "True" || v == "1"),
        cycle_threshold_min: pi("cycle_threshold_min").unwrap_or(5),
        min_edge_support_sv: pi("min_edge_support_sv"),
        all_seq_in_graph: m.get_flag("all_seq_in_graph"),
        clean_edges: !m.get_flag("clean_edges"), // store_false, default True
        aln: f("aln"),
        alr: f("alr").unwrap_or_else(|| "mafft".into()),
        codons: m.get_flag("codons"),
        strict_codons: m.get_flag("strict_codons"),
        core: pf("core").unwrap_or(0.95),
        subset: pi("subset"),
        hc_threshold: pf("hc_threshold"),
        // kept signed: joblib reads a negative `n_jobs` as "all cores" (-1), "all but
        // one" (-2), and so on, and the raw value also reaches the cd-hit command string.
        n_cpu: pi("n_cpu").unwrap_or(1),
        table: pi("table").unwrap_or(11),
        verbose: !m.get_flag("verbose"), // store_false, default True
    }
}

/// `__main__.py::main`
///
/// The whole pipeline, in order:
///
/// ```text
///  1. get_options / set_default_args
///  2. check_cdhit_version, check_aligner_install
///  3. expand a file-of-filenames; convert GenBank via create_temp_gff3
///  4. process_prokka_input        -> gene_data.csv, combined_{DNA,protein}_CDS.fasta
///  5. run_cdhit                   -> combined_protein_cdhit_out.txt.clstr
///  6. generate_network            -> G, centroid_contexts, seqid_to_centroid
///  7. collapse_paralogs
///  8. write_gml                   -> pre_filt_graph.gml
///  9. collapse_families(correct_mistranslations=True)
/// 10. collapse_families(correct_mistranslations=False)
/// 11. trim_low_support_trailing_ends
/// 12. find_missing + collapse_families again   (unless --refind_mode off)
/// 13. clean_misassembly_edges, merge_paralogs
/// 14. generate_roary_gene_presence_absence, generate_summary_stats,
///     generate_pan_genome_reference, generate_common_struct_presence_absence
/// 15. write_gml                   -> final_graph.gml
/// 16. generate_{pan,core}_genome_alignment     (if --alignment)
/// ```
///
/// Steps 8 and 15 mutate node attributes in place before writing — `centroid`, `dna` and
/// `protein` become `";".join(...)` strings. See PORTING_PLAN.md §7 for how the Rust
/// handles that without losing the types.
pub fn panaroo_main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut args = get_options(&argv);

    // `set_default_args` runs at the END of `get_options`, i.e. before any of the checks
    // below -- not after them.
    panaroo::set_default_args::set_default_args(&mut args);

    // Check cd-hit is installed
    panaroo::cdhit::check_cdhit_version("cd-hit");
    // Check aligner is installed
    if args.aln.is_some() {
        panaroo::generate_alignments::check_aligner_install(&args.alr);
        // Get number of isolates if input file is provided.
        //
        // With a single `-i` argument the Python treats it as a file-of-filenames and
        // counts its LINES, not the one argument. Using `input_files.len()` here would
        // always yield 1 and so silently skip `check_aligner_sanity`'s large-input
        // warnings. Note it counts every line, blank ones included, and does this before
        // the file is parsed below.
        let no_isolates = if args.input_files.len() == 1 {
            let text = std::fs::read_to_string(&args.input_files[0])
                .unwrap_or_else(|e| panic!("could not read {}: {e}", args.input_files[0]));
            count_python_lines(&text)
        } else {
            args.input_files.len()
        };
        panaroo::generate_alignments::check_aligner_sanity(
            &args.alr,
            args.codons || args.strict_codons,
            no_isolates,
        );
    }

    // create directory if it isn't present already
    //
    // `os.mkdir`, not `os.makedirs`: a missing PARENT is an error rather than being
    // created silently, so `-o a/b/c` with no `a/b` fails fast instead of succeeding.
    if !std::path::Path::new(&args.output_dir).exists() {
        std::fs::create_dir(&args.output_dir).unwrap_or_else(|e| {
            panic!(
                "FileNotFoundError: could not create {}: {e}",
                args.output_dir
            )
        });
    }
    // make sure trailing forward slash is present
    if !args.output_dir.ends_with('/') {
        args.output_dir.push('/');
    }

    // Create temporary directory
    //
    // A randomly named `mkdtemp` directory INSIDE the output directory, not a fixed
    // `tmp_panaroo/`: concurrent runs sharing an output directory must not collide, and a
    // crashed run's scratch must not be reused.
    let temp_dir = format!(
        "{}/",
        panaroo::support::pytempfile::mkdtemp(std::path::Path::new(&args.output_dir))
            .unwrap_or_else(|e| panic!("could not create temp dir in {}: {e}", args.output_dir))
            .display()
    );
    // Child processes (cd-hit, mafft, prank, clustal) inherit this and put their own
    // scratch inside the run directory rather than in the system /tmp.
    std::env::set_var("TMPDIR", &temp_dir);

    // check if input is a file containing filenames
    if args.input_files.len() == 1 {
        let text = std::fs::read_to_string(&args.input_files[0])
            .unwrap_or_else(|e| panic!("could not read {}: {e}", args.input_files[0]));
        let mut files = Vec::new();
        for line in text.lines() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            match parts.len() {
                1 => {
                    let ext = std::path::Path::new(parts[0])
                        .extension()
                        .map(|e| format!(".{}", e.to_string_lossy()))
                        .unwrap_or_default();
                    if [".gbk", ".gb", ".gbff"].contains(&ext.as_str()) {
                        files.push(panaroo::prokka::create_temp_gff3(parts[0], None, &temp_dir));
                    } else if [".gff", ".gff3"].contains(&ext.as_str()) {
                        files.push(parts[0].to_string());
                    } else {
                        panic!("RuntimeError: Invalid file extension! ({ext})");
                    }
                }
                2 => files.push(panaroo::prokka::create_temp_gff3(
                    parts[0],
                    Some(parts[1]),
                    &temp_dir,
                )),
                // Also catches the EMPTY line: `"".strip().split()` is `[]`, whose length
                // is neither 1 nor 2, so a blank line in the file-of-filenames is a hard
                // error upstream rather than something to skip over.
                _ => {
                    // `print("Problem reading input line: ", line)` where `line` has
                    // already been rebound to the split list -- so what is printed is the
                    // list repr, and `print`'s separator makes the two spaces.
                    println!(
                        "Problem reading input line:  {}",
                        panaroo::support::pyfmt::py_repr_str_list(&parts)
                    );
                    panic!("RuntimeError: Error reading files!");
                }
            }
        }
        args.input_files = files;
    }

    if args.verbose {
        println!("pre-processing gff3 files...");
    }

    // convert input GFF3 files into summary files
    panaroo::prokka::process_prokka_input(
        &args.input_files,
        &args.output_dir,
        args.filter_invalid,
        !args.verbose,
        args.n_cpu,
        args.table,
    );

    // Cluster protein sequences using cdhit
    let cd_hit_out = format!("{}combined_protein_cdhit_out.txt", args.output_dir);
    panaroo::cdhit::run_cdhit(
        &format!("{}combined_protein_CDS.fasta", args.output_dir),
        &cd_hit_out,
        args.id.unwrap(),
        args.n_cpu,
        panaroo::support::pyfmt::PyNum::Float(args.len_dif_percent.unwrap()),
        panaroo::support::pyfmt::PyNum::Float(0.0),
        99999999,
        panaroo::support::pyfmt::PyNum::Float(0.0),
        99999999,
        true,
        false,
        None,
        None,
        !args.verbose,
    );

    if args.verbose {
        println!("generating initial network...");
    }

    // generate network from clusters and adjacency information
    let net = panaroo::generate_network::generate_network(
        &format!("{cd_hit_out}.clstr"),
        &format!("{}gene_data.csv", args.output_dir),
        &format!("{}combined_protein_CDS.fasta", args.output_dir),
        args.all_seq_in_graph,
    );
    let mut g = net.graph;
    let mut centroid_contexts = net.centroid_context;
    let seqid_to_centroid = net.seqid_to_centroid;

    // merge paralogs
    if args.verbose {
        println!("Processing paralogs...");
    }
    panaroo::clean_network::collapse_paralogs(&mut g, &mut centroid_contexts, 5, !args.verbose);

    // write out pre-filter graph in GML format.
    //
    // main assigns size/genomeIDs/geneIDs/degrees to every node here. Those Python dict
    // insertions fix the attribute order in BOTH GML files, so the flag records which nodes
    // got them before `name` was assigned later -- see NodeAttrs::gml_late_attrs_before_name.
    for n in g.nodes() {
        g.node_mut(n).gml_late_attrs_before_name = true;
    }
    write_graph_gml(&g, &format!("{}pre_filt_graph.gml", args.output_dir), true);

    if args.verbose {
        println!("collapse mistranslations...");
    }

    // clean up translation errors
    panaroo::clean_network::collapse_families(
        &mut g,
        &seqid_to_centroid,
        &temp_dir,
        0.7,
        0.98,
        args.family_len_dif_percent,
        true,
        args.length_outlier_support_proportion,
        args.n_cpu,
        !args.verbose,
        None,
        None,
        &[1, 2, 3],
        None,
    );

    if args.verbose {
        println!("collapse gene families...");
    }

    // collapse gene families
    let (dist, c2i) = panaroo::clean_network::collapse_families(
        &mut g,
        &seqid_to_centroid,
        &temp_dir,
        args.family_threshold.unwrap(),
        0.99,
        args.family_len_dif_percent,
        false,
        args.length_outlier_support_proportion,
        args.n_cpu,
        !args.verbose,
        None,
        None,
        &[1, 2, 3],
        None,
    );

    if args.verbose {
        println!("trimming contig ends...");
    }

    // re-trim low support trailing ends
    panaroo::clean_network::trim_low_support_trailing_ends(
        &mut g,
        args.min_trailing_support.unwrap() as usize,
        args.trailing_recursive.unwrap() as usize,
    );

    if g.number_of_nodes() < 2 {
        panic!(
            "RuntimeError: Nearly all clusters have been trimmed! Try reducing your \
             clustering sequence identity thresholds and/or running Panaroo in sensitive mode."
        );
    }

    if args.refind_mode != "off" {
        let only_valid_genes = args.refind_mode == "strict";
        if args.verbose {
            println!("refinding genes...");
        }
        panaroo::find_missing::find_missing(
            &mut g,
            &args.input_files,
            &format!("{}combined_DNA_CDS.fasta", args.output_dir),
            &format!("{}combined_protein_CDS.fasta", args.output_dir),
            &format!("{}gene_data.csv", args.output_dir),
            args.family_threshold.unwrap().max(0.8),
            args.search_radius,
            args.refind_prop_match,
            args.id.unwrap(),
            args.n_cpu,
            args.remove_by_consensus.unwrap(),
            only_valid_genes,
            args.verbose,
        );

        // merge again in case refinding has resolved issues
        if args.verbose {
            println!("collapse gene families with refound genes...");
        }
        panaroo::clean_network::collapse_families(
            &mut g,
            &seqid_to_centroid,
            &temp_dir,
            args.family_threshold.unwrap(),
            0.99,
            args.family_len_dif_percent,
            false,
            args.length_outlier_support_proportion,
            args.n_cpu,
            !args.verbose,
            Some(dist),
            Some(c2i),
            &[1, 2, 3],
            None,
        );
    }

    if args.clean_edges {
        panaroo::clean_network::clean_misassembly_edges(
            &mut g,
            args.edge_support_threshold.unwrap(),
        );
    }

    // if requested merge paralogs
    if args.merge_paralogs {
        panaroo::clean_network::merge_paralogs(&mut g);
    }

    let isolate_names: Vec<String> = args
        .input_files
        .iter()
        .map(|x| {
            let base = x.rsplit('/').next().unwrap_or(x);
            match base.rfind('.') {
                Some(i) if i > 0 => base[..i].to_string(),
                _ => base.to_string(),
            }
        })
        .collect();
    g.graph.isolate_names = isolate_names.clone();
    let mut mems_to_isolates: panaroo::support::pydict::PyDict<usize, String> =
        panaroo::support::pydict::PyDict::new();
    for (i, iso) in isolate_names.iter().enumerate() {
        mems_to_isolates.insert(i, iso.clone());
    }

    if args.verbose {
        println!("writing output...");
    }

    // get original annotation IDs, lengths and whether an internal stop codon is present
    let mut orig_ids: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut ids_len_stop: std::collections::HashMap<String, (usize, bool, bool)> =
        std::collections::HashMap::new();
    {
        let text = std::fs::read_to_string(format!("{}gene_data.csv", args.output_dir)).unwrap();
        for line in text.lines().skip(1) {
            let f: Vec<&str> = line.split(',').collect();
            if f.len() < 6 {
                continue;
            }
            orig_ids.insert(f[2].to_string(), f[3].to_string());
            let prot = f[4];
            // `"*" in line[4][1:-3]`
            let inner = if prot.len() > 4 {
                &prot[1..prot.len() - 3]
            } else {
                ""
            };
            ids_len_stop.insert(
                f[2].to_string(),
                (
                    prot.len(),
                    inner.contains('*'),
                    panaroo::isvalid::is_valid_gene(f[5], f[4]),
                ),
            );
        }
    }

    panaroo::generate_output::generate_roary_gene_presence_absence(
        &mut g,
        &mems_to_isolates,
        &orig_ids,
        &ids_len_stop,
        &args.output_dir,
    );
    panaroo::generate_output::generate_summary_stats(&args.output_dir);
    panaroo::generate_output::generate_pan_genome_reference(
        &g,
        &args.output_dir,
        &ids_len_stop,
        false,
    );
    panaroo::generate_output::generate_common_struct_presence_absence(
        &g,
        &args.output_dir,
        &mems_to_isolates,
        args.min_edge_support_sv.unwrap() as usize,
    );

    // add helpful attributes and write out graph in GML format
    write_graph_gml(&g, &format!("{}final_graph.gml", args.output_dir), false);

    // Write out core/pan-genome alignments
    match args.aln.as_deref() {
        Some("pan") => {
            if args.verbose {
                println!("generating pan genome MSAs...");
            }
            panaroo::generate_alignments::check_resume_manifest_collision(&args.output_dir, false);
            panaroo::generate_alignments::write_resume_manifest(
                &args.output_dir,
                "pan",
                &args.alr,
                args.codons,
                args.strict_codons,
                args.core,
                None,
                false,
            );
            panaroo::generate_output::generate_pan_genome_alignment(
                &g,
                &temp_dir,
                &args.output_dir,
                args.n_cpu,
                &args.alr,
                args.codons,
                args.strict_codons,
                &isolate_names,
                false,
            );
            if args.alr != "none" {
                let core_nodes = panaroo::generate_output::get_core_gene_nodes(
                    &g,
                    args.core,
                    args.input_files.len(),
                    None,
                );
                let core_names: Vec<String> = core_nodes
                    .iter()
                    .map(|&x| g.node(x).name.clone().expect("node name set"))
                    .collect();
                panaroo::generate_output::concatenate_core_genome_alignments(
                    &core_names,
                    &args.output_dir,
                    args.hc_threshold,
                );
            }
        }
        Some("core") => {
            if args.verbose {
                println!("generating core genome MSAs...");
            }
            panaroo::generate_alignments::check_resume_manifest_collision(&args.output_dir, false);
            panaroo::generate_alignments::write_resume_manifest(
                &args.output_dir,
                "core",
                &args.alr,
                args.codons,
                args.strict_codons,
                args.core,
                args.subset,
                false,
            );
            panaroo::generate_output::generate_core_genome_alignment(
                &g,
                &temp_dir,
                &args.output_dir,
                args.n_cpu,
                &args.alr,
                &isolate_names,
                args.core,
                args.codons,
                args.strict_codons,
                args.input_files.len(),
                args.hc_threshold,
                args.subset.map(|x| {
                    usize::try_from(x).expect("RuntimeError: --core_subset must be non-negative")
                }),
                false,
            );
        }
        _ => {}
    }

    // remove temporary directory
    let _ = std::fs::remove_dir_all(&temp_dir);
}

/// Assemble the node/edge attribute lists for a GML write and hand them to
/// [`panaroo::support::graph::write_gml`].
///
/// Not a Python function: `__main__.py` mutates the node attributes in place before each
/// `nx.write_gml` call — turning `centroid`/`dna`/`protein` into `";".join(...)` strings for
/// the final write, and converting `members`/`seqIDs` to lists. Doing that to a typed struct
/// would mean making every field dynamically typed, so the stringification happens here at
/// write time instead. **This is the one sanctioned deviation from "same logic"**
/// (PORTING_PLAN.md §7); the attribute *order* below is the Python's assignment order and
/// must not be reordered.
fn write_graph_gml(g: &panaroo::support::graph::Graph, path: &str, pre_filt: bool) {
    use panaroo::isvalid::{custom_stringizer, GmlValue};
    use panaroo::support::graph::{write_gml, GmlAttr};

    let mut nodes = Vec::new();
    for n in g.nodes() {
        let nd = g.node(n);
        let mut kvs: Vec<(String, GmlAttr)> = Vec::new();
        let push =
            |kvs: &mut Vec<(String, GmlAttr)>, k: &str, v: GmlAttr| kvs.push((k.to_string(), v));

        push(&mut kvs, "size", GmlAttr::Int(nd.members.len() as i64));
        if pre_filt {
            push(
                &mut kvs,
                "centroid",
                GmlAttr::List(
                    nd.centroid
                        .iter()
                        .map(|c| GmlAttr::Str(c.clone()))
                        .collect(),
                ),
            );
        } else {
            push(&mut kvs, "centroid", GmlAttr::Str(nd.centroid.join(";")));
        }
        push(&mut kvs, "maxLenId", GmlAttr::Int(nd.max_len_id as i64));
        if pre_filt {
            push(
                &mut kvs,
                "members",
                GmlAttr::Opaque(GmlValue::IntBitSet(
                    nd.members.iter().map(|m| m as i64).collect(),
                )),
            );
            push(
                &mut kvs,
                "seqIDs",
                GmlAttr::Opaque(GmlValue::Set(
                    nd.seq_ids
                        .iter()
                        .map(|s| GmlValue::Str(s.clone()))
                        .collect(),
                )),
            );
        } else {
            push(
                &mut kvs,
                "members",
                GmlAttr::List(nd.members.iter().map(|m| GmlAttr::Int(m as i64)).collect()),
            );
            push(
                &mut kvs,
                "seqIDs",
                GmlAttr::List(nd.seq_ids.iter().map(|s| GmlAttr::Str(s.clone())).collect()),
            );
        }
        push(&mut kvs, "hasEnd", GmlAttr::Bool(nd.has_end));
        if pre_filt {
            push(
                &mut kvs,
                "protein",
                GmlAttr::List(nd.protein.iter().map(|c| GmlAttr::Str(c.clone())).collect()),
            );
            push(
                &mut kvs,
                "dna",
                GmlAttr::List(nd.dna.iter().map(|c| GmlAttr::Str(c.clone())).collect()),
            );
        } else {
            push(&mut kvs, "protein", GmlAttr::Str(nd.protein.join(";")));
            push(&mut kvs, "dna", GmlAttr::Str(nd.dna.join(";")));
        }
        push(&mut kvs, "annotation", GmlAttr::Str(nd.annotation.clone()));
        push(
            &mut kvs,
            "description",
            GmlAttr::Str(nd.description.clone()),
        );
        push(
            &mut kvs,
            "lengths",
            GmlAttr::List(nd.lengths.iter().map(|&l| GmlAttr::Int(l as i64)).collect()),
        );
        push(
            &mut kvs,
            "longCentroidID",
            GmlAttr::List(vec![
                GmlAttr::Int(nd.long_centroid_id.0 as i64),
                GmlAttr::Str(nd.long_centroid_id.1.clone()),
            ]),
        );
        push(&mut kvs, "paralog", GmlAttr::Bool(nd.paralog));
        push(&mut kvs, "mergedDNA", GmlAttr::Bool(nd.merged_dna));
        // Attribute order: `name` goes last for a node that already had
        // genomeIDs/geneIDs/degrees when generate_roary assigned it, and first for a node
        // created by a later merge. See NodeAttrs::gml_late_attrs_before_name.
        let late = |kvs: &mut Vec<(String, GmlAttr)>| {
            kvs.push((
                "genomeIDs".into(),
                GmlAttr::Str(
                    nd.members
                        .iter()
                        .map(|m| m.to_string())
                        .collect::<Vec<_>>()
                        .join(";"),
                ),
            ));
            kvs.push((
                "geneIDs".into(),
                GmlAttr::Str(nd.seq_ids.iter().cloned().collect::<Vec<_>>().join(";")),
            ));
            kvs.push(("degrees".into(), GmlAttr::Int(g.degree(n) as i64)));
        };
        let name_attr = if pre_filt { None } else { nd.name.clone() };

        if nd.gml_late_attrs_before_name {
            late(&mut kvs);
            if let Some(name) = name_attr {
                push(&mut kvs, "name", GmlAttr::Str(name));
            }
        } else {
            if let Some(name) = name_attr {
                push(&mut kvs, "name", GmlAttr::Str(name));
            }
            late(&mut kvs);
        }

        nodes.push((n.to_string(), kvs));
    }

    let index: std::collections::HashMap<usize, usize> = g
        .nodes()
        .into_iter()
        .enumerate()
        .map(|(i, n)| (n, i))
        .collect();
    let mut edges = Vec::new();
    for (u, v) in g.edges() {
        let e = g.edge(u, v);
        let mut kvs: Vec<(String, GmlAttr)> = vec![("size".into(), GmlAttr::Int(e.size as i64))];
        if pre_filt {
            kvs.push((
                "members".into(),
                GmlAttr::Opaque(GmlValue::IntBitSet(
                    e.members.iter().map(|m| m as i64).collect(),
                )),
            ));
        } else {
            kvs.push((
                "members".into(),
                GmlAttr::List(e.members.iter().map(|m| GmlAttr::Int(m as i64)).collect()),
            ));
        }
        kvs.push((
            "genomeIDs".into(),
            GmlAttr::Str(
                e.members
                    .iter()
                    .map(|m| m.to_string())
                    .collect::<Vec<_>>()
                    .join(";"),
            ),
        ));
        edges.push((index[&u], index[&v], kvs));
    }

    let graph_attrs: Vec<(String, GmlAttr)> = if pre_filt {
        Vec::new()
    } else {
        vec![(
            "isolateNames".into(),
            GmlAttr::List(
                g.graph
                    .isolate_names
                    .iter()
                    .map(|s| GmlAttr::Str(s.clone()))
                    .collect(),
            ),
        )]
    };

    let sz: Option<&dyn Fn(&GmlValue) -> String> = if pre_filt {
        Some(&custom_stringizer)
    } else {
        None
    };
    write_gml(path, &nodes, &edges, &graph_attrs, sz);
}

/// Number of lines Python's `for line in file` would yield.
///
/// Differs from `str::lines()` on a trailing newline: Python yields no extra empty line
/// for a file ending in `"\n"`, and `lines()` agrees -- but Python DOES yield a line for a
/// file that is just `"\n"`, and so does `lines()`. The one real difference is the empty
/// file, which both give 0. Kept as a named function so the intent is checkable.
fn count_python_lines(text: &str) -> usize {
    if text.is_empty() {
        return 0;
    }
    text.lines().count()
}

fn main() {
    panaroo_main()
}
