//! Translation of `panaroo/panaroo/biocode_convert.py`.
//!
//! # Vendored third-party code
//!
//! This whole module is **not Panaroo's own**: it is a vendored copy of the GenBank-to-GFF3
//! converter from [biocode](https://github.com/jorvis/biocode) (MIT), and it still imports
//! `biocode.annotation`, `biocode.things` and `biocode.utils` at line 21 — so `biocode` is a
//! hard runtime dependency of `python -m panaroo` despite appearances. See `NOTICE.md`.
//!
//! Only exercised when the input list names `.gbk` / `.gb` / `.gbff` files.

/// `biocode_convert.py::convert_gbk_gff3`
pub fn convert_gbk_gff3(input_file: &str, output_file: &str, fasta: bool) {
    use crate::support::biocode_gff3::{wrapped_fasta, Row};
    use crate::support::pydict::PyDict;
    use std::collections::HashMap;
    use std::io::Write;

    let f =
        std::fs::File::create(output_file).unwrap_or_else(|e| panic!("create {output_file}: {e}"));
    let mut ofh = std::io::BufWriter::new(f);
    writeln!(ofh, "##gff-version 3").unwrap();

    // assemblies is a plain dict -- insertion-ordered, and the ##FASTA block follows that
    // order.
    let mut assemblies: PyDict<String, String> = PyDict::new();
    let mut current_assembly;
    let mut seqs_pending_writes = false;
    let mut features_skipped_count = 0usize;

    let mut rna_count_by_gene: HashMap<String, usize> = HashMap::new();
    let mut exon_count_by_RNA: HashMap<String, usize> = HashMap::new();

    // The Python buffers one gene's rows in a `things.Gene` and flushes it when the next
    // `gene` feature arrives (and once more at the end).
    let mut gene_rows: Vec<Row> = Vec::new();
    let mut have_gene = false;
    let mut current_gene_lt = String::new();
    let mut current_rna_id: Option<String> = None;
    let mut pending_polypeptide: Option<crate::support::biocode_gff3::Row> = None;

    for gb_record in crate::support::seqio::parse_genbank_file(input_file) {
        let mol_id = gb_record.name.clone();
        if !assemblies.contains_key(&mol_id) {
            assemblies.insert(mol_id.clone(), String::new());
        }
        if !gb_record.seq.is_empty() {
            seqs_pending_writes = true;
            assemblies.insert(mol_id.clone(), gb_record.seq.clone());
        }
        current_assembly = mol_id.clone();

        for feat in &gb_record.features {
            let fmin = feat.start;
            let fmax = feat.end;
            let strand = match feat.strand {
                1 => '+',
                -1 => '-',
                _ => panic!("Exception: ERROR: unstranded feature encountered: {feat:?}"),
            };

            match feat.feature_type.as_str() {
                "source" => continue,
                "gene" => {
                    // print the previous gene (if there is one)
                    if have_gene {
                        flush(&mut ofh, &mut gene_rows);
                    }
                    let lt = feat
                        .qual("locus_tag")
                        .expect("KeyError: locus_tag")
                        .to_string();
                    gene_rows.push(Row {
                        seqid: current_assembly.clone(),
                        feature_type: "gene".into(),
                        fmin,
                        fmax,
                        strand,
                        phase: None,
                        attributes: vec![
                            ("ID".into(), lt.clone()),
                            ("locus_tag".into(), lt.clone()),
                        ],
                    });
                    current_gene_lt = lt;
                    have_gene = true;
                    current_rna_id = None;
                }
                "mRNA" => {
                    let lt = feat
                        .qual("locus_tag")
                        .expect("KeyError: locus_tag")
                        .to_string();
                    let n = rna_count_by_gene.entry(lt.clone()).or_insert(0);
                    *n += 1;
                    let feat_id = format!("{lt}.mRNA.{n}");
                    gene_rows.push(Row {
                        seqid: current_assembly.clone(),
                        feature_type: "mRNA".into(),
                        fmin,
                        fmax,
                        strand,
                        phase: None,
                        attributes: vec![
                            ("ID".into(), feat_id.clone()),
                            ("Parent".into(), current_gene_lt.clone()),
                        ],
                    });
                    if exon_count_by_RNA.contains_key(&feat_id) {
                        panic!(
                            "Exception: ERROR: two different RNAs found with same ID: {feat_id}"
                        );
                    }
                    exon_count_by_RNA.insert(feat_id.clone(), 0);
                    current_rna_id = Some(feat_id);
                }
                "tRNA" | "rRNA" => {
                    let kind = &feat.feature_type;
                    let lt = feat
                        .qual("locus_tag")
                        .expect("KeyError: locus_tag")
                        .to_string();
                    let n = rna_count_by_gene.entry(lt.clone()).or_insert(0);
                    *n += 1;
                    let feat_id = format!("{lt}.{kind}.{n}");
                    let mut attrs = vec![
                        ("ID".into(), feat_id.clone()),
                        ("Parent".into(), current_gene_lt.clone()),
                    ];
                    // tRNA carries the product as `anticodon`; rRNA as `product_name`.
                    if let Some(p) = feat.qual("product") {
                        let key = if kind == "tRNA" {
                            "anticodon"
                        } else {
                            "product_name"
                        };
                        attrs.push((key.into(), p.to_string()));
                    }
                    gene_rows.push(Row {
                        seqid: current_assembly.clone(),
                        feature_type: kind.clone(),
                        fmin,
                        fmax,
                        strand,
                        phase: None,
                        attributes: attrs,
                    });
                    if exon_count_by_RNA.contains_key(&feat_id) {
                        panic!(
                            "Exception: ERROR: two different RNAs found with same ID: {feat_id}"
                        );
                    }
                    exon_count_by_RNA.insert(feat_id.clone(), 0);
                    current_rna_id = Some(feat_id);
                }
                "CDS" => {
                    let lt = feat
                        .qual("locus_tag")
                        .expect("KeyError: locus_tag")
                        .to_string();
                    // In a prokaryotic GBK the CDS arrives before any mRNA, so one is made
                    // here. Note `rna_count_by_gene[lt]` is *not* incremented, which is why
                    // these ids end in `.0`.
                    if current_rna_id.is_none() {
                        let n = *rna_count_by_gene.get(&lt).unwrap_or(&0);
                        let feat_id = format!("{lt}.mRNA.{n}");
                        gene_rows.push(Row {
                            seqid: current_assembly.clone(),
                            feature_type: "mRNA".into(),
                            fmin,
                            fmax,
                            strand,
                            phase: None,
                            attributes: vec![
                                ("ID".into(), feat_id.clone()),
                                ("Parent".into(), current_gene_lt.clone()),
                            ],
                        });
                        exon_count_by_RNA.entry(feat_id.clone()).or_insert(0);
                        current_rna_id = Some(feat_id.clone());

                        // the polypeptide carries the functional annotation
                        let mut attrs = vec![
                            ("ID".into(), format!("{lt}.polypeptide.{n}")),
                            ("Parent".into(), feat_id),
                        ];
                        if let Some(g) = feat.qual("gene") {
                            attrs.push(("gene_symbol".into(), g.to_string()));
                        }
                        if let Some(p) = feat.qual("product") {
                            attrs.push(("product_name".into(), p.to_string()));
                        }
                        // deferred: the polypeptide row is emitted after the CDS/exon rows
                        pending_polypeptide = Some(Row {
                            seqid: current_assembly.clone(),
                            feature_type: "polypeptide".into(),
                            fmin,
                            fmax,
                            strand,
                            phase: None,
                            attributes: attrs,
                        });
                    }

                    let rna_id = current_rna_id.clone().expect("current RNA set");
                    let c = exon_count_by_RNA.entry(rna_id.clone()).or_insert(0);
                    *c += 1;
                    let cds_id = format!("{rna_id}.CDS.{c}");
                    let mut current_cds_phase: i64 = 0;

                    for &(subfmin, subfmax) in &feat.parts {
                        gene_rows.push(Row {
                            seqid: current_assembly.clone(),
                            feature_type: "CDS".into(),
                            fmin: subfmin,
                            fmax: subfmax,
                            strand,
                            phase: Some(current_cds_phase as u8),
                            attributes: vec![
                                ("ID".into(), cds_id.clone()),
                                ("Parent".into(), rna_id.clone()),
                            ],
                        });

                        // general: 3 - ((length - previous phase) % 3)
                        current_cds_phase = 3 - (((subfmax - subfmin) - current_cds_phase) % 3);
                        if current_cds_phase == 3 {
                            current_cds_phase = 0;
                        }

                        let n = *exon_count_by_RNA.get(&rna_id).unwrap();
                        gene_rows.push(Row {
                            seqid: current_assembly.clone(),
                            feature_type: "exon".into(),
                            fmin: subfmin,
                            fmax: subfmax,
                            strand,
                            phase: None,
                            attributes: vec![
                                ("ID".into(), format!("{rna_id}.exon.{n}")),
                                ("Parent".into(), rna_id.clone()),
                            ],
                        });
                        *exon_count_by_RNA.get_mut(&rna_id).unwrap() += 1;
                    }

                    if let Some(p) = pending_polypeptide.take() {
                        gene_rows.push(p);
                    }
                }
                _ => {
                    println!("WARNING: The following feature was skipped:\n{feat:?}");
                    features_skipped_count += 1;
                }
            }
        }
    }

    // don't forget to do the last gene, if there were any
    if have_gene {
        flush(&mut ofh, &mut gene_rows);
    }

    if fasta && seqs_pending_writes {
        writeln!(ofh, "##FASTA").unwrap();
        for (assembly_id, residues) in assemblies.items() {
            writeln!(ofh, ">{assembly_id}").unwrap();
            writeln!(ofh, "{}", wrapped_fasta(residues)).unwrap();
        }
    }

    if features_skipped_count > 0 {
        println!("Warning: {features_skipped_count} unsupported feature types were skipped");
    }
}

fn flush<W: std::io::Write>(ofh: &mut W, rows: &mut Vec<crate::support::biocode_gff3::Row>) {
    for r in rows.iter() {
        writeln!(ofh, "{}", r.render()).unwrap();
    }
    rows.clear();
}
