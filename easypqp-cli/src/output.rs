use std::fs::{OpenOptions, metadata};
use std::io::{Write, BufWriter};
use std::path::Path;

use easypqp_core::{InsilicoPQPSettings, PeptideProperties};
use redeem_properties::utils::peptdeep_utils::remove_mass_shift;
use sage_core::peptide::Peptide;
use anyhow::Result;
use csv::ReaderBuilder;
use std::fs::File;
use std::io::BufReader;
use report_builder::{Report, ReportSection};
use maud::{html, Markup};
use std::collections::{HashMap, HashSet};
use plotly::{Plot, Scatter, Bar, Layout};
use redeem_properties::utils::peptdeep_utils::ion_mobility_to_ccs_bruker;


pub fn write_assays_to_tsv<P: AsRef<Path>>(
    assays: &[PeptideProperties],
    peptides: &[Peptide],
    path: P,
    insilico_settings: &InsilicoPQPSettings,
) -> std::io::Result<()> {
    let path = path.as_ref();
    let write_header = !path.exists() || metadata(path)?.len() == 0;

    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    let mut writer = BufWriter::new(file);

    if write_header {
        writeln!(
            writer,
            "PrecursorMz\t\
            ProductMz\t\
            PrecursorCharge\t\
            ProductCharge\t\
            LibraryIntensity\t\
            NormalizedRetentionTime\t\
            PeptideSequence\t\
            ModifiedPeptideSequence\t\
            PeptideGroupLabel\t\
            LabelType\t\
            CompoundName\t\
            SumFormula\t\
            SMILES\t\
            Adducts\t\
            ProteinId\t\
            UniprotId\t\
            GeneName\t\
            FragmentType\t\
            FragmentSeriesNumber\t\
            Annotation\t\
            CollisionEnergy\t\
            PrecursorIonMobility\t\
            TransitionGroupId\t\
            TransitionId\t\
            Decoy\t\
            DetectingTransition\t\
            IdentifyingTransition\t\
            QuantifyingTransition\t\
            Peptidoforms"
        )?;
    }

    // Decide whether to scale normalized RT (0-1) up to 0-100 for the TSV.
    // Some models/consumers expect RT in a 0-100 normalized scale (AlphaPeptiDeep),
    // while others use 0-1. Apply a harmless output-only heuristic: if the
    // maximum predicted RT across assays is <= 1.0, scale by 100 for output.
    let max_rt = assays
        .iter()
        .map(|a| a.retention_time)
        .fold(0.0_f32, |m, v| if v > m { v } else { m });
    // let rt_scale = if max_rt > 0.0 && max_rt <= 1.0 { 100.0_f32 } else { 1.0_f32 };
    let rt_scale = 100.0_f32;

    for assay in assays {
        let peptide_idx = assay.peptide_index as usize;
        let decoy = peptides[peptide_idx].decoy;
        let modified_peptide = peptides[peptide_idx].to_string();
        let naked_peptide = remove_mass_shift(&modified_peptide);
        let protein = &peptides[peptide_idx].proteins;

        // Parse protein entries which may be in the form `db|ACCESSION|GENE_ID` (e.g. sp|P26196|DDX6_HUMAN)
        // We want three output columns: ProteinId (accession(s)), UniprotId (duplicate of accession(s)),
        // and GeneName (gene id(s)). If an entry doesn't follow the pipe-separated format, fall back
        // to using the original string for ProteinId/UniProtId and leave GeneName empty for that entry.
        let protein_accessions: Vec<String> = protein
            .iter()
            .map(|p| {
                if p.contains('|') {
                    let parts: Vec<&str> = p.split('|').collect();
                    if parts.len() >= 2 {
                        parts[1].to_string()
                    } else {
                        p.to_string()
                    }
                } else {
                    p.to_string()
                }
            })
            .collect();

        let protein_genes: Vec<String> = protein
            .iter()
            .map(|p| {
                if p.contains('|') {
                    let parts: Vec<&str> = p.split('|').collect();
                    if parts.len() >= 3 {
                        parts[2].to_string()
                    } else {
                        "".to_string()
                    }
                } else {
                    "".to_string()
                }
            })
            .filter(|s| !s.is_empty())
            .collect();

        let protein_id_field = if protein_accessions.is_empty() {
            "".to_string()
        } else {
            protein_accessions.join(";")
        };

        // UniProtId will usually be the same accession(s); duplicate the accession field.
        let uniprot_id_field = protein_id_field.clone();

        let gene_name_field = if protein_genes.is_empty() {
            "".to_string()
        } else {
            protein_genes.join(";")
        };

        let non_zero_indices: Vec<usize> = assay.product.intensity
            .iter()
            .enumerate()
            .filter(|(_, &intensity)| intensity > 0.0)
            .map(|(i, _)| i)
            .collect();

        if non_zero_indices.len() < insilico_settings.min_transitions as usize {
            continue;
        }

        let product_indices = if insilico_settings.max_transitions > 0 {
            let mut sorted = non_zero_indices.clone();
            sorted.sort_by(|&a, &b| {
                assay.product.intensity[b]
                    .partial_cmp(&assay.product.intensity[a])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            sorted.into_iter()
                .take(insilico_settings.max_transitions as usize)
                .collect()
        } else {
            non_zero_indices
        };

    for &i in &product_indices {
            let fragment_type = &assay.product.ion_type[i];
            let series_number = assay.product.ion_ordinal[i];
            let product_charge = assay.product.charge[i];
            let annotation = format!("{}{}^{}", fragment_type, series_number, product_charge);

            writeln!(
                writer,
                "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                assay.precursor.precursor_mz,
                assay.product.product_mz[i],
                assay.precursor.charge,
                product_charge,
                assay.product.intensity[i] * 10_000.0,
                assay.retention_time * rt_scale,
                naked_peptide,
                modified_peptide,
                "", "", "", "", "", "", // Placeholder columns
                protein_id_field,
                uniprot_id_field,
                gene_name_field,
                fragment_type,
                series_number,
                annotation,
                "", // CollisionEnergy
                assay.precursor.ion_mobility,
                assay.peptide_index,
                assay.peptide_index * 1000 + i as u32,
                decoy as i8,
                1, 0, 1,
                "" // Peptidoforms
            )?;
        }
    }

    Ok(())
}

/// Generate a simple HTML report from the TSV library using `report-builder`.
pub fn generate_html_report<P: AsRef<Path>>(tsv_path: P, output_html: &str) -> Result<()> {
    let path = tsv_path.as_ref();

    let file = File::open(path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(BufReader::new(file));

    let headers = rdr
        .headers()?
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();

    // Read all rows for correct statistics/plots, but we'll keep a small sample for the sample table
    let mut all_rows: Vec<Vec<String>> = Vec::new();
    for result in rdr.records() {
        let record = result?;
        all_rows.push(record.iter().map(|c| c.to_string()).collect());
    }
    // small sample for the data table
    let sample_rows: Vec<Vec<String>> = all_rows.iter().take(20).cloned().collect();

    // Map headers to indices for easy lookup
    let mut idx: HashMap<String, usize> = HashMap::new();
    for (i, h) in headers.iter().enumerate() {
        idx.insert(h.clone(), i);
    }

    // Helper to get column value by header name
    let get_str = |row: &Vec<String>, name: &str| -> Option<String> {
        idx.get(name).and_then(|&i| row.get(i)).cloned()
    };

    // Build TargetDecoy and typed columns
    let mut target_flags: Vec<String> = Vec::new();
    let mut lib_intensity: Vec<Option<f64>> = Vec::new();
    let mut rt_vals: Vec<Option<f64>> = Vec::new();
    let mut im_vals: Vec<Option<f64>> = Vec::new();
    let mut precursor_mz: Vec<Option<f64>> = Vec::new();
    let mut precursor_charge: Vec<Option<i32>> = Vec::new();
    let mut product_charge: Vec<Option<i32>> = Vec::new();
    let mut fragment_type: Vec<Option<String>> = Vec::new();
    let mut peptide_seq: Vec<Option<String>> = Vec::new();
    let mut protein_id: Vec<Option<String>> = Vec::new();
    let mut gene_name: Vec<Option<String>> = Vec::new();

    for row in &all_rows {
        // Decoy -> Target/Decoy
    let decoy = get_str(row, "Decoy").unwrap_or_else(|| "0".to_string());
    let td = if decoy == "1" { "Decoy" } else { "Target" };
        target_flags.push(td.to_string());

        lib_intensity.push(get_str(row, "LibraryIntensity").and_then(|s| s.parse::<f64>().ok()));
        rt_vals.push(get_str(row, "NormalizedRetentionTime").and_then(|s| s.parse::<f64>().ok()));
        im_vals.push(get_str(row, "PrecursorIonMobility").and_then(|s| s.parse::<f64>().ok()));
        precursor_mz.push(get_str(row, "PrecursorMz").and_then(|s| s.parse::<f64>().ok()));
        precursor_charge.push(get_str(row, "PrecursorCharge").and_then(|s| s.parse::<i32>().ok()));
        product_charge.push(get_str(row, "ProductCharge").and_then(|s| s.parse::<i32>().ok()));
        fragment_type.push(get_str(row, "FragmentType").map(|s| s.to_string()));
        peptide_seq.push(get_str(row, "PeptideSequence").map(|s| s.to_string()));
        protein_id.push(get_str(row, "ProteinId").map(|s| s.to_string()));
        gene_name.push(get_str(row, "GeneName").map(|s| s.to_string()));
    }

    // Basic statistics
    let mut unique_target_peptides: HashSet<String> = HashSet::new();
    let mut unique_decoy_peptides: HashSet<String> = HashSet::new();
    let mut unique_target_proteins: HashSet<String> = HashSet::new();
    let mut unique_decoy_proteins: HashSet<String> = HashSet::new();
    let mut target_transitions = 0usize;
    let mut decoy_transitions = 0usize;

    for (i, td) in target_flags.iter().enumerate() {
        let pep = peptide_seq[i].clone().unwrap_or_default();
        let prot = protein_id[i].clone().unwrap_or_default();
        if td == "Target" {
            target_transitions += 1;
            if !pep.is_empty() { unique_target_peptides.insert(pep); }
            if !prot.is_empty() { unique_target_proteins.insert(prot); }
        } else {
            decoy_transitions += 1;
            if !pep.is_empty() { unique_decoy_peptides.insert(pep); }
            if !prot.is_empty() { unique_decoy_proteins.insert(prot); }
        }
    }

    // Intensity stats by Target/Decoy
    fn mean_sd(v: &Vec<f64>) -> (f64, f64) {
        let n = v.len() as f64;
        if n == 0.0 { return (f64::NAN, f64::NAN); }
        let mean = v.iter().sum::<f64>() / n;
        let var = v.iter().map(|x| (x-mean).powi(2)).sum::<f64>() / n;
        (mean, var.sqrt())
    }

    let mut intensity_target: Vec<f64> = Vec::new();
    let mut intensity_decoy: Vec<f64> = Vec::new();
    for (i, opt) in lib_intensity.iter().enumerate() {
        if let Some(val) = opt {
            if target_flags[i] == "Target" { intensity_target.push(*val); } else { intensity_decoy.push(*val); }
        }
    }

    let (mean_t, sd_t) = mean_sd(&intensity_target);
    let (mean_d, sd_d) = mean_sd(&intensity_decoy);

    // Build Plotly plots
    // Normalize RT values to 0-100 if they look like absolute minutes or wider-range values.
    // This mirrors the Rscript / tuning behavior which maps raw RTs into a 0-100 normalized scale.
    let rt_numeric: Vec<f64> = rt_vals.iter().filter_map(|o| *o).collect();
    let mut normalized_rt_vals: Vec<Option<f64>> = Vec::with_capacity(rt_vals.len());
    if !rt_numeric.is_empty() {
        let rt_min = rt_numeric.iter().cloned().fold(std::f64::INFINITY, f64::min);
        let rt_max = rt_numeric.iter().cloned().fold(std::f64::NEG_INFINITY, f64::max);
        let rt_range = rt_max - rt_min;
        // If the range is significant (values > 1 or range > 1), normalize to 0-100.
        if rt_max > 1.0 || rt_range > 1.0 {
            for v in &rt_vals {
                if let Some(x) = v {
                    if rt_range > 0.0 {
                        normalized_rt_vals.push(Some(100.0 * (*x - rt_min) / rt_range));
                    } else {
                        normalized_rt_vals.push(Some(50.0));
                    }
                } else {
                    normalized_rt_vals.push(None);
                }
            }
        } else {
            // Already normalized between 0-1 or small; scale to 0-100 for plotting
            for v in &rt_vals {
                normalized_rt_vals.push(v.map(|x| x * 100.0));
            }
        }
    } else {
        normalized_rt_vals = rt_vals.clone();
    }
    // Library intensity histograms (log10) — mirror R script: bins = 50 and separate panels for Target/Decoy
    let x_t: Vec<f64> = intensity_target.iter().map(|v| v.log10()).collect();
    let x_d: Vec<f64> = intensity_decoy.iter().map(|v| v.log10()).collect();

    // Combined intensity histogram: target and decoy overlaid, colored
    let mut p_intensity = Plot::new();
    let mut hist_t = plotly::Histogram::new(x_t.clone()).name("Target");
    hist_t = hist_t.n_bins_x(50).opacity(0.6);
    let mut hist_d = plotly::Histogram::new(x_d.clone()).name("Decoy");
    hist_d = hist_d.n_bins_x(50).opacity(0.6);
    p_intensity.add_trace(hist_t);
    p_intensity.add_trace(hist_d);
    // set overlay mode by setting barmode overlay if available
    p_intensity.set_layout(Layout::new().title("Library Intensity Distribution (log10)").clone());

    // RT histograms (normalized) — bins = 50 and separate panels for Target/Decoy
    let rt_t: Vec<f64> = normalized_rt_vals.iter().enumerate().filter_map(|(i, o)| o.map(|v| if target_flags[i]=="Target" { v } else { std::f64::NAN })).filter(|v| v.is_finite()).collect();
    let rt_d: Vec<f64> = normalized_rt_vals.iter().enumerate().filter_map(|(i, o)| o.map(|v| if target_flags[i]=="Decoy" { v } else { std::f64::NAN })).filter(|v| v.is_finite()).collect();

    // Combined RT histogram: overlay Target/Decoy in same plot
    let mut p_rt = Plot::new();
    let mut rt_hist_t = plotly::Histogram::new(rt_t.clone()).name("Target");
    rt_hist_t = rt_hist_t.n_bins_x(50).opacity(0.6);
    let mut rt_hist_d = plotly::Histogram::new(rt_d.clone()).name("Decoy");
    rt_hist_d = rt_hist_d.n_bins_x(50).opacity(0.6);
    p_rt.add_trace(rt_hist_t);
    p_rt.add_trace(rt_hist_d);
    p_rt.set_layout(Layout::new().title("Normalized Retention Time Distribution"));

    // Ion mobility scatter — if values look like mobility (very small), convert back to CCS for plotting
    let mut p_im = Plot::new();
    let mut x_t_mz: Vec<f64> = Vec::new();
    let mut y_t_im: Vec<f64> = Vec::new();
    let mut x_d_mz: Vec<f64> = Vec::new();
    let mut y_d_im: Vec<f64> = Vec::new();

    // Determine if IM values are small (mobility) or large (CCS). If max < 0.01 assume mobility and convert to CCS
    let max_im = im_vals.iter().filter_map(|o| *o).fold(0.0_f64, |m, v| if v > m { v } else { m });
    let convert_to_ccs = max_im > 0.0 && max_im < 0.01;

    for i in 0..all_rows.len() {
        if let (Some(mz), Some(im)) = (precursor_mz[i], im_vals[i]) {
            if convert_to_ccs {
                // require precursor charge to compute CCS
                if let Some(ch) = precursor_charge[i] {
                    let ccs = ion_mobility_to_ccs_bruker(im as f64, ch as i32, mz as f64);
                    let ccs_f = ccs as f64;
                    if target_flags[i] == "Target" { x_t_mz.push(mz); y_t_im.push(ccs_f); } else { x_d_mz.push(mz); y_d_im.push(ccs_f); }
                }
            } else {
                if target_flags[i] == "Target" { x_t_mz.push(mz); y_t_im.push(im); } else { x_d_mz.push(mz); y_d_im.push(im); }
            }
        }
    }

    if convert_to_ccs {
        p_im.add_trace(Scatter::new(x_t_mz.clone(), y_t_im.clone()).name("Target").mode(plotly::common::Mode::Markers));
        p_im.add_trace(Scatter::new(x_d_mz.clone(), y_d_im.clone()).name("Decoy").mode(plotly::common::Mode::Markers));
        p_im.set_layout(Layout::new().title("Ion Mobility (converted to CCS) vs Precursor m/z"));
    } else {
        p_im.add_trace(Scatter::new(x_t_mz.clone(), y_t_im.clone()).name("Target").mode(plotly::common::Mode::Markers));
        p_im.add_trace(Scatter::new(x_d_mz.clone(), y_d_im.clone()).name("Decoy").mode(plotly::common::Mode::Markers));
        p_im.set_layout(Layout::new().title("Ion Mobility vs Precursor m/z"));
    }

    // Fragment type counts
    let mut frag_counts_target: HashMap<String, usize> = HashMap::new();
    let mut frag_counts_decoy: HashMap<String, usize> = HashMap::new();
    for i in 0..all_rows.len() {
        if let Some(ft) = &fragment_type[i] {
            if target_flags[i] == "Target" { *frag_counts_target.entry(ft.clone()).or_default() += 1; } else { *frag_counts_decoy.entry(ft.clone()).or_default() += 1; }
        }
    }
    // Build bar plot of fragment types (targets stacked with decoys)
    let mut frag_types: Vec<String> = frag_counts_target.keys().chain(frag_counts_decoy.keys()).map(|s| s.clone()).collect();
    frag_types.sort();
    frag_types.dedup();
    let x_frag = frag_types.clone();
    let y_t_frag: Vec<i32> = frag_types.iter().map(|k| *frag_counts_target.get(k).unwrap_or(&0) as i32).collect();
    let y_d_frag: Vec<i32> = frag_types.iter().map(|k| *frag_counts_decoy.get(k).unwrap_or(&0) as i32).collect();
    let mut p_frag = Plot::new();
    p_frag.add_trace(Bar::new(x_frag.clone(), y_t_frag.clone()).name("Target"));
    p_frag.add_trace(Bar::new(x_frag.clone(), y_d_frag.clone()).name("Decoy"));
    p_frag.set_layout(Layout::new().title("Fragment Type Distribution"));

    // Target/Decoy intensity stats table
    let stats_table = vec![
        ("Target", intensity_target.len(), mean_t, sd_t),
        ("Decoy", intensity_decoy.len(), mean_d, sd_d),
    ];

    // Top proteins by number of transitions
    let mut protein_map: HashMap<(String,String), (usize, HashSet<String>)> = HashMap::new();
    for i in 0..all_rows.len() {
        let pid = protein_id[i].clone().unwrap_or_default();
        let gname = gene_name[i].clone().unwrap_or_default();
        if pid.is_empty() { continue; }
        let key = (pid.clone(), gname.clone());
        let entry = protein_map.entry(key).or_insert((0usize, HashSet::new()));
        entry.0 += 1;
        if let Some(p) = &peptide_seq[i] { if !p.is_empty() { entry.1.insert(p.clone()); } }
    }
    let mut protein_stats: Vec<(String,String,usize,usize)> = protein_map.into_iter().map(|((pid,gname),(count,peps))| (pid,gname,count,peps.len())).collect();
    protein_stats.sort_by(|a,b| b.2.cmp(&a.2));

    // Build report sections
    let mut report = Report::new("easypqp", clap::crate_version!(), None, "In-Silico Library Report");

    // Basic Stats section
    let mut s_basic = ReportSection::new("Basic Statistics");
    let basic_markup: Markup = html! {
        p { (format!("Number of unique target peptides: {} and decoy peptides: {}", unique_target_peptides.len(), unique_decoy_peptides.len())) }
        p { (format!("Number of unique target proteins: {} and decoy proteins: {}", unique_target_proteins.len(), unique_decoy_proteins.len())) }
    p { (format!("Number of transitions: {}", all_rows.len())) }
        p { (format!("Number of target transitions: {}", target_transitions)) }
        p { (format!("Number of decoy transitions: {}", decoy_transitions)) }
    };
    s_basic.add_content(basic_markup);
    report.add_section(s_basic);

    // Plots section
    let mut s_plots = ReportSection::new("Plots");
    // Intensity: combined Target/Decoy overlay
    s_plots.add_plot(p_intensity);
    // RT: combined Target/Decoy overlay
    s_plots.add_plot(p_rt);
    // Ion mobility / CCS
    s_plots.add_plot(p_im);
    // Fragment type distribution
    s_plots.add_plot(p_frag);
    report.add_section(s_plots);

    // Stats table section
    let mut s_stats = ReportSection::new("Target/Decoy Statistics");
    let stats_markup: Markup = html! {
        table { thead { tr { th { "Group" } th { "N" } th { "Mean Intensity" } th { "SD Intensity" } } } tbody {
            @for (grp,n,mean,sd) in &stats_table {
                tr { td { (grp) } td { (n.to_string()) } td { (format!("{:.3}", mean)) } td { (format!("{:.3}", sd)) } }
            }
        } }
    };
    s_stats.add_content(stats_markup);
    report.add_section(s_stats);

    // Top proteins
    let mut s_prot = ReportSection::new("Top Proteins by Number of Transitions");
    let prot_markup: Markup = html! {
        table { thead { tr { th { "ProteinId" } th { "GeneName" } th { "Transitions" } th { "UniquePeptides" } } } tbody {
            @for row in protein_stats.iter().take(20) {
                tr { td { (&row.0) } td { (&row.1) } td { (row.2) } td { (row.3) } }
            }
        } }
    };
    s_prot.add_content(prot_markup);
    report.add_section(s_prot);

    // Sample data
    let mut s_sample = ReportSection::new("Sample of Library Data");
    // reuse previous table building approach for sample rows
    let sample_markup: Markup = html! {
        table class="display" id="dataTable" {
            thead { tr { @for h in &headers { th { (h) } } } }
            tbody { @for row in sample_rows.iter() { tr { @for cell in row { td { (cell) } } } } }
        }
    };
    s_sample.add_content(sample_markup);
    report.add_section(s_sample);

    report.save_to_file(output_html)?;

    Ok(())
}