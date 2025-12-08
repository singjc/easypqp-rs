use std::fs::{OpenOptions, metadata};
use std::io::{Write, BufWriter};
use std::path::Path;

use easypqp_core::{InsilicoPQPSettings, PeptideProperties};
use redeem_properties::utils::peptdeep_utils::remove_mass_shift;
use sage_core::peptide::Peptide;


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