use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use redeem_properties::utils::data_handling::PeptideData;
use redeem_properties::utils::peptdeep_utils::{
    get_modification_indices, get_modification_string, ion_mobility_to_ccs_bruker,
    load_modifications, remove_mass_shift,
};
use regex::{Captures, Regex};

fn canonicalize_mass_shift_annotations(sequence: &str) -> String {
    let re = Regex::new(r"\[([+-]?\d*\.?\d+)\]").unwrap();
    re.replace_all(sequence, |caps: &Captures| {
        let mass = caps.get(1).map(|m| m.as_str()).unwrap_or_default();
        if mass.starts_with('+') || mass.starts_with('-') {
            format!("[{}]", mass)
        } else {
            format!("[+{}]", mass)
        }
    })
    .into_owned()
}

fn normalize_sequence_for_mod_parsing(sequence: &str) -> String {
    let canonical_sequence = canonicalize_mass_shift_annotations(sequence.trim().trim_matches('.'));
    let n_term_separator = Regex::new(r"^((?:\[[+-]?\d*\.?\d+\]|\(UniMod:\d+\))+)-").unwrap();
    n_term_separator
        .replace(&canonical_sequence, "$1")
        .into_owned()
}

pub fn read_peptide_data_from_tsv<P: AsRef<Path>>(
    path: P,
    nce: i32,
    instrument: &str,
) -> std::io::Result<Vec<PeptideData>> {
    let modifications = load_modifications()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;

    let file = File::open(&path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read header and map column names
    let header = match lines.next() {
        Some(h) => h?,
        None => return Ok(Vec::new()),
    };

    let header_columns: Vec<&str> = header.split('\t').collect();
    let mut column_indices = HashMap::new();

    // Flexible column name matching
    for (idx, col) in header_columns.iter().enumerate() {
        match col.to_lowercase().as_str() {
            s if s.contains("modifiedpeptide") || s.contains("fullpeptidename") => {
                column_indices.insert("sequence", idx);
            }
            s if s.contains("precursormz") || s.contains("precursor_mz") => {
                column_indices.insert("precursor_mz", idx);
            }
            s if s.contains("precursorcharge")
                || s.contains("charge")
                || s.contains("precursor_charge") =>
            {
                column_indices.insert("precursor_charge", idx);
            }
            s if s.contains("fragmenttype") || s.contains("fragment_type") => {
                column_indices.insert("fragment_type", idx);
            }
            s if s.contains("fragmentseriesnumber") || s.contains("fragment_series_number") => {
                column_indices.insert("fragment_series_number", idx);
            }
            s if s.contains("productcharge") || s.contains("product_charge") => {
                column_indices.insert("product_charge", idx);
            }
            s if s.contains("libraryintensity") || s.contains("intensity") => {
                column_indices.insert("intensity", idx);
            }
            s if s.contains("normalizedretentiontime")
                || s.contains("rt")
                || s.contains("retention_time") =>
            {
                column_indices.insert("retention_time", idx);
            }
            s if s.contains("precursorionmobility")
                || s.contains("im")
                || s.contains("ion_mobility") =>
            {
                column_indices.insert("ion_mobility", idx);
            }
            // s if s.contains("collisionenergy") || s.contains("nce") => {
            //     column_indices.insert("nce", idx);
            // }
            // s if s.contains("instrument") => {
            //     column_indices.insert("instrument", idx);
            // }
            _ => (),
        }
    }

    // Verify required columns
    let required_columns = [
        "sequence",
        "precursor_charge",
        "intensity",
        "retention_time",
    ];
    for col in required_columns {
        if !column_indices.contains_key(col) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Missing required column: {}", col),
            ));
        }
    }

    // First pass - collect all retention times for normalization
    let mut all_retention_times = Vec::new();

    for line in lines.by_ref() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if let Some(rt_idx) = column_indices.get("retention_time") {
            if let Some(rt_str) = fields.get(*rt_idx) {
                if let Ok(rt) = rt_str.parse::<f32>() {
                    all_retention_times.push(rt);
                }
            }
        }
    }

    // Calculate global RT min/max (skip if no valid RTs found)
    let (rt_min, rt_max) = if !all_retention_times.is_empty() {
        let min = all_retention_times
            .iter()
            .fold(f32::INFINITY, |a, &b| a.min(b));
        let max = all_retention_times
            .iter()
            .fold(f32::NEG_INFINITY, |a, &b| a.max(b));
        (min, max)
    } else {
        (0.0, 1.0) // fallback if no RT data
    };
    let rt_range = rt_max - rt_min;

    // Reset reader for second pass
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let _ = lines.next(); // skip header

    // Group transitions by peptide AND charge
    let mut peptide_map = HashMap::new();
    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        let get_field = |name: &str| -> std::io::Result<&str> {
            column_indices
                .get(name)
                .and_then(|&idx| fields.get(idx))
                .map(|&s| s)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Missing field: {}", name),
                    )
                })
        };

        let sequence = get_field("sequence")?.to_string();
        let charge = match get_field("precursor_charge")?.parse() {
            Ok(c) => c,
            Err(_) => continue,
        };
        let precursor_mz = get_field("precursor_mz")?.parse::<f32>().unwrap_or(0.0_f32);
        let fragment_type = get_field("fragment_type")?;
        let series_number = get_field("fragment_series_number")?.parse().unwrap_or(0);
        let product_charge = column_indices
            .get("product_charge")
            .and_then(|&idx| fields.get(idx))
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);
        let intensity = get_field("intensity")?.parse().unwrap_or(0.0);
        let raw_rt = column_indices
            .get("retention_time")
            .and_then(|&idx| fields.get(idx))
            .and_then(|s| s.parse().ok());
        // Normalize RT to 0-100 scale (iRT-like)
        let normalized_rt = raw_rt.map(|rt: f32| {
            if rt_range > 0.0 {
                100.0 * (rt - rt_min) / rt_range
            } else {
                50.0 // fallback for zero range
            }
        });
        let ion_mobility = column_indices
            .get("ion_mobility")
            .and_then(|&idx| fields.get(idx).and_then(|s| s.parse::<f32>().ok()));
        let ccs = ion_mobility_to_ccs_bruker(
            ion_mobility.unwrap_or(0.0_f32) as f64,
            charge,
            precursor_mz as f64,
        );

        // Create unique key combining sequence and charge
        let peptide_key = (sequence.clone(), charge);

        // Get or create peptide entry
        let entry = peptide_map.entry(peptide_key).or_insert_with(|| {
            // Estimate peptide length from sequence (remove modifications first)
            let canonical_sequence = normalize_sequence_for_mod_parsing(&sequence);
            let naked_seq = remove_mass_shift(&canonical_sequence);
            let mods = get_modification_string(&canonical_sequence, &modifications);
            let mod_sites = get_modification_indices(&canonical_sequence);
            let peptide_len = naked_seq.len();

            // Initialize with empty intensity matrix
            PeptideData::new(
                &sequence,                                 // modified_sequence
                &naked_seq,                                // naked_sequence
                &mods,                                     // mods
                &mod_sites,                                // mod_sites
                Some(charge),                              // charge
                Some(precursor_mz),                        // precursor_mass
                Some(nce),                                 // nce
                Some(instrument),                          // instrument
                normalized_rt,                             // retention_time
                ion_mobility,                              // ion_mobility
                Some(ccs),                                 // ccs
                Some(vec![vec![0.0; 8]; peptide_len - 1]), // ms2_intensities
            )
        });

        // Process intensity into the matrix format
        if let Some(ref mut intensities) = entry.ms2_intensities {
            let col = match (fragment_type, product_charge) {
                ("b", 1) => 0, // b_z1
                ("b", 2) => 1, // b_z2
                ("y", 1) => 2, // y_z1
                ("y", 2) => 3, // y_z2
                _ => continue, // Skip unsupported fragment types/charges
            };

            let row = (series_number - 1) as usize; // Convert to zero-based index
            if row < intensities.len() && col < intensities[0].len() {
                intensities[row][col] = intensity;
            }
        }
    }

    // Apply min-max normalization to each peptide's intensities
    for peptide_data in peptide_map.values_mut() {
        if let Some(ref mut intensities) = peptide_data.ms2_intensities {
            let min_int = intensities
                .iter()
                .flatten()
                .fold(f32::INFINITY, |a, &b| a.min(b));
            let max_int = intensities
                .iter()
                .flatten()
                .fold(f32::NEG_INFINITY, |a, &b| a.max(b));
            let range = max_int - min_int;

            if range > 0.0 {
                for row in intensities.iter_mut() {
                    for val in row.iter_mut() {
                        *val = (*val - min_int) / range;
                    }
                }
            }
        }
    }

    Ok(peptide_map.into_values().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::get_test_file;

    #[test]
    fn test_read_peptide_data_from_tsv() {
        let test_file = get_test_file("fine_tune_data.tsv");
        let nce = 20;
        let instrument = "QE";

        let result = read_peptide_data_from_tsv(test_file, nce, instrument);
        assert!(result.is_ok());

        let peptide_data = result.unwrap();
        assert!(!peptide_data.is_empty());
        assert!(peptide_data.len() == 500);
    }

    #[test]
    fn test_modification_parsing_for_tuning_data() {
        let modifications = load_modifications().unwrap();

        let mass_shift = canonicalize_mass_shift_annotations("PEPT[79.9663]IDE");
        assert_eq!(mass_shift, "PEPT[+79.9663]IDE");
        assert_eq!(
            get_modification_string(&mass_shift, &modifications),
            "Phospho@T"
        );
        assert_eq!(get_modification_indices(&mass_shift), "4");

        let unimod = "(UniMod:1)M(UniMod:35)PEPTIDE";
        assert_eq!(
            get_modification_string(unimod, &modifications),
            "Acetyl@Protein_N-term;Oxidation@M"
        );
        assert_eq!(get_modification_indices(unimod), "0;1");

        let dotted_unimod = normalize_sequence_for_mod_parsing(".(UniMod:1)M(UniMod:35)PEPTIDE");
        assert_eq!(dotted_unimod, "(UniMod:1)M(UniMod:35)PEPTIDE");
        assert_eq!(
            get_modification_string(&dotted_unimod, &modifications),
            "Acetyl@Protein_N-term;Oxidation@M"
        );
        assert_eq!(get_modification_indices(&dotted_unimod), "0;1");

        let sage_mass_shift = normalize_sequence_for_mod_parsing("[+42.0106]-MPEPT[+79.9663]IDE");
        assert_eq!(sage_mass_shift, "[+42.0106]MPEPT[+79.9663]IDE");
        assert_eq!(
            get_modification_string(&sage_mass_shift, &modifications),
            "Acetyl@Protein_N-term;Phospho@T"
        );
        assert_eq!(get_modification_indices(&sage_mass_shift), "0;5");
    }
}
