use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use redeem_properties::utils::data_handling::PeptideData;
use redeem_properties::utils::peptdeep_utils::{ion_mobility_to_ccs_bruker, remove_mass_shift};

pub fn read_peptide_data_from_tsv<P: AsRef<Path>>(
    path: P,
    nce: i32,
    instrument: &str,
) -> std::io::Result<Vec<PeptideData>> {
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
            s if s.contains("precursormz") || s.contains("precursor_mz")  => {
                column_indices.insert("precursor_mz", idx);
            }
            s if s.contains("precursorcharge") || s.contains("charge") || s.contains("precursor_charge") => {
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
            s if s.contains("normalizedretentiontime") || s.contains("rt")  || s.contains("retention_time") => {
                column_indices.insert("retention_time", idx);
            }
            s if s.contains("precursorionmobility") || s.contains("im") || s.contains("ion_mobility") => {
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
    let required_columns = ["sequence", "precursor_charge", "intensity", "retention_time"];
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
    let mut line_count = 0;
    
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
        line_count += 1;
    }

    // Calculate global RT min/max (skip if no valid RTs found)
    let (rt_min, rt_max) = if !all_retention_times.is_empty() {
        let min = all_retention_times.iter().fold(f32::INFINITY, |a, &b| a.min(b));
        let max = all_retention_times.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
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
        let precursor_mz = get_field("precursor_mz")?.parse().unwrap_or(0.0);
        let fragment_type = get_field("fragment_type")?;
        let series_number = get_field("fragment_series_number")?.parse().unwrap_or(0);
        let product_charge = column_indices.get("product_charge")
            .and_then(|&idx| fields.get(idx))
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);
        let intensity = get_field("intensity")?.parse().unwrap_or(0.0);
        let raw_rt = column_indices.get("retention_time")
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
        let ion_mobility = column_indices.get("ion_mobility").and_then(|&idx| fields.get(idx).and_then(|s| s.parse().ok()));
        let ccs = ion_mobility_to_ccs_bruker(ion_mobility.unwrap_or(0.0), charge, precursor_mz);

        // Create unique key combining sequence and charge
        let peptide_key = (sequence.clone(), charge);

        // Get or create peptide entry
        let entry = peptide_map.entry(peptide_key).or_insert_with(|| {
            // Estimate peptide length from sequence (remove modifications first)
            let naked_seq = remove_mass_shift(&sequence);
            let peptide_len = naked_seq.len();
            
            // Initialize with empty intensity matrix
            PeptideData::new(
                &sequence,
                Some(charge),
                Some(nce),
                Some(instrument),
                normalized_rt,
                Some(ccs),
                Some(vec![vec![0.0; 8]; peptide_len - 1]), // Same format as experimental_to_predicted_format
            )
        });

        // Process intensity into the matrix format
        if let Some(ref mut intensities) = entry.ms2_intensities {
            let col = match (fragment_type, product_charge) {
                ("b", 1) => 0,  // b_z1
                ("b", 2) => 1,  // b_z2
                ("y", 1) => 2,  // y_z1
                ("y", 2) => 3,  // y_z2
                _ => continue,   // Skip unsupported fragment types/charges
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
            let min_int = intensities.iter().flatten().fold(f32::INFINITY, |a, &b| a.min(b));
            let max_int = intensities.iter().flatten().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
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