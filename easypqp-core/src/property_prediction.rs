use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use core::f32;
use anyhow::Result;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use rustyms::fragment::FragmentKind;
use rustyms::FragmentationModel;
use rustyms::system::e;
use rustyms::system::usize::Charge;
use sage_core::database::IndexedDatabase;
use redeem_properties::models::model_interface::{PredictionValue, PredictionResult};
use redeem_properties::utils::logging::Progress;
use redeem_properties::{
    models::model_interface::DLModels,
    utils::peptdeep_utils::{
        ccs_to_mobility_bruker, get_modification_string,
        remove_mass_shift, ModificationMap,
    },
};

use crate::{is_allowed_fragment, select_model, InsilicoPQPSettings, PeptideProperties, PrecursorProperties, ProductProperties};


pub struct PropertyPrediction<'db, 'params> {
    pub db: &'db IndexedDatabase,
    pub insilico_settings: &'params InsilicoPQPSettings,
    pub dl_models: DLModels,
    pub modifications: HashMap<(String, Option<char>), ModificationMap>,
    pub batch_size: usize,
    pub predictions:
        Arc<Mutex<HashMap<(u32, u8), (Option<f32>, Option<f32>, Option<PredictionValue>)>>>,
}

impl<'db, 'params> PropertyPrediction<'db, 'params> {
    pub fn new(
        db: &'db IndexedDatabase,
        insilico_settings: &'params InsilicoPQPSettings,
        dl_models: DLModels,
        modifications: HashMap<(String, Option<char>), ModificationMap>,
        batch_size: usize,
    ) -> Self {
        PropertyPrediction {
            db,
            insilico_settings,
            dl_models,
            modifications,
            batch_size,
            predictions: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    pub fn predict_properties(&mut self) -> Result<Vec<PeptideProperties>> {
        // let predictions = Arc::clone(&self.predictions);

        let total = if self.db.peptides.len() > self.batch_size {
            self.db.peptides.len() / self.batch_size
        } else {
            1 // Ensure at least one batch if batch_size is larger than total length
        };

        let progress = Progress::new(total, "Predicting properties...");

        // Lock the model once to reduce lock contention
        let rt_model_guard = self
            .dl_models
            .rt_model
            .as_ref()
            .and_then(|model| model.lock().ok());

        let ccs_model_guard = self
            .dl_models
            .ccs_model
            .as_ref()
            .and_then(|model| model.lock().ok());

        let ms2_model_guard = self
            .dl_models
            .ms2_model
            .as_ref()
            .and_then(|model| model.lock().ok());

        // Use a Vec to collect thread-local HashMaps
        let predictions: Vec<_> = self.db.peptides
            .iter()
            .enumerate() // Attach global indices (0..N)
            .collect::<Vec<_>>() // Collect into Vec first
            .par_chunks(self.batch_size) // Use Rayon's parallel chunks
            .map(|batch| {
                // Prepare batched inputs
                let mut peptide_sequences = Vec::new();
                let mut mods = Vec::new();
                let mut mod_sites = Vec::new();
                let mut charges = Vec::new();
                let mut peptide_indices = Vec::new();
        
                for (global_idx, peptide) in batch {
                    for &charge in &self.insilico_settings.precursor_charge {
                        // let peptide = &std::str::from_utf8(&unique_peptide.sequence).unwrap();
                        let peptide_string = peptide.to_string();
                        let naked_peptide = remove_mass_shift(&peptide_string);
                        let naked_peptide = naked_peptide.trim_start_matches("-");
        
                        let modified_indices: Vec<String> = peptide
                            .modifications
                            .iter()
                            .enumerate()
                            .filter_map(|(index, &value)| {
                                if value > 0.0 {
                                    Some(index.to_string())
                                } else {
                                    None
                                }
                            })
                            .collect();
        
                        let mod_str = get_modification_string(&peptide.to_string(), &self.modifications);
                        let mod_site_str = modified_indices.join(";");
        
                        peptide_sequences.push(naked_peptide.to_string());
                        mods.push(mod_str);
                        mod_sites.push(mod_site_str);
                        charges.push(charge as i32);
                        peptide_indices.push((
                            *global_idx as u32, // Dereference the global_idx
                            charge as u8,
                        ));
                    }
                }

                // Perform predictions 
                let rt_preds: Option<PredictionResult> = rt_model_guard.as_ref().and_then(|model| {
                    model.predict(&peptide_sequences, &mods, &mod_sites).or_else(|err| {
                        eprintln!("Error predicting RT: {:?}", err);
                        Err(None::<PredictionResult>)
                    }).ok()
                });

                let ccs_preds: Option<PredictionResult> = ccs_model_guard.as_ref().and_then(|model| {
                    model.predict(
                        &peptide_sequences,
                        &mods,
                        &mod_sites,
                        charges.clone(),
                    ).or_else(|err| {
                        eprintln!("Error predicting CCS: {:?}", err);
                        Err(None::<PredictionResult>)
                    }).ok()
                });

                let ms2_preds: Option<PredictionResult> = ms2_model_guard.as_ref().and_then(|model| {
                    model.predict(
                        &peptide_sequences,
                        &mods,
                        &mod_sites,
                        charges.clone(),
                        vec![
                            self.dl_models.params.as_ref().unwrap().nce as i32;
                            peptide_sequences.len()
                        ],
                        vec![
                            self.dl_models
                                .params
                                .as_ref()
                                .unwrap()
                                .instrument
                                .clone();
                            peptide_sequences.len()
                        ],
                    ).or_else(|err| {
                        eprintln!("Error predicting MS2: {:?}", err);
                        Err(None::<PredictionResult>)
                    }).ok()
                });

                // Create a thread-local HashMap to store predictions
                let mut thread_predictions = HashMap::new();

                // Insert predictions into the thread-local HashMap
                for (i, peptide_idx) in peptide_indices.iter().enumerate() {
                    let rt_pred = rt_preds.as_ref().map(|v| v.get_prediction_entry(i)[0]);
                    let ccs_pred = ccs_preds.as_ref().map(|v| v.get_prediction_entry(i)[0]);
                    let ms2_pred = ms2_preds.as_ref().map(|v| v.get_prediction_entry(i).clone());
                    thread_predictions.insert(*peptide_idx, (Some(rt_pred.unwrap_or_default()), Some(ccs_pred.unwrap_or_default()), Some(ms2_pred.unwrap())));
                }

                progress.inc(); // Thread-safe progress increment

                thread_predictions // Return the thread-local HashMap
            })
            .collect(); // Collect all thread-local HashMaps into a Vec

        // Combine all thread-local HashMaps into the final shared HashMap
        {
            let mut final_predictions = self.predictions.lock().unwrap();
            for thread_predictions in predictions {
                for (key, value) in thread_predictions {
                    final_predictions.insert(key, value);
                }
            }
        } // Scope locking to release the lock

        progress.finish(); // Ensure all updates are processed before exiting

        // Clone the predictions HashMap to avoid locking in the loop
        let preds_clone: HashMap<(u32, u8), (Option<f32>, Option<f32>, Option<PredictionValue>)> = {
            let preds = self.predictions.lock().unwrap();
            preds.clone()
        };

        // Assign predictions
        let mut assays = Vec::new();
        let total = preds_clone.len();
        let progress = Progress::new(
            total,
            "Structuring and storing predictions...",
        );
        // Iterate over preds_clone
        for ((peptide_idx, charge), (rt_pred, ccs_pred, ms2_pred)) in preds_clone.iter() {
            // Check if the predictions are None
            if rt_pred.is_none() && ccs_pred.is_none() && ms2_pred.is_none() {
                continue; // Skip to the next iteration
            }

            let precursor_mz = (self.db.peptides[*peptide_idx as usize].monoisotopic / (*charge as f32));

            log::trace!(
                "Peptide: {}, m/z: {}, Charge: {}, RT: {:?}, CCS: {:?}",
                self.db.peptides[*peptide_idx as usize].to_string(),
                precursor_mz,
                charge,
                rt_pred,
                ccs_pred
            );

            let peptide = rustyms::peptidoform::PeptidoformIon::pro_forma(&self.db.peptides[*peptide_idx as usize].to_string(), None)?;
            let model = select_model(&self.insilico_settings.fragmentation_model, FragmentationModel::cid_hcd());

            let fragments = peptide.generate_theoretical_fragments(Charge::new::<e>(self.insilico_settings.max_fragment_charge), model);

            // filter fragments for ion kind b or y
            let fragments: Vec<_> = fragments
                .into_iter()
                .filter(|fragment| {
                    is_allowed_fragment(
                        fragment.ion.kind(),
                        &self.insilico_settings.allowed_fragment_types
                    )
                })
                .collect();

            // filter for fragments with neutral_loss vector length = 0, current MS2 prediction model does not support neutral loss
            let fragments: Vec<_> = fragments
                .into_iter()
                .filter(|fragment| fragment.neutral_loss.len() == 0)
                .collect();
            
            // println!("Fragments 0 {:?}{:?}^{:?} mz: {:?}", fragments[0].ion.kind(), fragments[0].ion.position().unwrap().series_number, fragments[0].charge.value, fragments[0].mz(rustyms::MassMode::Monoisotopic).unwrap().value);

            let mut ion_type = Vec::new();
            let mut ion_ordinal = Vec::new();
            let mut charges = Vec::new();
            let mut product_mz = Vec::new();
            let mut intensities = Vec::new();

            for fragment in fragments {
                let ion = &fragment.ion;
                let mz = fragment.mz(rustyms::MassMode::Monoisotopic).unwrap().value;

                let (row, col) = match (ion.kind(), fragment.charge.value) {
                    (FragmentKind::b, 1) => (ion.position().unwrap().series_number as u8 - 1, 0), // b_z1
                    (FragmentKind::b, 2) => (ion.position().unwrap().series_number as u8 - 1, 1), // b_z2
                    (FragmentKind::y, 1) => (ion.position().unwrap().series_number as u8 - 1, 2), // y_z1
                    (FragmentKind::y, 2) => (ion.position().unwrap().series_number as u8 - 1, 3), // y_z2
                    _ => continue, // Skip unsupported fragment types or charges
                };

                // Check if ms2_pred is not None and get the intensity, else return 0.0
                let intensity = if let Some(ms2_pred) = ms2_pred {
                    *ms2_pred.get(row as usize, col as usize).unwrap()
                } else {
                    0.0 // Default value if intensity is not available
                };

                ion_type.push(ion.kind());
                ion_ordinal.push(ion.position().unwrap().series_number as u8);
                charges.push(fragment.charge.value as u8);
                product_mz.push(mz);
                intensities.push(intensity);
            }
            let product = ProductProperties {
                peptide_index: *peptide_idx,
                ion_type,
                ion_ordinal,
                charge: charges,
                product_mz: product_mz,
                intensity: intensities,
            };

            
            let assay = PeptideProperties {
                peptide_index: *peptide_idx,
                retention_time: rt_pred.unwrap_or_default(),
                precursor: PrecursorProperties {
                    peptide_index: *peptide_idx,
                    charge: *charge,
                    precursor_mz: precursor_mz,
                    ion_mobility: ccs_pred.map(|pred| ccs_to_mobility_bruker(
                        pred as f64,
                        *charge as f64,
                        precursor_mz as f64,
                    )).unwrap_or(f64::NAN),
                },
                product: product,
            };

            assays.push(assay);

            progress.inc(); // Thread-safe progress increment
        }
        progress.finish(); // Ensure all updates are processed before exiting

        Ok(assays)
    }
}