use anyhow::Result;
use core::f32;
use itertools::multiunzip;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use redeem_properties::models::model_interface::PredictionValue;
use redeem_properties::utils::logging::Progress;
use redeem_properties::{
    models::model_interface::DLModels,
    utils::peptdeep_utils::{
        ccs_to_mobility_bruker, get_modification_string, remove_mass_shift, ModificationMap,
    },
};
use rustyms::fragment::FragmentKind;
use rustyms::system::e;
use rustyms::system::usize::Charge;
use rustyms::FragmentationModel;
use sage_core::peptide::Peptide;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use crate::{
    is_allowed_fragment, select_model, InsilicoPQPSettings, PeptideProperties, PrecursorProperties,
    ProductProperties,
};

pub struct PropertyPrediction<'db, 'params> {
    pub peptides: &'db Vec<Peptide>,
    pub insilico_settings: &'params InsilicoPQPSettings,
    pub dl_models: DLModels,
    pub modifications: HashMap<(String, Option<char>), ModificationMap>,
    pub batch_size: usize,
    pub predictions:
        Arc<Mutex<HashMap<(u32, u8), (Option<f32>, Option<f32>, Option<PredictionValue>)>>>,
}

impl<'db, 'params> PropertyPrediction<'db, 'params> {
    pub fn new(
        peptides: &'db Vec<Peptide>,
        insilico_settings: &'params InsilicoPQPSettings,
        dl_models: DLModels,
        modifications: HashMap<(String, Option<char>), ModificationMap>,
        batch_size: usize,
    ) -> Self {
        PropertyPrediction {
            peptides,
            insilico_settings,
            dl_models,
            modifications,
            batch_size,
            predictions: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    pub fn predict_properties(&mut self) -> Result<Vec<PeptideProperties>> {
        let rt_model = self.dl_models.rt_model.as_ref().cloned();
        let ccs_model = self.dl_models.ccs_model.as_ref().cloned();
        let ms2_model = self.dl_models.ms2_model.as_ref().cloned();

        // Step 1: Parallelized peptide + charge expansion with progress
        let total_peptides = self.peptides.len();
        let progress = Progress::new(total_peptides, "Expanding peptides...");
        let (send, recv) = crossbeam_channel::bounded(1024);

        // Spawn thread for progress bar
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let expanded: Vec<_> = self
            .peptides
            .par_iter()
            .enumerate()
            .map_init(
                || send.clone(),
                |sender, (idx, peptide)| {
                    let peptide_string = peptide.to_string();
                    let naked_peptide = remove_mass_shift(&peptide_string)
                        .trim_start_matches("-")
                        .to_string();
                    let mod_str = get_modification_string(&peptide_string, &self.modifications);
                    let mod_site_str = peptide
                        .modifications
                        .iter()
                        .enumerate()
                        .filter_map(|(i, &v)| if v > 0.0 { Some(i.to_string()) } else { None })
                        .collect::<Vec<_>>()
                        .join(";");

                    let _ = sender.send(()); // update progress bar

                    self.insilico_settings
                        .precursor_charge
                        .iter()
                        .map(move |&charge| {
                            (
                                naked_peptide.clone(),
                                mod_str.clone(),
                                mod_site_str.clone(),
                                charge as i32,
                                (idx as u32, charge as u8),
                            )
                        })
                        .collect::<Vec<_>>()
                },
            )
            .flatten()
            .collect();

        drop(send);
        progress_thread.join().unwrap();

        let (peptide_sequences, mods, mod_sites, charges, peptide_indices): (
            Vec<String>,
            Vec<String>,
            Vec<String>,
            Vec<i32>,
            Vec<(u32, u8)>,
        ) = multiunzip(expanded);

        // Step 2: Create chunks for prediction
        let total_batches = (peptide_sequences.len() + self.batch_size - 1) / self.batch_size;
        let progress = Progress::new(total_batches, "Predicting properties...");
        let (send, recv) = crossbeam_channel::bounded(1024);

        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let batches: Vec<_> = (0..peptide_sequences.len())
            .step_by(self.batch_size)
            .map(|start| {
                let end = (start + self.batch_size).min(peptide_sequences.len());
                (
                    peptide_sequences[start..end].to_vec(),
                    mods[start..end].to_vec(),
                    mod_sites[start..end].to_vec(),
                    charges[start..end].to_vec(),
                    peptide_indices[start..end].to_vec(),
                )
            })
            .collect();

        let params = self
            .dl_models
            .params
            .as_ref()
            .expect("DL model parameters must be present")
            .clone(); // clone once here

        let predictions: Vec<_> = batches
            .par_iter()
            .map_init(
                || {
                    send.clone()
                },
                |sender, (seqs, mods, sites, chgs, indices)| {
                    let start_time = std::time::Instant::now();
                    let rt_preds = rt_model
                        .as_ref()
                        .and_then(|m| m.predict(seqs, mods, sites).ok());
                    log::debug!("RT prediction time: {:?}", start_time.elapsed());
                    let start_time = std::time::Instant::now();
                    let ccs_preds = ccs_model
                        .as_ref()
                        .and_then(|m| m.predict(seqs, mods, sites, chgs.clone()).ok());
                    log::debug!("CCS prediction time: {:?}", start_time.elapsed());
                    let start_time = std::time::Instant::now();
                    let ms2_preds = ms2_model.as_ref().and_then(|m| {
                        m.predict(
                            seqs,
                            mods,
                            sites,
                            chgs.clone(),
                            vec![params.nce as i32; seqs.len()],
                            vec![params.instrument.clone(); seqs.len()],
                        )
                        .ok()
                    });
                    log::debug!("MS2 prediction time: {:?}", start_time.elapsed());

                    let start_time = std::time::Instant::now();
                    let mut result = Vec::new();
                    for (i, (peptide_idx, charge)) in indices.iter().enumerate() {
                        result.push((
                            (*peptide_idx, *charge),
                            (
                                rt_preds.as_ref().map(|v| v.get_prediction_entry(i)[0]),
                                ccs_preds.as_ref().map(|v| v.get_prediction_entry(i)[0]),
                                ms2_preds
                                    .as_ref()
                                    .map(|v| v.get_prediction_entry(i).clone()),
                            ),
                        ));
                    }
                    log::debug!(
                        "Collecting prediction results time: {:?}",
                        start_time.elapsed()
                    );

                    let _ = sender.send(());
                    result
                },
            )
            .flatten()
            .collect();

        drop(send);
        progress_thread.join().unwrap();

        {
            let mut final_predictions = self.predictions.lock().unwrap();
            for (key, value) in predictions {
                final_predictions.insert(key, value);
            }
        }

        // Step 3: Build assays from predictions
        let preds_clone = {
            let preds = self.predictions.lock().unwrap();
            preds.clone()
        };

        let total_assays = preds_clone.len();
        let progress = Progress::new(total_assays, "Creating PQP assays...");
        let (send, recv) = crossbeam_channel::bounded(1024);

        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let model = select_model(
            &self.insilico_settings.fragmentation_model,
            FragmentationModel::cid_hcd(),
        );

        let assays: Vec<_> = preds_clone
            .par_iter()
            .map_init(
                || (model, send.clone()),
                |(model, sender), ((peptide_idx, charge), (rt_pred, ccs_pred, ms2_pred))| {
                    let peptide = &self.peptides[*peptide_idx as usize];
                    let precursor_mz = peptide.monoisotopic / (*charge as f32);

                    let peptidoform = rustyms::peptidoform::PeptidoformIon::pro_forma(
                        &peptide.to_string(),
                        None,
                    )?;

                    let fragments = peptidoform.generate_theoretical_fragments(
                        Charge::new::<e>(self.insilico_settings.max_fragment_charge),
                        *model,
                    );

                    let mut product = ProductProperties {
                        peptide_index: *peptide_idx,
                        ion_type: Vec::new(),
                        ion_ordinal: Vec::new(),
                        charge: Vec::new(),
                        product_mz: Vec::new(),
                        intensity: Vec::new(),
                    };

                    for fragment in fragments {
                        if !is_allowed_fragment(
                            fragment.ion.kind(),
                            &self.insilico_settings.allowed_fragment_types,
                        ) || !fragment.neutral_loss.is_empty()
                        {
                            continue;
                        }

                        let ion = &fragment.ion;
                        let mz = fragment.mz(rustyms::MassMode::Monoisotopic).unwrap().value;

                        let (row, col) = match (ion.kind(), fragment.charge.value) {
                            (FragmentKind::b, 1) => {
                                (ion.position().unwrap().series_number as u8 - 1, 0)
                            }
                            (FragmentKind::b, 2) => {
                                (ion.position().unwrap().series_number as u8 - 1, 1)
                            }
                            (FragmentKind::y, 1) => {
                                (ion.position().unwrap().series_number as u8 - 1, 2)
                            }
                            (FragmentKind::y, 2) => {
                                (ion.position().unwrap().series_number as u8 - 1, 3)
                            }
                            _ => continue,
                        };

                        let intensity = ms2_pred
                            .as_ref()
                            .and_then(|pred| pred.get(row as usize, col as usize))
                            .copied()
                            .unwrap_or(0.0);

                        product.ion_type.push(ion.kind());
                        product
                            .ion_ordinal
                            .push(ion.position().unwrap().series_number as u8);
                        product.charge.push(fragment.charge.value as u8);
                        product.product_mz.push(mz);
                        product.intensity.push(intensity);
                    }

                    let _ = sender.send(());

                    Ok(PeptideProperties {
                        peptide_index: *peptide_idx,
                        retention_time: rt_pred.unwrap_or_default(),
                        precursor: PrecursorProperties {
                            peptide_index: *peptide_idx,
                            charge: *charge,
                            precursor_mz,
                            ion_mobility: ccs_pred
                                .map(|ccs| {
                                    ccs_to_mobility_bruker(
                                        ccs as f64,
                                        *charge as f64,
                                        precursor_mz as f64,
                                    )
                                })
                                .unwrap_or(f64::NAN),
                        },
                        product,
                    })
                },
            )
            .collect::<Result<Vec<_>>>()?;

        drop(send);
        progress_thread.join().unwrap();

        Ok(assays)
    }
}
