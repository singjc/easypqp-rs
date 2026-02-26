use anyhow::Result;
use chrono::Utc;
use core::f32;
use itertools::multiunzip;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use redeem_properties::models::model_interface::PredictionValue;
use redeem_properties::utils::logging::Progress;
use redeem_properties::{
    models::model_interface::DLModels,
    utils::peptdeep_utils::{
        ccs_to_mobility_bruker, ModificationMap,
    },
};
use rustyms::fragment::FragmentKind;
use rustyms::system::e;
use rustyms::system::usize::Charge;
use rustyms::{Chemical, FragmentationModel, MolecularCharge};
use sage_core::peptide::Peptide;
use std::collections::HashMap;
use std::sync::Arc;

use crate::{
    is_allowed_fragment, select_model, InsilicoPQPSettings, PeptideProperties, PrecursorProperties,
    ProductProperties,
};


/// A simplified feature representation of a peptide used for model input.
///
/// Contains the stripped amino acid sequence (`naked_sequence`), a
/// stringified representation of modifications (`mod_string`), and a
/// semicolon-separated list of modification site indices (`mod_sites`).
#[derive(Clone)]
pub struct PeptideFeatures {
    /// Amino acid sequence with mass shifts removed (e.g., "ACDEFGHIK").
    pub naked_sequence: String,

    /// Modification string formatted for input into deep learning models (e.g., "Oxidation@M").
    pub mod_string: String,

    /// Semicolon-separated modification site indices (e.g., "3;6").
    pub mod_sites: String,
}


/// Holds configuration and data required for peptide property prediction.
///
/// This struct encapsulates all components necessary to predict
/// retention time (RT), collisional cross-section (CCS), and MS2 fragment
/// intensities using deep learning models for a given list of peptides.
pub struct PropertyPrediction<'db, 'params> {
    /// Slice of peptides to predict properties for.
    pub peptides: &'db [Peptide],

    /// Settings that control prediction and fragmentation behavior.
    pub insilico_settings: &'params InsilicoPQPSettings,

    /// Pre-loaded deep learning models for RT, CCS, and MS2 prediction.
    pub dl_models: DLModels,

    /// Mapping of modification identifiers used to convert peptide mods to model-compatible strings.
    pub modifications: HashMap<(String, Option<char>), ModificationMap>,

    /// Number of peptide inputs to include per model batch.
    pub batch_size: usize,
}


impl<'db, 'params> PropertyPrediction<'db, 'params> {
    /// Constructs a new `PropertyPrediction` instance for predicting peptide properties.
    ///
    /// # Arguments
    /// * `peptides` - A slice of peptide structures to process.
    /// * `insilico_settings` - Configuration parameters for fragmentation and prediction.
    /// * `dl_models` - Loaded prediction models for RT, CCS, and MS2.
    /// * `modifications` - Map of modification names and sites to use in feature encoding.
    /// * `batch_size` - Number of inputs per prediction batch.
    ///
    /// # Returns
    /// A new `PropertyPrediction` instance with all fields initialized.
    pub fn new(
        peptides: &'db [Peptide],
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
        }
    }

    /// Predicts peptide properties (RT, CCS, MS2) and builds PQP assays.
    ///
    /// This is the main entry point for generating predicted peptide assays.
    /// It performs the following steps:
    /// 1. Extracts sequence and modification features from peptides.
    /// 2. Predicts retention times (RT) in batches.
    /// 3. Expands peptides across precursor charge states.
    /// 4. Predicts CCS and MS2 intensities in batches.
    /// 5. Builds assays from the combined predictions.
    ///
    /// Returns a vector of `PeptideProperties` ready for export.
    pub fn predict_properties(&mut self) -> Result<Vec<PeptideProperties>> {
        // Step 1: Extract features once per peptide
        let peptide_features = self.extract_peptide_features();

        // Step 2: Predict RT in batches
        let rt_predictions = self.predict_rt(&peptide_features);

        // Step 3: Expand peptides by charge
        let expanded = self.expand_peptides_by_charge(&peptide_features);

        // Step 4: Predict CCS + MS2
        let predictions = self.predict_ccs_and_ms2(expanded, rt_predictions);

        // Step 4: Build assays from predictions
        let assays = self.build_assays(&predictions)?;

        Ok(assays)
    }

    /// Extracts RT-relevant features from each peptide sequence using fast u8-based access.
    ///
    /// This version avoids expensive string conversions by operating directly on
    /// the `Peptide` struct's internal representation (`sequence` and `modifications`).
    ///
    /// Returns a vector of `PeptideFeatures` containing:
    /// - the naked amino acid sequence,
    /// - a model-compatible modification string,
    /// - and a semicolon-delimited list of modified site indices.
    fn extract_peptide_features(&self) -> Vec<PeptideFeatures> {
        let total_peptides = self.peptides.len();
        let timestamp = Utc::now()
            .format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]")
            .to_string();
        let description = format!("{} Preparing features...", timestamp);
        let progress = Progress::new(total_peptides, &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let features: Vec<PeptideFeatures> = self
            .peptides
            .par_iter()
            .map_init(
                || send.clone(),
                |sender, peptide| {
                    let _ = sender.send(());
                    peptide_feature_from(peptide, &self.modifications)
                },
            )
            .collect();

        drop(send);
        progress_thread.join().unwrap();
        features
    }

    /// Predicts retention time (RT) in batches using the RT model.
    ///
    /// Accepts a slice of `PeptideFeatures` and performs batch inference using
    /// the configured RT model. Progress is logged with timestamps.
    ///
    /// Returns a shared vector (`Arc<Vec<...>>`) of optional RT predictions,
    /// indexed by peptide.
    fn predict_rt(&self, peptide_features: &[PeptideFeatures]) -> Arc<Vec<Option<f32>>> {
        let batched_inputs: Vec<(Vec<Arc<[u8]>>, Vec<Arc<[u8]>>, Vec<Arc<[u8]>>)> = peptide_features
            .chunks(self.batch_size)
            .map(|chunk| {
                (
                    chunk
                        .iter()
                        .map(|f| Arc::from(f.naked_sequence.as_bytes().to_vec().into_boxed_slice()))
                        .collect(),
                    chunk
                        .iter()
                        .map(|f| Arc::from(f.mod_string.as_bytes().to_vec().into_boxed_slice()))
                        .collect(),
                    chunk
                        .iter()
                        .map(|f| Arc::from(f.mod_sites.as_bytes().to_vec().into_boxed_slice()))
                        .collect(),
                )
            })
            .collect();

        let total_batches = batched_inputs.len();
        let timestamp = Utc::now()
            .format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]")
            .to_string();
        let description = format!("{} Predicting RT...", timestamp);
        let progress = Progress::new(total_batches, &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let rt_model = self.dl_models.rt_model.as_ref().cloned();

        let rt_predictions: Arc<Vec<Option<f32>>> = Arc::new(
            batched_inputs
                .par_iter()
                .map_init(
                    || send.clone(),
                    |sender, (seqs, mods, sites)| {
                            let preds = rt_model
                                .as_ref()
                                .and_then(|m| m.predict(seqs.as_slice(), mods.as_slice(), sites.as_slice()).ok());
                        let results = (0..seqs.len())
                            .map(|i| preds.as_ref().map(|p| p.get_prediction_entry(i)[0]))
                            .collect::<Vec<_>>();
                        let _ = sender.send(());
                        results
                    },
                )
                .flatten()
                .collect(),
        );
        drop(send);
        progress_thread.join().unwrap();
        rt_predictions
    }

    /// Expands each peptide to all specified precursor charge states.
    ///
    /// For each peptide feature, it generates a new entry for each charge
    /// defined in `insilico_settings.precursor_charge`.
    ///
    /// Returns a flattened vector of tuples containing sequence info, charge,
    /// and the (peptide_index, charge) key used to link back predictions.
    fn expand_peptides_by_charge(
        &self,
        peptide_features: &[PeptideFeatures],
    ) -> Vec<(String, String, String, i32, (u32, u8))> {
        let timestamp = Utc::now()
            .format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]")
            .to_string();
        let description = format!("{} Expanding features...", timestamp);
        let progress = Progress::new(self.peptides.len(), &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let expanded: Vec<_> = peptide_features
            .par_iter()
            .enumerate()
            .map_init(
                || send.clone(),
                |sender, (idx, f)| {
                    let _ = sender.send(());
                    self.insilico_settings
                        .precursor_charge
                        .iter()
                        .map(move |&charge| {
                            (
                                f.naked_sequence.clone(),
                                f.mod_string.clone(),
                                f.mod_sites.clone(),
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
        expanded
    }

    /// Predicts CCS and MS2 intensities in batches for expanded peptide-charge pairs.
    ///
    /// Performs batched inference using the CCS and MS2 models.
    /// RT values are retrieved using the original peptide indices.
    ///
    /// Returns a vector of predictions keyed by (peptide_index, charge).
    /// Each value is a tuple of optional RT, CCS, and MS2 predictions.
    fn predict_ccs_and_ms2(
        &self,
        expanded: Vec<(String, String, String, i32, (u32, u8))>,
        rt_predictions: Arc<Vec<Option<f32>>>,
    ) -> Vec<(
        (u32, u8),
        (Option<f32>, Option<f32>, Option<PredictionValue>),
    )> {
        let (seqs, mods, sites, charges, peptide_indices): (
            Vec<String>,
            Vec<String>,
            Vec<String>,
            Vec<i32>,
            Vec<(u32, u8)>,
        ) = multiunzip(expanded);

        let param = self
            .dl_models
            .params
            .as_ref()
            .expect("DL model parameters must be present")
            .clone();

        let batches: Vec<_> = (0..seqs.len())
            .step_by(self.batch_size)
            .map(|start| {
                let end = (start + self.batch_size).min(seqs.len());
                (
                    seqs[start..end].to_vec(),
                    mods[start..end].to_vec(),
                    sites[start..end].to_vec(),
                    charges[start..end].to_vec(),
                    peptide_indices[start..end].to_vec(),
                )
            })
            .collect();

        let total_batches = batches.len();
        let timestamp = Utc::now()
            .format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]")
            .to_string();
        let description = format!("{} Predicting CCS and MS2...", timestamp);
        let progress = Progress::new(total_batches, &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let ccs_model = self.dl_models.ccs_model.as_ref().cloned();
        let ms2_model = self.dl_models.ms2_model.as_ref().cloned();

        let predictions: Vec<_> = batches
            .par_iter()
            .map_init(
                || (send.clone(), Arc::clone(&rt_predictions)),
                |(sender, rt_preds), (seqs, mods, sites, chgs, indices)| {
                    // Convert string inputs to Arc<[u8]> for redeem API
                    let seqs_arc: Vec<Arc<[u8]>> = seqs
                        .iter()
                        .map(|s| Arc::from(s.as_bytes().to_vec().into_boxed_slice()))
                        .collect();
                    let mods_arc: Vec<Arc<[u8]>> = mods
                        .iter()
                        .map(|s| Arc::from(s.as_bytes().to_vec().into_boxed_slice()))
                        .collect();
                    let sites_arc: Vec<Arc<[u8]>> = sites
                        .iter()
                        .map(|s| Arc::from(s.as_bytes().to_vec().into_boxed_slice()))
                        .collect();

                    let ccs_preds = ccs_model.as_ref().and_then(|m| {
                        m.predict(
                            seqs_arc.as_slice(),
                            mods_arc.as_slice(),
                            sites_arc.as_slice(),
                            chgs.clone(),
                        )
                        .ok()
                    });

                    let instrument_arc: Arc<[u8]> = Arc::from(param.instrument.as_bytes().to_vec().into_boxed_slice());

                    let ms2_preds = ms2_model.as_ref().and_then(|m| {
                        m.predict(
                            seqs_arc.as_slice(),
                            mods_arc.as_slice(),
                            sites_arc.as_slice(),
                            chgs.clone(),
                            vec![param.nce as i32; seqs_arc.len()],
                            vec![Some(instrument_arc.clone()); seqs_arc.len()],
                        )
                        .ok()
                    });

                    let result: Vec<_> = indices
                        .iter()
                        .enumerate()
                        .map(|(i, (peptide_idx, charge))| {
                            (
                                (*peptide_idx, *charge),
                                (
                                    rt_preds.get(*peptide_idx as usize).copied().flatten(),
                                    ccs_preds.as_ref().map(|p| p.get_prediction_entry(i)[0]),
                                    ms2_preds
                                        .as_ref()
                                        .map(|p| p.get_prediction_entry(i).clone()),
                                ),
                            )
                        })
                        .collect();
                    let _ = sender.send(());
                    result
                },
            )
            .flatten()
            .collect();

        drop(send);
        progress_thread.join().unwrap();
        predictions
    }

    /// Constructs PQP assays from the full set of predicted properties.
    ///
    /// For each (peptide_index, charge) entry, it reconstructs theoretical
    /// fragments and associates MS2 intensities based on model outputs.
    /// Only allowed fragment types are used.
    ///
    /// Returns a vector of `PeptideProperties` representing the final assay library.
    fn build_assays(
        &self,
        predictions: &Vec<(
            (u32, u8),
            (Option<f32>, Option<f32>, Option<PredictionValue>),
        )>,
    ) -> Result<Vec<PeptideProperties>> {
        let total_assays = predictions.len();
        let timestamp = Utc::now()
            .format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]")
            .to_string();
        let description = format!("{} Creating PQP assays...", timestamp);
        let progress = Progress::new(total_assays, &description);
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

        let assays: Vec<_> = predictions
            .par_iter()
            .map_init(
                || (model, send.clone()),
                |(model, sender), ((peptide_idx, charge), (rt_pred, ccs_pred, ms2_pred))| {
                    let peptide = &self.peptides[*peptide_idx as usize];

                    let peptidoform = rustyms::peptidoform::PeptidoformIon::pro_forma(
                        &peptide.to_string(),
                        None,
                    )?;

                    // Compute precursor m/z using rustyms:
                    // neutral_formula + z protons, then divide by charge
                    let neutral_formula = peptidoform.formulas().first().cloned().unwrap_or_default();
                    let proton_formula = MolecularCharge::proton(*charge as isize).formula();
                    let charged_mass = (neutral_formula + &proton_formula)
                        .monoisotopic_mass();
                    let precursor_mz = (charged_mass
                        / rustyms::system::f64::Charge::new::<e>(*charge as f64))
                        .value as f32;

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
                                .unwrap_or(-1.0),
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


/// Constructs `PeptideFeatures` for a single peptide by extracting:
/// - the unmodified sequence as a `String`,
/// - the modification string (e.g., "Carbamidomethyl@C"),
/// - and a semicolon-separated list of modified site indices.
///
/// # Example
/// ```rust
/// use std::sync::Arc;
/// use std::collections::HashMap;
/// use sage_core::peptide::Peptide;
/// use sage_core::enzyme::Position;
/// use redeem_properties::utils::peptdeep_utils::{load_modifications, ModificationMap};
/// use easypqp_core::property_prediction::peptide_feature_from;
/// let peptide = Peptide {
///     decoy: false,
///     sequence: Arc::from(b"ACDEFGHIK".to_vec().into_boxed_slice()),
///     modifications: vec![0.0, 57.02146, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
///     nterm: None,
///     cterm: None,
///     monoisotopic: 0.0,
///     missed_cleavages: 0,
///     semi_enzymatic: false,
///     position: Position::Internal,
///     proteins: vec![Arc::from("P12345")],
/// };
/// let mod_map = load_modifications().unwrap();
/// let features = peptide_feature_from(&peptide, &mod_map);
/// assert_eq!(features.naked_sequence, "ACDEFGHIK");
/// assert_eq!(features.mod_string, "Carbamidomethyl@C");
/// assert_eq!(features.mod_sites, "1");
/// ```
pub fn peptide_feature_from(
    peptide: &Peptide,
    mod_map: &HashMap<(String, Option<char>), ModificationMap>,
) -> PeptideFeatures {
    let naked_sequence: String = peptide.sequence.iter().map(|&b| b as char).collect();

    let mod_string = peptide
        .sequence
        .iter()
        .zip(peptide.modifications.iter())
        .enumerate()
        .filter_map(|(_, (&aa, &mass))| {
            if mass != 0.0 {
                let key = (format!("{:.4}", mass), Some(aa as char));
                mod_map.get(&key).map(|m| m.name.clone())
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
        .join(";");

    let mod_sites = peptide
        .modifications
        .iter()
        .enumerate()
        .filter_map(|(i, &mass)| (mass > 0.0).then(|| i.to_string()))
        .collect::<Vec<_>>()
        .join(";");

    PeptideFeatures {
        naked_sequence,
        mod_string,
        mod_sites,
    }
}