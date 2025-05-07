
use anyhow::{Context, Result};
use easypqp_core::{property_prediction::PropertyPrediction, tuning_data::read_peptide_data_from_tsv, PeptideProperties};
use log::{info, warn};
use sage_core::{database, peptide::Peptide};
use std::{path::Path, sync::{Arc, Mutex}, time::Instant};

use crate::{input::InsilicoPQP, output::write_assays_to_tsv};


// use easypqp_core::property_prediction::{extract_unique_peptide_info, format_feature_for_fine_tuning, PropertyPrediction};
use redeem_properties::{models::{
    ccs_model::{load_collision_cross_section_model, CCSModelWrapper}, model_interface::DLModels, ms2_model::{load_ms2_model, MS2ModelWrapper}, rt_model::{load_retention_time_model, RTModelWrapper}}, utils::data_handling::PeptideData};
use redeem_properties::utils::{peptdeep_utils::load_modifications, utils::get_device};
// use easypqp_core::property_prediction::DLFeatureScores;


struct PropertyPredictionScores<'a> {
    parameters: &'a InsilicoPQP, 
    peptides: &'a Vec<Peptide>, 
    fine_tune_data: Option<Vec<PeptideData>>
} 

impl<'a> PropertyPredictionScores<'a> {
    pub fn new(parameters: &'a InsilicoPQP, peptides: &'a Vec<Peptide>) -> Self {
        let fine_tune_data = if parameters.dl_feature_generators.fine_tune_config.fine_tune {
            Some(read_peptide_data_from_tsv(
                parameters.dl_feature_generators.fine_tune_config.train_data_path.clone(),
                parameters.dl_feature_generators.nce as i32,
                &parameters.dl_feature_generators.instrument,
            ).unwrap())
        } else {
            None
        };

        Self { parameters, peptides, fine_tune_data }
    }

    fn initialize_dl_models(&self) -> Result<DLModels> {
        let mut dl_models = DLModels::new();

        if self.parameters.dl_feature_generators.retention_time.is_not_empty() {
            dl_models.rt_model = Some(self.load_rt_model()?);
        }

        if self.parameters.dl_feature_generators.ion_mobility.is_not_empty() {
            dl_models.ccs_model = Some(self.load_ccs_model()?);
        }

        if self.parameters.dl_feature_generators.ms2_intensity.is_not_empty() {
            dl_models.ms2_model = Some(self.load_ms2_model()?);
        }

        if self.parameters.dl_feature_generators.is_not_empty() {
            dl_models.params = Some(self.create_dl_params());
        }

        Ok(dl_models)
    }

    fn create_dl_params(&self) -> redeem_properties::models::model_interface::Parameters {
        redeem_properties::models::model_interface::Parameters::new(
            &self.parameters.dl_feature_generators.instrument,
            self.parameters.dl_feature_generators.nce,
        )
    }

    fn remove_extension(&self, model_path: &str) -> String {
        let path = Path::new(model_path);
        match path.file_stem() {
            Some(stem) => stem.to_string_lossy().to_string(),
            None => model_path.to_string(),
        }
    }

    fn load_rt_model(&self) -> Result<RTModelWrapper> {
        let mut model = load_retention_time_model(
            &self.parameters.dl_feature_generators.retention_time.model_path,
            &self.parameters.dl_feature_generators.retention_time.constants_path,
            &self.parameters.dl_feature_generators.retention_time.architecture,
            get_device(&self.parameters.dl_feature_generators.device)?,
        )?;

        // Fine-tune the models if requested
        if self.parameters.dl_feature_generators.fine_tune_config.fine_tune {
            let n_fine_tune_data = self.fine_tune_data.as_ref().unwrap().len();
            if n_fine_tune_data < 10 {
                warn!("Insufficient training data: only {} samples. Skipping fine-tuning.", n_fine_tune_data);
            } else {

                // Load modifications map
                let modifications = load_modifications()?;

                // Fine-tune the model
                model.fine_tune(&self.fine_tune_data.as_ref().unwrap(), modifications, self.parameters.dl_feature_generators.fine_tune_config.batch_size, self.parameters.dl_feature_generators.fine_tune_config.learning_rate, self.parameters.dl_feature_generators.fine_tune_config.epochs)?;

                if self.parameters.dl_feature_generators.fine_tune_config.save_model {
                    let model_path = self.parameters.dl_feature_generators.retention_time.model_path.clone();
                    let model_path = self.remove_extension(&model_path);
                    let model_path = format!("{}_fine_tuned.safetensors", model_path);
                    info!("Saving fine-tuned RT model to {}", model_path);
                    model.save(&model_path)?;
                }
            }
        }

        model.set_evaluation_mode();
        Ok(model)
    }

    fn load_ccs_model(&self) -> Result<CCSModelWrapper> {
        let mut model = load_collision_cross_section_model(
            &self.parameters.dl_feature_generators.ion_mobility.model_path,
            &self.parameters.dl_feature_generators.ion_mobility.constants_path,
            &self.parameters.dl_feature_generators.ion_mobility.architecture,
            get_device(&self.parameters.dl_feature_generators.device)?,
        )?;

        // Fine-tune the models if requested
        if self.parameters.dl_feature_generators.fine_tune_config.fine_tune {
            let n_fine_tune_data = self.fine_tune_data.as_ref().unwrap().len();
            if n_fine_tune_data < 10 {
                warn!("Insufficient training data: only {} samples. Skipping fine-tuning.", n_fine_tune_data);
            } else {
                // Load modifications map
                let modifications = load_modifications()?;

                // Fine-tune the model
                model.fine_tune(&self.fine_tune_data.as_ref().unwrap(), modifications, self.parameters.dl_feature_generators.fine_tune_config.batch_size, self.parameters.dl_feature_generators.fine_tune_config.learning_rate, self.parameters.dl_feature_generators.fine_tune_config.epochs)?;

                if self.parameters.dl_feature_generators.fine_tune_config.save_model {
                    let model_path = self.parameters.dl_feature_generators.ion_mobility.model_path.clone();
                    let model_path = self.remove_extension(&model_path);
                    let model_path = format!("{}_fine_tuned.safetensors", model_path);
                    info!("Saving fine-tuned CCS model to {}", model_path);
                    model.save(&model_path)?;
                }
            }

        }

        Ok(model)
    }

    fn load_ms2_model(&self) -> Result<MS2ModelWrapper> {
        let mut model = load_ms2_model(
            &self.parameters.dl_feature_generators.ms2_intensity.model_path,
            &self.parameters.dl_feature_generators.ms2_intensity.constants_path,
            &self.parameters.dl_feature_generators.ms2_intensity.architecture,
            get_device(&self.parameters.dl_feature_generators.device)?,
        )?;

        // Fine-tune the models if requested
        if self.parameters.dl_feature_generators.fine_tune_config.fine_tune {
            let n_fine_tune_data = self.fine_tune_data.as_ref().unwrap().len();
            if n_fine_tune_data < 10 {
                warn!("Insufficient training data: only {} samples. Skipping fine-tuning.", n_fine_tune_data);
            } else {
                // Load modifications map
                let modifications = load_modifications()?;

                // Fine-tune the model
                model.fine_tune(&self.fine_tune_data.as_ref().unwrap(), modifications, self.parameters.dl_feature_generators.fine_tune_config.batch_size, self.parameters.dl_feature_generators.fine_tune_config.learning_rate, self.parameters.dl_feature_generators.fine_tune_config.epochs)?;

                if self.parameters.dl_feature_generators.fine_tune_config.save_model {
                    let model_path = self.parameters.dl_feature_generators.ms2_intensity.model_path.clone();
                    let model_path = self.remove_extension(&model_path);
                    let model_path = format!("{}_fine_tuned.safetensors", model_path);
                    info!("Saving fine-tuned MS2 Intensity model to {}", model_path);
                    model.save(&model_path)?;
                }
            }
        }

        Ok(model)
    }

    pub fn predict_properties(&mut self) -> Result<Vec<PeptideProperties>> {
        let dl_models = self.initialize_dl_models()?;
        
        // Load modifications map
        let modifications = load_modifications()?;

        // Create a PropertyPrediction instance and predict properties
        let mut property_prediction = PropertyPrediction::new(&self.peptides, &self.parameters.insilico_settings, dl_models, modifications, self.parameters.dl_feature_generators.batch_size);

        
        let assays = property_prediction.predict_properties()?;

        Ok(assays)
    }
}


pub struct Runner {
    peptides: Vec<Peptide>,
    parameters: InsilicoPQP
}

impl Runner {
    pub fn new(parameters: InsilicoPQP) -> anyhow::Result<Self> {
        let start = Instant::now();
        let fasta = easypqp_core::util::read_fasta(
            &parameters.database.fasta,
            &parameters.database.decoy_tag,
            parameters.database.generate_decoys,
        )
        .with_context(|| {
            format!(
                "Failed to build database from `{}`",
                parameters.database.fasta
            )
        })?;

        // let database = parameters.database.clone().build(fasta);
        let peptides = parameters.database.clone().digest(&fasta);

        info!(
            "generated {} peptides in {}ms",
            peptides.len(),
            (Instant::now() - start).as_millis()
        );

        Ok(Self {
            peptides,
            parameters
        })
    }

    pub fn run(mut self) -> anyhow::Result<()> {

        let start_time = Instant::now();
        let mut property_prediction_scores = PropertyPredictionScores::new(&self.parameters,  &self.peptides);
        let assays: Vec<PeptideProperties> = property_prediction_scores.predict_properties()?;

        write_assays_to_tsv(&assays, &self.peptides, &self.parameters.output_file, &self.parameters.insilico_settings)?;

        let execution_time = Instant::now() - start_time;
        log::info!("Insilico library generation: {:8} min", execution_time.as_secs() / 60);

        let path = "easypqp_insilico.json";
        println!("{}", serde_json::to_string_pretty(&self.parameters.as_serializable())?);
        let bytes = serde_json::to_vec_pretty(&self.parameters.as_serializable())?;
        write_bytes_to_file(path, &bytes)?;

        Ok(())
    }
}

use std::fs::File;
use std::io::Write;

fn write_bytes_to_file(path: &str, bytes: &[u8]) -> std::io::Result<()> {
    let path = Path::new(path);
    let mut file = File::create(path)?;
    file.write_all(bytes)?;
    Ok(())
}