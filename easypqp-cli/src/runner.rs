use anyhow::{Context, Result};
use log::{info, warn};
use sage_core::peptide::Peptide;
use std::{path::Path, time::Instant};
use redeem_properties::utils::{peptdeep_utils::load_modifications, utils::get_device};
use redeem_properties::{
    models::{
        ccs_model::{load_collision_cross_section_model, CCSModelWrapper},
        model_interface::DLModels,
        ms2_model::{load_ms2_model, MS2ModelWrapper},
        rt_model::{load_retention_time_model, RTModelWrapper},
    },
    utils::data_handling::{PeptideData, TargetNormalization},
};
use easypqp_core::{
    property_prediction::PropertyPrediction, tuning_data::read_peptide_data_from_tsv,
    util::write_bytes_to_file, PeptideProperties,
};

use crate::input::InsilicoPQP;
use crate::output::write_assays_to_tsv;
use easypqp_core::unimod::UnimodDb;

struct PropertyPredictionScores<'a> {
    parameters: &'a InsilicoPQP,
    peptides: &'a [Peptide],
    fine_tune_data: Option<Vec<PeptideData>>,
}

impl<'a> PropertyPredictionScores<'a> {
    pub fn new(parameters: &'a InsilicoPQP, peptides: &'a [Peptide]) -> Self {
        let fine_tune_data = if parameters.dl_feature_generators.fine_tune_config.fine_tune {
            Some(
                read_peptide_data_from_tsv(
                    parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .train_data_path
                        .clone(),
                    parameters.dl_feature_generators.nce as i32,
                    &parameters.dl_feature_generators.instrument,
                )
                .unwrap(),
            )
        } else {
            None
        };

        Self {
            parameters,
            peptides,
            fine_tune_data,
        }
    }

    fn initialize_dl_models(&self) -> Result<DLModels> {
        let mut dl_models = DLModels::new();

        if self
            .parameters
            .dl_feature_generators
            .retention_time
            .is_not_empty()
        {
            dl_models.rt_model = Some(self.load_rt_model()?);
        }

        if self
            .parameters
            .dl_feature_generators
            .ion_mobility
            .is_not_empty()
        {
            dl_models.ccs_model = Some(self.load_ccs_model()?);
        }

        if self
            .parameters
            .dl_feature_generators
            .ms2_intensity
            .is_not_empty()
        {
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
        let model_path = self
            .parameters
            .dl_feature_generators
            .retention_time
            .model_path
            .clone();
        let const_path_str = self
            .parameters
            .dl_feature_generators
            .retention_time
            .constants_path
            .clone();
        let constants_opt = if const_path_str.is_empty() {
            None
        } else {
            Some(const_path_str)
        };

        let mut model = load_retention_time_model(
            &model_path,
            constants_opt.as_ref(),
            &self
                .parameters
                .dl_feature_generators
                .retention_time
                .architecture,
            get_device(&self.parameters.dl_feature_generators.device)?,
        )?;

        // Fine-tune the models if requested
        if self
            .parameters
            .dl_feature_generators
            .fine_tune_config
            .fine_tune
        {
            let n_fine_tune_data = self.fine_tune_data.as_ref().unwrap().len();
            if n_fine_tune_data < 10 {
                warn!(
                    "Insufficient training data: only {} samples. Skipping fine-tuning.",
                    n_fine_tune_data
                );
            } else {
                // Load modifications map
                let modifications = load_modifications()?;

                // Fine-tune the model
                model.fine_tune(
                    &self.fine_tune_data.as_ref().unwrap(),
                    modifications,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .batch_size,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .learning_rate,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .epochs,
                    TargetNormalization::None,
                )?;

                if self
                    .parameters
                    .dl_feature_generators
                    .fine_tune_config
                    .save_model
                {
                    let model_path = self
                        .parameters
                        .dl_feature_generators
                        .retention_time
                        .model_path
                        .clone();
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
        let model_path = self
            .parameters
            .dl_feature_generators
            .ion_mobility
            .model_path
            .clone();
        let const_path_str = self
            .parameters
            .dl_feature_generators
            .ion_mobility
            .constants_path
            .clone();

        let mut model = load_collision_cross_section_model(
            &model_path,
            &const_path_str,
            &self
                .parameters
                .dl_feature_generators
                .ion_mobility
                .architecture,
            get_device(&self.parameters.dl_feature_generators.device)?,
        )?;

        // Fine-tune the models if requested
        if self
            .parameters
            .dl_feature_generators
            .fine_tune_config
            .fine_tune
        {
            let n_fine_tune_data = self.fine_tune_data.as_ref().unwrap().len();
            if n_fine_tune_data < 10 {
                warn!(
                    "Insufficient training data: only {} samples. Skipping fine-tuning.",
                    n_fine_tune_data
                );
            } else {
                // Load modifications map
                let modifications = load_modifications()?;

                // Fine-tune the model
                model.fine_tune(
                    &self.fine_tune_data.as_ref().unwrap(),
                    modifications,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .batch_size,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .learning_rate,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .epochs,
                    TargetNormalization::None,
                )?;

                if self
                    .parameters
                    .dl_feature_generators
                    .fine_tune_config
                    .save_model
                {
                    let model_path = self
                        .parameters
                        .dl_feature_generators
                        .ion_mobility
                        .model_path
                        .clone();
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
        let model_path = self
            .parameters
            .dl_feature_generators
            .ms2_intensity
            .model_path
            .clone();
        let const_path_str = self
            .parameters
            .dl_feature_generators
            .ms2_intensity
            .constants_path
            .clone();

        let mut model = load_ms2_model(
            &model_path,
            &const_path_str,
            &self
                .parameters
                .dl_feature_generators
                .ms2_intensity
                .architecture,
            get_device(&self.parameters.dl_feature_generators.device)?,
        )?;

        // Fine-tune the models if requested
        if self
            .parameters
            .dl_feature_generators
            .fine_tune_config
            .fine_tune
        {
            let n_fine_tune_data = self.fine_tune_data.as_ref().unwrap().len();
            if n_fine_tune_data < 10 {
                warn!(
                    "Insufficient training data: only {} samples. Skipping fine-tuning.",
                    n_fine_tune_data
                );
            } else {
                // Load modifications map
                let modifications = load_modifications()?;

                // Fine-tune the model
                model.fine_tune(
                    &self.fine_tune_data.as_ref().unwrap(),
                    modifications,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .batch_size,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .learning_rate,
                    self.parameters
                        .dl_feature_generators
                        .fine_tune_config
                        .epochs,
                    TargetNormalization::None,
                )?;

                if self
                    .parameters
                    .dl_feature_generators
                    .fine_tune_config
                    .save_model
                {
                    let model_path = self
                        .parameters
                        .dl_feature_generators
                        .ms2_intensity
                        .model_path
                        .clone();
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
        let mut property_prediction = PropertyPrediction::new(
            &self.peptides,
            &self.parameters.insilico_settings,
            dl_models,
            modifications,
            self.parameters.dl_feature_generators.batch_size,
        );

        let assays = property_prediction.predict_properties()?;

        Ok(assays)
    }
}

pub struct Runner {
    peptides: Vec<Peptide>,
    parameters: InsilicoPQP,
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
            parameters,
        })
    }

    pub fn run(self) -> anyhow::Result<()> {
        let max_chunk_size = self.parameters.peptide_chunking.resolve(self.parameters.dl_feature_generators.batch_size);
    
        // Only write log if max_chunk_size is less than the number of peptides
        if max_chunk_size < self.peptides.len() {
            log::info!("Processing max {} peptides per chunk", max_chunk_size);
        }

        if Path::new(&self.parameters.output_file).exists() {
            log::warn!(
                "There is an existing file: {}. It will be overwritten.",
                self.parameters.output_file
            );
            std::fs::remove_file(&self.parameters.output_file)?;
        }
    
        let start_time = Instant::now();
        
        // Determine output paths based on parquet flag
        let (output_path, parquet_path) = if self.parameters.parquet_output {
            // When parquet is enabled, derive .parquet filename from output_file
            let base = self.parameters.output_file.trim_end_matches(".tsv");
            let parquet_file = format!("{}.parquet", base);
            (None, Some(parquet_file))
        } else {
            // TSV output
            (Some(self.parameters.output_file.clone()), None)
        };
        
        // Prepare optional streaming parquet writer (created only when feature enabled).
        #[cfg(feature = "parquet")]
        let mut parquet_writer = if let Some(ref p) = parquet_path {
            // instantiate parquet writer
            Some(crate::output::ParquetChunkWriter::try_new(p, &self.parameters.insilico_settings, &self.parameters.database.decoy_tag)?)
        } else {
            None
        };

        // Build UniMod database once for TSV output reannotation
        let unimod_db = if self.parameters.insilico_settings.unimod_annotation {
            let db_result = if let Some(ref xml_path) = self.parameters.insilico_settings.unimod_xml_path {
                info!("Loading custom UniMod XML from '{}'", xml_path);
                UnimodDb::from_file(xml_path, self.parameters.insilico_settings.max_delta_unimod)
            } else {
                UnimodDb::from_embedded(self.parameters.insilico_settings.max_delta_unimod)
            };
            match db_result {
                Ok(db) => {
                    info!("Loaded UniMod database for modification reannotation (max_delta={} Da)", self.parameters.insilico_settings.max_delta_unimod);
                    Some(db)
                }
                Err(e) => {
                    warn!("Failed to load UniMod database, skipping reannotation: {}", e);
                    None
                }
            }
        } else {
            None
        };

        // No in-memory report accumulation: when Parquet output + report is requested
        // we'll generate the report directly from the Parquet file after closing the writer.

        let mut tsv_peptide_offset: usize = 0;

        for (i, peptide_chunk) in self.peptides.chunks(max_chunk_size).enumerate() {
            if max_chunk_size < self.peptides.len() {
                log::info!(
                    "Processing chunk {} of {} with {} peptides",
                    i + 1,
                    (self.peptides.len() + max_chunk_size - 1) / max_chunk_size,
                    peptide_chunk.len()
                );
            }

            let mut predictor = PropertyPredictionScores::new(&self.parameters, peptide_chunk);
            let assays = predictor.predict_properties()?;

            if parquet_path.is_some() {
                #[cfg(feature = "parquet")]
                {
                    if let Some(ref mut w) = parquet_writer {
                        w.write_chunk(&assays, peptide_chunk, &self.parameters.insilico_settings)?;
                    }
                }

                #[cfg(not(feature = "parquet"))]
                {
                    warn!("Parquet support not enabled in this build. Rebuild with feature 'parquet' to use --parquet.");
                }

                // No-op for report accumulation here when using Parquet; report will be
                // generated from the written Parquet file after the writer is closed.
            } else {
                write_assays_to_tsv(
                    &assays,
                    peptide_chunk, // use the chunk peptides, not full self.peptides
                    output_path.as_ref().unwrap(),
                    &self.parameters.insilico_settings,
                    unimod_db.as_ref(),
                    &self.parameters.database.decoy_tag,
                    tsv_peptide_offset,
                )?;
            }

            // Advance offset for globally unique IDs across chunks
            tsv_peptide_offset += peptide_chunk.len();
        }

        // Close parquet writer if present
        #[cfg(feature = "parquet")]
        if let Some(w) = parquet_writer.take() {
            w.close()?;
        }

        // Generate HTML report: if Parquet output was requested, generate report from
        // the Parquet file; otherwise use the TSV-based generator on the output TSV.
        if self.parameters.write_report {
            if let Some(ref p) = parquet_path {
                let report_path = format!("{}.html", p.trim_end_matches(".parquet"));
                match crate::output::generate_html_report_from_parquet(p, &report_path) {
                    Ok(_) => info!("Generated HTML report: {}", report_path),
                    Err(e) => warn!("Failed to generate HTML report from Parquet: {}", e),
                }
            } else {
                let report_path = format!("{}.html", self.parameters.output_file.trim_end_matches(".tsv"));
                match crate::output::generate_html_report(&self.parameters.output_file, &report_path) {
                    Ok(_) => info!("Generated HTML report: {}", report_path),
                    Err(e) => warn!("Failed to generate HTML report: {}", e),
                }
            }
        }
    
        let execution_time = Instant::now() - start_time;
        log::info!(
            "Insilico library generation: {:8} min",
            execution_time.as_secs() / 60
        );
    
        let path = "easypqp_insilico.json";
        let json = serde_json::to_string_pretty(&self.parameters.as_serializable())?;
        println!("{}", json);
        let bytes = serde_json::to_vec_pretty(&self.parameters.as_serializable())?;
        write_bytes_to_file(path, &bytes)?;
    
        Ok(())
    }
    
}
