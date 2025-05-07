use std::collections::HashMap;

use anyhow::{ensure, Context, Result};
use clap::ArgMatches;
use serde::{Deserialize, Serialize};
use sage_core::{database::{Builder, EnzymeBuilder, Parameters}, modification::ModificationSpecificity};
use easypqp_core::InsilicoPQPSettings;


#[derive(Serialize, Clone)]
/// Actual insilico parameters - may include overrides or default values not set by user
pub struct InsilicoPQP {
    pub version: String,
    pub database: Parameters,
    pub insilico_settings: InsilicoPQPSettings,
    pub dl_feature_generators: DLFeatureGeneratorSettings,
    pub output_file: String,
}

#[derive(Serialize, Deserialize, Default, Debug, Clone)]
pub struct DLModel {
    pub model_path: String,
    pub constants_path: String,
    pub architecture: String,
}

impl DLModel {
    pub fn is_not_empty(&self) -> bool {
        !self.model_path.is_empty()
            && !self.constants_path.is_empty()
            && !self.architecture.is_empty()
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FineTuneOptions {
    pub fine_tune: Option<bool>,
    pub train_data_path: Option<String>,
    pub batch_size: Option<usize>,
    pub epochs: Option<usize>,
    pub learning_rate: Option<f64>,
    pub save_model: Option<bool>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct FineTuneSettings {
    pub fine_tune: bool,
    pub train_data_path: String,
    #[serde(default = "FineTuneSettings::default_batch_size")]
    pub batch_size: usize,
    #[serde(default = "FineTuneSettings::default_epochs")]
    pub epochs: usize,
    #[serde(default = "FineTuneSettings::default_learning_rate")]
    pub learning_rate: f64,
    #[serde(default = "FineTuneSettings::default_save_model")]
    pub save_model: bool,
}

impl FineTuneSettings {
    fn default_batch_size() -> usize {
        256
    }
    fn default_epochs() -> usize {
        3
    }
    fn default_learning_rate() -> f64 {
        0.001
    }
    fn default_save_model() -> bool {
        false
    }

    pub fn default() -> Self {
        Self {
            fine_tune: false,
            train_data_path: String::new(),
            batch_size: Self::default_batch_size(),
            epochs: Self::default_epochs(),
            learning_rate: Self::default_learning_rate(),
            save_model: Self::default_save_model(),
        }
    }
}

impl Default for FineTuneSettings {
    fn default() -> Self {
        Self {
            fine_tune: false,
            train_data_path: String::new(),
            batch_size: 256,
            epochs: 3,
            learning_rate: 0.001,
            save_model: false,
        }
    }
}

impl From<FineTuneOptions> for FineTuneSettings {
    fn from(value: FineTuneOptions) -> Self {
        Self {
            fine_tune: value.fine_tune.unwrap_or(false),
            train_data_path: value.train_data_path.unwrap_or(String::new()),
            batch_size: value.batch_size.unwrap_or(256),
            epochs: value.epochs.unwrap_or(3),
            learning_rate: value.learning_rate.unwrap_or(0.001),
            save_model: value.save_model.unwrap_or(false),
        }
    }
}

impl FineTuneSettings {
    pub fn is_not_empty(&self) -> bool {
        self.fine_tune
    }
}

#[derive(Serialize, Deserialize, Default, Debug, Clone)]
pub struct DLFeatureGenerators {
    pub predict_properties: Option<bool>,
    pub retention_time: Option<DLModel>,
    pub ion_mobility: Option<DLModel>,
    pub ms2_intensity: Option<DLModel>,
    pub device: Option<String>,
    pub fine_tune_config: Option<FineTuneSettings>,
    pub instrument: Option<String>,
    pub nce: Option<f32>,
    pub batch_size: Option<usize>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct DLFeatureGeneratorSettings {
    pub predict_properties: bool,
    pub retention_time: DLModel,
    pub ion_mobility: DLModel,
    pub ms2_intensity: DLModel,
    pub device: String,
    pub fine_tune_config: FineTuneSettings,
    pub instrument: String,
    pub nce: f32,
    pub batch_size: usize,
}

impl Default for DLFeatureGeneratorSettings {
    fn default() -> Self {
        Self {
            predict_properties: false,
            retention_time: DLModel {
                model_path: String::new(),
                constants_path: String::new(),
                architecture: String::new(),
            },
            ion_mobility: DLModel {
                model_path: String::new(),
                constants_path: String::new(),
                architecture: String::new(),
            },
            ms2_intensity: DLModel {
                model_path: String::new(),
                constants_path: String::new(),
                architecture: String::new(),
            },
            device: "cpu".into(),
            fine_tune_config: FineTuneSettings::default(),
            instrument: "timsTOF".into(),
            nce: 20.0,
            batch_size: 64,
        }
    }
}

impl DLFeatureGeneratorSettings {
    pub(crate) fn is_not_empty(&self) -> bool {
        self.retention_time.is_not_empty()
            || self.ion_mobility.is_not_empty()
            || self.ms2_intensity.is_not_empty()
    }
}

impl From<DLFeatureGenerators> for DLFeatureGeneratorSettings {
    fn from(value: DLFeatureGenerators) -> Self {
        let mut retention_time_model_config = DLModel {
            model_path: String::new(),
            constants_path: String::new(),
            architecture: String::new(),
        };
        let mut ion_mobility_model_config = DLModel {
            model_path: String::new(),
            constants_path: String::new(),
            architecture: String::new(),
        };
        let mut ms2_intensity_model_config = DLModel {
            model_path: String::new(),
            constants_path: String::new(),
            architecture: String::new(),
        };

        if value.predict_properties.unwrap_or(false)
            && value.retention_time.clone().is_none()
            && value.ion_mobility.clone().is_none()
            && value.ms2_intensity.clone().is_none()
        {
            log::info!("Predicting properties is enabled, but no model configurations are provided. Will attempt to retrieve and use AlphaPeptDeep generic pretrained models.");
            let _ = redeem_properties::utils::peptdeep_utils::download_pretrained_models_exist();
            retention_time_model_config = DLModel {
                model_path: "data/peptdeep_generic_pretrained_models/generic/rt.pth".to_string(),
                constants_path:
                    "data/peptdeep_generic_pretrained_models/generic/rt.pth.model_const.yaml"
                        .to_string(),
                architecture: "rt_cnn_lstm".to_string(),
            };

            // Note: Probably only want to use and set CCS model if instrument is TIMSTOF, if the model config is not explicitly set
            if value
                .instrument
                .clone()
                .unwrap_or("QE".into())
                .to_uppercase()
                == "TIMSTOF".to_string()
            {
                ion_mobility_model_config = DLModel {
                    model_path: "data/peptdeep_generic_pretrained_models/generic/ccs.pth"
                        .to_string(),
                    constants_path:
                        "data/peptdeep_generic_pretrained_models/generic/ccs.pth.model_const.yaml"
                            .to_string(),
                    architecture: "ccs_cnn_lstm".to_string(),
                };
            }

            ms2_intensity_model_config = DLModel {
                model_path: "data/peptdeep_generic_pretrained_models/generic/ms2.pth".to_string(),
                constants_path:
                    "data/peptdeep_generic_pretrained_models/generic/ms2.pth.model_const.yaml"
                        .to_string(),
                architecture: "ms2_bert".to_string(),
            };
        }

        Self {
            predict_properties: value.predict_properties.unwrap_or(false),
            retention_time: value
                .retention_time
                .unwrap_or_else(|| retention_time_model_config),
            ion_mobility: value
                .ion_mobility
                .unwrap_or_else(|| ion_mobility_model_config),
            ms2_intensity: value
                .ms2_intensity
                .unwrap_or_else(|| ms2_intensity_model_config),
            device: value.device.unwrap_or("cpu".into()),
            fine_tune_config: FineTuneSettings::from(value.fine_tune_config.unwrap_or_default()),
            instrument: value.instrument.unwrap_or("QE".into()),
            nce: value.nce.unwrap_or(20.0),
            batch_size: value.batch_size.unwrap_or(64),
        }
    }
}

#[derive(Deserialize)]
pub struct Input {
    pub database: Builder,
    pub insilico_settings: InsilicoPQPSettings,
    pub dl_feature_generators: Option<DLFeatureGenerators>,
    pub output_file: Option<String>,
}

impl Input {
    /// Load parameters from a JSON file and validate them.
    pub fn from_arguments(matches: ArgMatches) -> Result<Self> {
        let path = matches
            .get_one::<String>("parameters")
            .expect("required parameters");

        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

        // Handle JSON configuration overrides
        if let Some(fasta) = matches.get_one::<String>("fasta") {
            input.database.fasta = Some(fasta.into());
        }
        // Handle JSON configuration overrides
        if let Some(output_file) = matches.get_one::<String>("output_file") {
            input.output_file = Some(output_file.into());
        }

        // avoid to later panic if these parameters are not set (but doesn't check if files exist)

        ensure!(
            input.database.fasta.is_some(),
            "`database.fasta` must be set. For more information try '--help'"
        );

        Ok(input)
    }

    /// Load parameters from a JSON file and validate them. This function is used for the python bindings.
    pub fn from_passed_arguments(
        path: &str,
        fasta: Option<String>,
        output_file: Option<String>,
        generate_decoys: Option<bool>,
        decoy_tag: Option<String>,
        precursor_charge: Option<Vec<u8>>,
        max_fragment_charge: Option<usize>,
        min_transitions: Option<u8>,
        max_transitions: Option<u8>,
        fragmentation_model: Option<String>,
        allowed_fragment_types: Option<Vec<String>>,
        fine_tune: Option<bool>,
        train_data_path: Option<String>,
        batch_size: Option<usize>,
    ) -> Result<Self> {
        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

        // Handle JSON configuration overrides
        if let Some(fasta) = fasta {
            input.database.fasta = Some(fasta.into());
        }
        if let Some(output_file) = output_file {
            input.output_file = Some(output_file.into());
        }
        if let Some(generate_decoys) = generate_decoys {
            input.database.generate_decoys = Some(generate_decoys);
        }
        if let Some(decoy_tag) = decoy_tag {
            input.database.decoy_tag = Some(decoy_tag);
        }
        if let Some(precursor_charge) = precursor_charge {
            input.insilico_settings.precursor_charge = precursor_charge;
        }
        if let Some(max_fragment_charge) = max_fragment_charge {
            input.insilico_settings.max_fragment_charge = max_fragment_charge;
        }
        if let Some(min_transitions) = min_transitions {
            input.insilico_settings.min_transitions = min_transitions;
        }
        if let Some(max_transitions) = max_transitions {
            input.insilico_settings.max_transitions = max_transitions;
        }
        if let Some(fragmentation_model) = fragmentation_model {
            input.insilico_settings.fragmentation_model = fragmentation_model;
        }
        if let Some(allowed_fragment_types) = allowed_fragment_types {
            input.insilico_settings.allowed_fragment_types = allowed_fragment_types;
        }
        if let Some(fine_tune) = fine_tune {
            input.dl_feature_generators.clone().unwrap().fine_tune_config.unwrap().fine_tune = fine_tune;
        }
        if let Some(train_data_path) = train_data_path {
            input.dl_feature_generators.clone().unwrap().fine_tune_config.unwrap().train_data_path = train_data_path;
        }
        if let Some(batch_size) = batch_size {
            input.dl_feature_generators.clone().unwrap().fine_tune_config.unwrap().batch_size = batch_size;
        }

        // avoid to later panic if these parameters are not set (but doesn't check if files exist)
        ensure!(
            input.database.fasta.is_some(),
            "`fasta` must be provided. For more information try '--help'"
        );

        Ok(input)
    }

    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        easypqp_core::util::read_json(path).map_err(anyhow::Error::from)
    }

    pub fn build(mut self) -> anyhow::Result<InsilicoPQP> {
        let database = self.database.make_parameters();

        Ok(InsilicoPQP {
            version: clap::crate_version!().into(),
            database,
            insilico_settings: self.insilico_settings.clone(),
            dl_feature_generators: self
                .dl_feature_generators
                .map(Into::into)
                .unwrap_or_default(),
            output_file: self
                .output_file
                .clone()
                .unwrap_or_else(|| "insilico_library.tsv".into()),
        })
    }
}


/// The below is for data transfer to avoid writing out fields specific for sage

// 1. Create a DTO for Sage's Parameters that skips unwanted fields
#[derive(Serialize)]
pub struct SageParametersDto<'a> {
    // Fields you want to keep
    pub enzyme: &'a EnzymeBuilder,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub static_mods: &'a HashMap<ModificationSpecificity, f32>,
    pub variable_mods: &'a HashMap<ModificationSpecificity, Vec<f32>>,
    pub max_variable_mods: usize,
    pub decoy_tag: &'a String,
    pub generate_decoys: bool,
    pub fasta: &'a String,
}

// 2. Create a DTO for your InsilicoPQP struct
#[derive(Serialize)]
pub struct InsilicoPQPDto<'a> {
    pub version: &'a String,
    pub database: SageParametersDto<'a>,
    pub insilico_settings: &'a InsilicoPQPSettings,
    pub dl_feature_generators: &'a DLFeatureGeneratorSettings,
    pub output_file: &'a String,
}

// 3. Implement conversion methods
impl InsilicoPQP {
    pub fn as_serializable(&self) -> InsilicoPQPDto {
        InsilicoPQPDto {
            version: &self.version,
            database: SageParametersDto {
                enzyme: &self.database.enzyme,
                peptide_min_mass: self.database.peptide_min_mass,
                peptide_max_mass: self.database.peptide_max_mass,
                static_mods: &self.database.static_mods,
                variable_mods: &self.database.variable_mods,
                max_variable_mods: self.database.max_variable_mods,
                decoy_tag: &self.database.decoy_tag,
                generate_decoys: self.database.generate_decoys,
                fasta: &self.database.fasta
            },
            insilico_settings: &self.insilico_settings,
            dl_feature_generators: &self.dl_feature_generators,
            output_file: &self.output_file,
        }
    }
}