use std::{collections::HashMap, path::Path};

use anyhow::{ensure, Context, Result};
use clap::ArgMatches;
use serde::{Deserialize, Serialize};
use schemars::{JsonSchema, schema_for};
use sage_core::{database::{Builder, EnzymeBuilder, Parameters}, modification::ModificationSpecificity};
use easypqp_core::{util::auto_chunk_size, InsilicoPQPSettings};


#[derive(Serialize, Deserialize, Clone, Copy, Debug, JsonSchema)]
#[schemars(description = "Peptide chunking strategy for memory management")]
pub struct ChunkingStrategy(pub usize);


impl Default for ChunkingStrategy {
    fn default() -> Self {
        Self(0)
    }
}

impl ChunkingStrategy {
    pub fn resolve(self, batch_size: usize) -> usize {
        if self.0 == 0 {
            auto_chunk_size(1024 * batch_size, 0.5)
        } else {
            self.0
        }
    }
}


#[derive(Serialize, Clone)]
/// Actual insilico parameters - may include overrides or default values not set by user
pub struct InsilicoPQP {
    pub version: String,
    pub database: Parameters,
    pub insilico_settings: InsilicoPQPSettings,
    pub dl_feature_generators: DLFeatureGeneratorSettings,
    pub peptide_chunking: ChunkingStrategy,
    pub output_file: String,
    pub write_report: bool,
    pub parquet_output: bool,
}

#[derive(Serialize, Deserialize, Default, Debug, Clone, JsonSchema)]
#[schemars(description = "Deep learning model file paths")]
pub struct DLModel {
    /// Path to model weights file (.pth or .safetensors)
    pub model_path: String,
    
    /// Path to model constants YAML file
    pub constants_path: String,
    
    /// Path to model architecture Python file
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

#[derive(Serialize, Deserialize, Debug, Clone, JsonSchema)]
#[schemars(description = "Fine-tuning configuration for transfer learning")]
pub struct FineTuneSettings {
    /// Enable fine-tuning
    pub fine_tune: bool,
    
    /// Path to training data TSV file. Required columns: 'sequence' (modified sequence annotated with square bracket mass shift, e.g., MGC[+57.0215]AAR), 'precursor_charge', 'retention_time', 'ion_mobility' (if using timsTOF),
    /// 'fragment_type', 'fragment_series_number', 'product_charge', 'intensity'
    #[schemars(description = "TSV file with columns: sequence, precursor_charge, retention_time, ion_mobility (if applicable), fragment_type, fragment_series_number, product_charge, intensity")]
    pub train_data_path: String,
    
    /// Batch size for training
    #[serde(default = "FineTuneSettings::default_batch_size")]
    pub batch_size: usize,
    
    /// Number of training epochs
    #[serde(default = "FineTuneSettings::default_epochs")]
    pub epochs: usize,
    
    /// Learning rate
    #[serde(default = "FineTuneSettings::default_learning_rate")]
    pub learning_rate: f64,
    
    /// Save fine-tuned model to disk
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

#[derive(Serialize, Deserialize, Default, Debug, Clone, JsonSchema)]
#[schemars(description = "Deep learning model configuration for property prediction")]
pub struct DLFeatureGenerators {
    /// Retention time prediction model
    #[schemars(description = "Custom RT model (model_path, constants_path, architecture)")]
    pub retention_time: Option<DLModel>,
    
    /// Ion mobility prediction model
    #[schemars(description = "Custom IM/CCS model (model_path, constants_path, architecture)")]
    pub ion_mobility: Option<DLModel>,
    
    /// MS2 intensity prediction model
    #[schemars(description = "Custom MS2 model (model_path, constants_path, architecture)")]
    pub ms2_intensity: Option<DLModel>,
    
    /// Compute device
    #[schemars(description = "Device for inference: 'cpu', 'cuda', or 'mps' (default: 'cpu')")]
    pub device: Option<String>,
    
    /// Fine-tuning configuration
    #[schemars(description = "Optional fine-tuning settings for transfer learning")]
    pub fine_tune_config: Option<FineTuneSettings>,
    
    /// Instrument type
    #[schemars(description = "Instrument type: 'QE' or 'timsTOF' (default: 'timsTOF')")]
    pub instrument: Option<String>,
    
    /// Normalized collision energy
    #[schemars(description = "NCE value for fragmentation (default: 20.0)")]
    pub nce: Option<f32>,
    
    /// Inference batch size
    #[schemars(description = "Batch size for model inference (default: 64)")]
    pub batch_size: Option<usize>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct DLFeatureGeneratorSettings {
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

        // If no models are explicitly configured, download and use pretrained models
        if value.retention_time.clone().is_none()
            && value.ion_mobility.clone().is_none()
            && value.ms2_intensity.clone().is_none()
        {
            log::info!("No model configurations provided. Will attempt to retrieve and use pretrained models.");
            let _ = redeem_properties::utils::peptdeep_utils::download_pretrained_models_exist();
            retention_time_model_config = DLModel {
                model_path: "data/pretrained_models/pretrained_models/redeem/20251205_100_epochs_min_max_rt_cnn_tf.safetensors".to_string(),
                constants_path:
                    "data/pretrained_models/pretrained_models/alphapeptdeep/generic/rt.pth.model_const.yaml"
                        .to_string(),
                architecture: "rt_cnn_tf".to_string(),
            };
            log::info!("Pre-trained retention time model (architecture: {}): {}", retention_time_model_config.architecture,retention_time_model_config.model_path);

            // Note: Probably only want to use and set CCS model if instrument is TIMSTOF, if the model config is not explicitly set
            if value
                .instrument
                .clone()
                .unwrap_or("QE".into())
                .to_uppercase()
                == "TIMSTOF".to_string()
            {
                ion_mobility_model_config = DLModel {
                    model_path: "data/pretrained_models/pretrained_models/redeem/20251205_500_epochs_early_stopped_100_min_max_ccs_cnn_tf.safetensors"
                        .to_string(),
                    constants_path:
                        "data/pretrained_models/pretrained_models/alphapeptdeep/generic/ccs.pth.model_const.yaml"
                            .to_string(),
                    architecture: "ccs_cnn_tf".to_string(),
                };
                log::info!("Pre-trained ion mobility model (architecture: {}): {}", ion_mobility_model_config.architecture,ion_mobility_model_config.model_path);
            }

            ms2_intensity_model_config = DLModel {
                model_path: "data/pretrained_models/pretrained_models/alphapeptdeep/generic/ms2.pth".to_string(),
                constants_path:
                    "data/pretrained_models/pretrained_models/alphapeptdeep/generic/ms2.pth.model_const.yaml"
                        .to_string(),
                architecture: "ms2_bert".to_string(),
            };
            log::info!("Pre-trained PeptDeep MS2 intensity model (architecture: {}): {}", ms2_intensity_model_config.architecture,ms2_intensity_model_config.model_path);
        }

        Self {
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
    pub peptide_chunking: ChunkingStrategy,
    pub output_file: Option<String>,
    pub write_report: Option<bool>,
    pub parquet_output: Option<bool>,
}

/// Documentation-only struct for generating JSON schema help
/// This mirrors the Input struct but only includes fields we control
#[derive(JsonSchema)]
#[schemars(title = "EasyPQP Configuration", description = "JSON configuration file for in-silico peptide library generation")]
#[allow(dead_code)]
struct InputSchema {
    /// Database and digestion parameters (REQUIRED)
    /// Fields: fasta (required), enzyme, peptide_min_mass, peptide_max_mass, 
    /// static_mods, variable_mods, max_variable_mods, generate_decoys, decoy_tag
    #[schemars(description = "FASTA database path and digestion settings (enzyme, modifications, decoys)")]
    database: DatabaseSchema,
    
    /// In-silico library generation settings (REQUIRED)
    #[schemars(description = "Precursor charges, fragment charges, transition limits, fragmentation model, and RT scaling")]
    insilico_settings: InsilicoPQPSettings,
    
    /// Deep learning feature prediction models (OPTIONAL)
    #[schemars(description = "Custom model paths for RT/IM/MS2 prediction. If omitted, pretrained AlphaPeptDeep models will be auto-downloaded")]
    dl_feature_generators: Option<DLFeatureGenerators>,
    
    /// Peptide chunking strategy for memory management
    #[schemars(description = "Number of peptides per chunk (default: 0 = auto-calculate based on available memory)")]
    peptide_chunking: ChunkingStrategy,
    
    /// Output file path
    #[schemars(description = "Path for output TSV file (default: 'insilico_library.tsv')")]
    output_file: Option<String>,
    
    /// Generate HTML report
    #[schemars(description = "Whether to generate an HTML quality report (default: true)")]
    write_report: Option<bool>,
    
    /// Output in Parquet format
    #[schemars(description = "Generate output in Parquet format instead of TSV (default: false)")]
    parquet_output: Option<bool>,
}

/// Schema representation of database configuration
#[derive(JsonSchema)]
#[schemars(description = "Database and protein digestion configuration")]
#[allow(dead_code)]
struct DatabaseSchema {
    /// Path to FASTA protein database file (REQUIRED)
    fasta: String,
    
    /// Generate decoy peptides (default: true)
    #[schemars(description = "Auto-generate decoy sequences")]
    generate_decoys: Option<bool>,
    
    /// Decoy protein tag/prefix (default: "rev_")
    #[schemars(description = "Prefix for decoy protein names")]
    decoy_tag: Option<String>,
    
    /// Minimum peptide mass in Daltons (default: 500.0)
    peptide_min_mass: Option<f32>,
    
    /// Maximum peptide mass in Daltons (default: 5000.0)
    peptide_max_mass: Option<f32>,
    
    /// Maximum number of variable modifications per peptide (default: 2)
    max_variable_mods: Option<usize>,
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

        // CLI flag to disable writing HTML report (default is to write)
        if matches.get_flag("no-write-report") {
            input.write_report = Some(false);
        }

        // CLI flag to enable Parquet output format
        if matches.get_flag("parquet") {
            input.parquet_output = Some(true);
        }

        // avoid to later panic if these parameters are not set (but doesn't check if files exist)

        ensure!(
            input.database.fasta.is_some(),
            "`database.fasta` must be set. For more information try '--help'"
        );

        Ok(input)
    }

    #[allow(dead_code)]
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
        rt_scale: Option<f32>,
        fine_tune: Option<bool>,
        train_data_path: Option<String>,
        save_model: Option<bool>,
        instrument: Option<String>,
        nce: Option<f32>,
        batch_size: Option<usize>,
        write_report: Option<bool>,
        parquet_output: Option<bool>,
    ) -> Result<Self> {
        // Check if it's a path to an existing file
        let mut input: Input = if Path::new(path).exists() {
            Input::load(path)
                .with_context(|| format!("Failed to read parameters from `{path}`"))?
        } else {
            serde_json::from_str(path)
                .with_context(|| "Failed to parse JSON configuration from string")?
        };

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
        if let Some(rt_scale) = rt_scale {
            input.insilico_settings.rt_scale = rt_scale;
        }
        if let Some(fine_tune) = fine_tune {
            input.dl_feature_generators.clone().unwrap().fine_tune_config.unwrap().fine_tune = fine_tune;
        }
        if let Some(train_data_path) = train_data_path {
            input.dl_feature_generators.clone().unwrap().fine_tune_config.unwrap().train_data_path = train_data_path;
        }
        if let Some(save_model) = save_model {
            input.dl_feature_generators.clone().unwrap().fine_tune_config.unwrap().save_model = save_model;
        }
        if let Some(instrument) = instrument {
            input.dl_feature_generators.clone().unwrap().instrument = Some(instrument);
        }
        if let Some(nce) = nce {
            input.dl_feature_generators.clone().unwrap().nce = Some(nce);
        }
        if let Some(batch_size) = batch_size {
            input.dl_feature_generators.clone().unwrap().fine_tune_config.unwrap().batch_size = batch_size;
        }
        if let Some(write_report) = write_report {
            input.write_report = Some(write_report);
        }
        if let Some(parquet_output) = parquet_output {
            input.parquet_output = Some(parquet_output);
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

    pub fn build(self) -> anyhow::Result<InsilicoPQP> {
        let database = self.database.make_parameters();
        
        let parquet_output = self.parquet_output.unwrap_or(false);
        let default_extension = if parquet_output { "parquet" } else { "tsv" };
        let output_file = self.output_file.clone().unwrap_or_else(|| {
            format!("insilico_library.{}", default_extension)
        });

        Ok(InsilicoPQP {
            version: clap::crate_version!().into(),
            database,
            insilico_settings: self.insilico_settings.clone(),
            dl_feature_generators: self
                .dl_feature_generators
                .map(Into::into)
                .unwrap_or_default(),
            peptide_chunking: self.peptide_chunking.clone(),
            output_file,
            write_report: self.write_report.unwrap_or(true),
            parquet_output,
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
    pub peptide_chunking: ChunkingStrategy,
    pub output_file: &'a String,
    pub write_report: bool,
}

// 3. Implement conversion methods
impl InsilicoPQP {
    pub fn as_serializable(&self) -> InsilicoPQPDto<'_> {
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
            peptide_chunking: self.peptide_chunking.clone(),
            output_file: &self.output_file,
            write_report: self.write_report,
        }
    }
}
/// Print detailed help about the JSON configuration file structure
/// This automatically generates help from the struct definitions and doc comments
pub fn print_config_help() {
    let schema = schema_for!(InputSchema);
    
    println!("═══════════════════════════════════════════════════════════════════════════════");
    println!("                  JSON CONFIGURATION FILE STRUCTURE");
    println!("═══════════════════════════════════════════════════════════════════════════════\n");
    
    if let Some(description) = &schema.schema.metadata.as_ref().and_then(|m| m.description.as_ref()) {
        println!("{}\n", description);
    }
    
    // Print schema in a readable format
    let json = serde_json::to_string_pretty(&schema).unwrap();
    println!("JSON Schema:");
    println!("{}\n", json);
    
    println!("═══════════════════════════════════════════════════════════════════════════════");
    println!("MINIMAL EXAMPLE:");
    println!("═══════════════════════════════════════════════════════════════════════════════");
    println!("{{");
    println!("  \"database\": {{");
    println!("    \"fasta\": \"proteins.fasta\"");
    println!("  }},");
    println!("  \"insilico_settings\": {{");
    println!("    \"precursor_charge\": [2, 3]");
    println!("  }}");
    println!("}}\n");
    
    println!("═══════════════════════════════════════════════════════════════════════════════");
    println!("FULL EXAMPLE:");
    println!("═══════════════════════════════════════════════════════════════════════════════");
    println!("{{");
    println!("  \"database\": {{");
    println!("    \"fasta\": \"proteins.fasta\",");
    println!("    \"generate_decoys\": true,");
    println!("    \"static_mods\": {{ \"C\": 57.021464 }}");
    println!("  }},");
    println!("  \"insilico_settings\": {{");
    println!("    \"precursor_charge\": [2, 3, 4],");
    println!("    \"fragmentation_model\": \"cid_hcd\",");
    println!("    \"rt_scale\": 100.0");
    println!("  }},");
    println!("  \"dl_feature_generators\": {{");
    println!("    \"instrument\": \"timsTOF\",");
    println!("    \"nce\": 25.0");
    println!("  }},");
    println!("  \"output_file\": \"my_library.tsv\"");
    println!("}}\n");
    
    println!("For more information: https://github.com/singjc/easypqp-rs\n");
}
