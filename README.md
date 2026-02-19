# easypqp-rs

---

[![Rust](https://github.com/singjc/easypqp-rs/actions/workflows/rust.yml/badge.svg?branch=master)](https://github.com/singjc/easypqp-rs/actions/workflows/rust.yml)

easypqp-rs is a Rust library for in-silico peptide library generation, with  Python bindings for integration with the python [EasyPQP](https://github.com/grosenberger/easypqp) library.

## Features

* Fast in-silico library generation using Rust

* Includes a command-line tool for batch library generation

* Python bindings for integration within the easypqp Python package

* Configurable via JSON for fine-tuning predictions, fragmentation settings, and NCE/instrument profiles

## Rust Binary CLI Example

easypqp-rs has an optional standalone command-line interface (CLI) binary for generating in-silico libraries. This can be used independently of the EasyPQP Python package if you prefer.

```bash

easypqp-insilico ./config.json
```

## Configuration Reference

The tool is configured via a JSON file. Below is a comprehensive guide to all available parameters.

### Complete Example Configuration

<details>
<summary>Click to expand full example config.json</summary>

```json
{
  "database": {
    "fasta": "path/to/proteins.fasta",
    "enzyme": {
      "name": "Trypsin/P",
      "cleave_at": "KR",
      "restrict": "P",
      "c_terminal": null,
      "min_len": 7,
      "max_len": 50,
      "missed_cleavages": 2
    },
    "peptide_min_mass": 500.0,
    "peptide_max_mass": 5000.0,
    "generate_decoys": true,
    "decoy_tag": "rev_",
    "static_mods": {
      "C": 57.0215
    },
    "variable_mods": {
      "M": [15.9949],
      "[": [42.0106]
    },
    "max_variable_mods": 2
  },
  "insilico_settings": {
    "precursor_charge": [2, 3, 4],
    "max_fragment_charge": 2,
    "min_transitions": 6,
    "max_transitions": 6,
    "fragmentation_model": "HCD",
    "allowed_fragment_types": ["b", "y"],
    "rt_scale": 100.0
  },
  "dl_feature_generators": {
    "retention_time": {
      "model_path": "path/to/rt_model.safetensors",
      "constants_path": "path/to/rt_model_const.yaml",
      "architecture": "rt_cnn_tf"
    },
    "ion_mobility": {
      "model_path": "path/to/ccs_model.safetensors",
      "constants_path": "path/to/ccs_model_const.yaml",
      "architecture": "ccs_cnn_tf"
    },
    "ms2_intensity": {
      "model_path": "path/to/ms2_model.pth",
      "constants_path": "path/to/ms2_model_const.yaml",
      "architecture": "ms2_bert"
    },
    "device": "cpu",
    "instrument": "timsTOF",
    "nce": 20.0,
    "batch_size": 64,
    "fine_tune_config": {
      "fine_tune": false,
      "train_data_path": "",
      "batch_size": 256,
      "epochs": 3,
      "learning_rate": 0.001,
      "save_model": false
    }
  },
  "peptide_chunking": 0,
  "output_file": "insilico_library.tsv",
  "write_report": true,
  "parquet_output": false
}
```

</details>

### Configuration Sections

#### 1. Database Settings (REQUIRED)

<details>
<summary><code>database</code> - FASTA file, enzyme, modifications, and decoy generation</summary>

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fasta` | string | **REQUIRED** | Path to FASTA protein database file |
| `generate_decoys` | boolean | `true` | Auto-generate decoy sequences by reversing protein sequences |
| `decoy_tag` | string | `"rev_"` | Prefix added to decoy protein names |
| `peptide_min_mass` | number | `500.0` | Minimum peptide mass in Daltons |
| `peptide_max_mass` | number | `5000.0` | Maximum peptide mass in Daltons |
| `max_variable_mods` | integer | `2` | Maximum number of variable modifications per peptide |

**Enzyme Configuration:**
```json
"enzyme": {
  "name": "Trypsin/P",          // Enzyme name (for reference)
  "cleave_at": "KR",             // Amino acids where enzyme cleaves
  "restrict": "P",               // Amino acid that prevents cleavage if following cleavage site
  "c_terminal": true,            // Cleavage occurs C-terminal to the cleavage site
  "min_len": 7,                  // Minimum peptide length
  "max_len": 50,                 // Maximum peptide length
  "missed_cleavages": 2          // Number of allowed missed cleavages
}
```

**Static Modifications:**
```json
"static_mods": {
  "C": 57.0215    // Carbamidomethylation of Cysteine
}
```

**Variable Modifications:**
```json
"variable_mods": {
  "M": [15.9949],    // Oxidation of Methionine
  "[": [42.0106]     // N-terminal Acetylation
}
```

Common modification masses:
- Carbamidomethyl (C): `57.0215`
- Oxidation (M): `15.9949`
- Phosphorylation (STY): `79.9663`
- N-terminal Acetylation: `42.0106`
- Deamidation (NQ): `0.9840`

</details>

#### 2. In-Silico Library Settings (REQUIRED)

<details>
<summary><code>insilico_settings</code> - Precursor/fragment charges, transitions, fragmentation model</summary>

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `precursor_charge` | array[int] | `[2, 3, 4]` | Precursor charge states to generate |
| `max_fragment_charge` | integer | `2` | Maximum fragment ion charge |
| `min_transitions` | integer | `6` | Minimum number of transitions per precursor |
| `max_transitions` | integer | `6` | Maximum number of transitions per precursor |
| `fragmentation_model` | string | `"HCD"` | Fragmentation type: `"HCD"`, `"CID"`, or `"ETD"` |
| `allowed_fragment_types` | array[string] | `["b", "y"]` | Allowed fragment ion types: `"b"`, `"y"` |
| `rt_scale` | number | `100.0` | Retention time scaling factor (multiplies predicted RT) |

> [!NOTE]
> The current MS2 intensity prediction models only support `"b"` and `"y"` fragment ions. 

**Example:**
```json
"insilico_settings": {
  "precursor_charge": [2, 3],
  "max_fragment_charge": 1,
  "min_transitions": 6,
  "max_transitions": 12,
  "fragmentation_model": "HCD",
  "allowed_fragment_types": ["b", "y"],
  "rt_scale": 1.0
}
```

</details>

#### 3. Deep Learning Models (OPTIONAL)

> [!NOTE]
> If no `retention_time`, `ion_mobility`, or `ms2_intensity` fields are provided under `dl_feature_generators`, pretrained models will be automatically downloaded and used. The current default pretrained models used are:
> - RT: `rt_cnn_tf` - A CNN-Transformer model trained on the [ProteomicsML repository RT dataset](https://proteomicsml.org/datasets/retentiontime/ProteomeTools_RT.html). This model is based on AlphaPeptDeep's CNN-LSTM implementation, with the LSTM replaced by a Transformer encoder.
> - CCS: `ccs_cnn_tf` - A CNN-Transformer model trained on the [ProteomicsML repository CCS dataset](https://proteomicsml.org/datasets/ionmobility/Meier_TIMS.html). This model is also based on AlphaPeptDeep's CNN-LSTM implementation, with the LSTM replaced by a Transformer encoder.
> - MS2: `ms2_bert` - A BERT-based model retreived from AlphaPeptDeep's pretrained models.

<details>
<summary><code>dl_feature_generators</code> - Custom or pretrained RT/IM/MS2 prediction models</summary>

If this section is **omitted** or **empty**, pretrained AlphaPeptDeep models will be automatically downloaded and used.

**Model Configuration:**

Each model (RT, IM, MS2) requires three files:
```json
{
  "model_path": "path/to/model.safetensors",     // Model weights (.pth or .safetensors)
  "constants_path": "path/to/model_const.yaml",  // Model configuration constants
  "architecture": "model_architecture_name"       // Architecture identifier
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `retention_time` | object | *pretrained* | Custom RT prediction model |
| `ion_mobility` | object | *pretrained* | Custom IM/CCS prediction model (timsTOF only) |
| `ms2_intensity` | object | *pretrained* | Custom MS2 intensity prediction model |
| `device` | string | `"cpu"` | Compute device: `"cpu"`, `"cuda"`, or `"mps"` (Apple Silicon) |
| `instrument` | string | `"timsTOF"` | Instrument type: `"QE"` or `"timsTOF"` |
| `nce` | number | `20.0` | Normalized collision energy for fragmentation |
| `batch_size` | integer | `64` | Batch size for model inference |
| `fine_tune_config` | object | *see below* | Optional fine-tuning configuration |

**Supported Architectures:**
- RT: `"rt_cnn_tf"`, `"rt_lstm"`, `"rt_transformer"`
- IM/CCS: `"ccs_cnn_tf"`, `"ccs_lstm"`  
- MS2: `"ms2_bert"`, `"ms2_transformer"`

</details>

#### 4. Fine-Tuning (OPTIONAL)

<details>
<summary><code>fine_tune_config</code> - Transfer learning on experimental data</summary>

Fine-tune pretrained models on your own experimental data for improved accuracy.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fine_tune` | boolean | `false` | Enable fine-tuning |
| `train_data_path` | string | `""` | Path to training data TSV file |
| `batch_size` | integer | `256` | Training batch size |
| `epochs` | integer | `3` | Number of training epochs |
| `learning_rate` | number | `0.001` | Learning rate for optimizer |
| `save_model` | boolean | `false` | Save fine-tuned model weights to disk |

**Training Data Format (TSV):**

Required columns:
- `sequence`: Modified sequence with square bracket notation (e.g., `MGC[+57.0215]AAR`)
- `precursor_charge`: Precursor charge state
- `retention_time`: Experimental retention time
- `ion_mobility`: CCS value (only if using timsTOF)
- `fragment_type`: Fragment ion type (`b`, `y`, etc.)
- `fragment_series_number`: Fragment position
- `product_charge`: Fragment charge
- `intensity`: Normalized fragment intensity

**Example:**
```json
"fine_tune_config": {
  "fine_tune": true,
  "train_data_path": "experimental_data.tsv",
  "batch_size": 256,
  "epochs": 5,
  "learning_rate": 0.0001,
  "save_model": true
}
```

</details>

#### 5. Output Settings (OPTIONAL)

<details>
<summary>Output file format, reporting, and memory management</summary>

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | string | `"insilico_library.tsv"` | Path for output library file |
| `write_report` | boolean | `true` | Generate HTML quality control report |
| `parquet_output` | boolean | `false` | Output in Parquet format instead of TSV |
| `peptide_chunking` | integer | `0` | Peptides per chunk (0 = auto-calculate based on memory) |

**Peptide Chunking:**
- `0` (default): Automatically calculate chunk size based on available memory (recommended)
- `> 0`: Manual chunk size for processing large FASTA files with limited RAM
- Larger chunks = faster processing but more memory usage

</details>

### Minimal Configuration

The minimum required configuration only needs a FASTA file:

```json
{
  "database": {
    "fasta": "proteins.fasta"
  }
}
```

All other parameters will use sensible defaults and pretrained models will be auto-downloaded.

### Command-Line Overrides

You can override JSON configuration values via command-line arguments:

```bash
easypqp-insilico config.json \
  --fasta my_proteins.fasta \
  --output_file my_library.tsv \
  --no-write-report \
  --parquet
```

**Available flags:**
- `--fasta <PATH>`: Override database FASTA file
- `--output_file <PATH>`: Override output file path  
- `--no-write-report`: Disable HTML report generation
- `--parquet`: Output in Parquet format instead of TSV


