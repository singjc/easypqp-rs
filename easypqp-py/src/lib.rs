use pyo3::prelude::*;
use easypqp_cli::input::Input;
use easypqp_cli::runner::Runner;

/// Generate an in-silico library using EasyPQP-rs
/// 
/// # Arguments
/// * `parameters`: JSON string containing the parameters for the library generation.
/// * `fasta`: Optional FASTA file path. Overrides the FASTA file specified in the parameters.
/// * `output_file`: Optional output file path. Overrides the output directory specified in the parameters.
/// * `generate_decoys`: Optional flag to generate decoy peptides. Overrides the JSON field parameter.
/// * `decoy_tag`: Optional tag for decoy peptides. Overrides the JSON field parameter.
/// * `precursor_charge`: Optional list of precursor charges. Overrides the JSON field parameter.
/// * `max_fragment_charge`: Optional maximum fragment charge. Overrides the JSON field parameter.
/// * min_transitions: Optional minimum number of transitions. Overrides the JSON field parameter.
/// * max_transitions: Optional maximum number of transitions. Overrides the JSON field parameter.
/// * `fragmentation_model`: Optional fragmentation model. Overrides the JSON field parameter.
/// * `allowed_fragment_types`: Optional list of allowed fragment types. Overrides the JSON field parameter.
/// * `fine_tune`: Optional flag to enable fine-tuning. Overrides the JSON field parameter.
/// * `train_data_path`: Optional path to training data. Overrides the JSON field parameter.
/// * `batch_size`: Optional batch size for parallel processing. Overrides the JSON field parameter.
/// 
/// # Returns
/// * `Ok(())` if the library generation is successful.
#[pyfunction]
fn generate_insilico_library(
    parameters: String,
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
) -> PyResult<()> {
    // Parse parameters from JSON string
    let params: Input = Input::from_passed_arguments(
        &parameters,
        fasta,
        output_file,
        generate_decoys,
        decoy_tag,
        precursor_charge,
        max_fragment_charge,
        min_transitions,
        max_transitions,
        fragmentation_model,
        allowed_fragment_types,
        fine_tune,
        train_data_path,
        batch_size,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Invalid parameters: {}", e))
    })?;

    // Create and run the runner
    let runner: Runner = Runner::new(params.build()?)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    runner.run()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn easypqp_rs(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(generate_insilico_library, m)?)?;
    Ok(())
}
