use pyo3::prelude::*;
// use easypqp_core::{InsilicoPQP, runner::Runner};
use easypqp_cli::input::Input;
use easypqp_cli::runner::Runner;

#[pyfunction]
fn generate_insilico_library(parameters: String, batch_size: Option<usize>) -> PyResult<()> {
    // Parse parameters from JSON string
    let params: Input = Input::load(&parameters)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Invalid parameters: {}", e)))?;

    // Create and run the runner
    let runner: Runner = Runner::new(params.build()?)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    
    let parallel: usize = batch_size.unwrap_or_else(|| num_cpus::get() / 2);
    runner.run(parallel, false)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn easypqp_rs(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(generate_insilico_library, m)?)?;
    Ok(())
}