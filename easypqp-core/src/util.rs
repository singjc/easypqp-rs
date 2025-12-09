use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use anyhow::Error;
use sysinfo::System;

/// Reads a FASTA file and returns a Fasta object.
/// TODO: Replace with rustyms?
pub fn read_fasta<S>(
    path: S,
    decoy_tag: S,
    generate_decoys: bool,
) -> Result<sage_core::fasta::Fasta, Error>
where
    S: AsRef<str>,
{
    let path = Path::new(path.as_ref());
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    
    Ok(sage_core::fasta::Fasta::parse(
        contents,
        decoy_tag.as_ref(),
        generate_decoys,
    ))
}


pub fn read_json<S, T>(path: S) -> Result<T, Error>
where
    S: AsRef<str>,
    T: for<'de> serde::Deserialize<'de>,
{
    let path = Path::new(path.as_ref());
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(serde_json::from_str(&contents)?)
}



pub fn write_bytes_to_file(path: &str, bytes: &[u8]) -> std::io::Result<()> {
    let path = Path::new(path);
    let mut file = File::create(path)?;
    file.write_all(bytes)?;
    Ok(())
}


pub fn auto_chunk_size(peptide_bytes_estimate: usize, safety_ratio: f64) -> usize {
    let mut sys = System::new_all();
    sys.refresh_memory();

    let available_bytes = sys.free_memory() * 1024; 
    let safe_bytes = (available_bytes as f64 * safety_ratio) as usize;

    let chunk = safe_bytes / peptide_bytes_estimate;
    chunk.clamp(1000, 10_000_000) // minimum of 1000 peptides per chunk
}


pub fn get_test_file(name: &str) -> PathBuf {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));

    // Walk up to workspace root
    let workspace_root = manifest_dir.ancestors().nth(1).unwrap(); 

    workspace_root.join("test-data").join(name)
}