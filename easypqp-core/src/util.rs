use std::fs::File;
use std::io::Read;
use std::path::Path;
use anyhow::Error;
use std::path::PathBuf;


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


pub fn get_test_file(name: &str) -> PathBuf {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));

    // Walk up to workspace root
    let workspace_root = manifest_dir.ancestors().nth(1).unwrap(); 

    workspace_root.join("test-data").join(name)
}