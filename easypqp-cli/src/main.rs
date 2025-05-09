use clap::{Arg, Command, ValueHint};
use input::Input;

mod input;
mod runner;
pub mod output;


fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("EASYPQP_LOG", "error,easypqp=info"))
        .init();

    let matches = Command::new("easypqp")
        .version(clap::crate_version!())
        .author("Justin Sing <justincsing@gmail.com>")
        .about("\u{1F52E} EasyPQP \u{1F9D9} - In-silico Peptide query parameter generation")
        .arg(
            Arg::new("parameters")
                .required(true)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Path to configuration parameters (JSON file)")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to FASTA database. Overrides the FASTA file \
                     specified in the configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("output_file")
                .short('o')
                .long("output_file")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "File path that in-silico library will be written to. \
                     Overrides the directory specified in the configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    // Not supported yet
    // let parquet = matches.get_one::<bool>("parquet").copied().unwrap_or(false);

    let input = Input::from_arguments(matches)?;

    let runner = input.build().and_then(runner::Runner::new)?;

    let _ = runner.run()?;

    Ok(())
}
