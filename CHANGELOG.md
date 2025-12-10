# Changelog

All notable changes to this project will be documented in this file.

## [0.1.2] - 2025-12-10

### 🐛 Bug Fixes

- Add manylinux support and allow manual workflow dispatch for PyPI publishing
- Update package version to 0.1.2 in Cargo.toml and pyproject.toml files
- Improve sccache setup and build steps for better clarity and organization

### ⚙️ Miscellaneous Tasks

- Update CHANGELOG.md
- Update CHANGELOG.md

## [0.1.1] - 2025-12-10

### 🚀 Features

- Enhance project metadata in pyproject.toml
- Add PyPI version badge to README.md
- Enhance GitHub Actions workflow for release tagging and binary uploads

### 🐛 Bug Fixes

- Update binary name for Windows build in rust-release workflow
- Switch Rust toolchain to stable and refine build/test commands
- Update Python interpreter versions in build args and classifiers
- Bump version to 0.1.1 in Cargo.toml and pyproject.toml files
- Refine conditions for release job and asset uploads in GitHub Actions workflow
- Add tag reference for binary uploads in GitHub Actions workflow
- Move OpenSSL environment variables to non-musl upload step
- Set OpenSSL environment variables for musl builds in GitHub Actions workflow
- Remove OpenSSL environment setup for musl and add vendored feature in Cargo.toml
- Add package reference for uploaded binaries in rust-release workflow
- Update model paths and architectures for pretrained models in DLFeatureGeneratorSettings
- Update development status classifier from 1 to 3 - Alpha in pyproject.toml
- Configure sccache for faster Rust compilation in Python wheel build

### ⚙️ Miscellaneous Tasks

- Update CHANGELOG.md
- Update CHANGELOG.md

## [0.1.0-alpha] - 2025-12-10

### 🚀 Features

- Add easypqp-py module for Python bindings
- Add save_model, instrument, and nce options for generating in-silico libraries
- Update peptide data length assertion in tests
- Add test data
- Add ChunkingStrategy for peptide chunking in input.rs
- Enable peptide chunking in config.json
- Add GitHub Actions badge to README.md
- Add HTML report generation after TSV output in runner
- Add option to generate HTML report after TSV output
- Change CLI flag to disable HTML report generation
- Add configurable retention time scale to InsilicoPQPSettings and update output functions
- Update Parquet output support and modify related CLI flags
- Add default retention time scale to InsilicoPQPSettings
- Integrate schemars for JSON schema generation in core and CLI modules
- Enhance documentation for training data path in FineTuneSettings
- Add GitHub Actions workflow for automatic changelog generation
- Add optional parameters for RT scaling and report generation in insilico library
- Refactor GitHub Actions workflow for building and publishing Python wheels
- Add threads parameter for parallel processing in generate_insilico_library function
- Add README.md with installation and development instructions for easypqp-py
- Add concurrency settings to GitHub Actions workflows for better job management
- Add Dockerfile and .dockerignore for building and publishing easypqp-insilico

### 🐛 Bug Fixes

- Remove emoji from about description in easypqp CLI
- Enhance musl-tools installation step with additional dependencies
- Update Docker image tag prefix and enhance Python wheel build arguments
- Add dead code allowance for InputSchema and DatabaseSchema structs
- Set OPENSSL_STATIC and OPENSSL_VENDORED environment variables for musl targets in release workflow
- Add missing newline at end of Cargo.toml and ensure proper formatting
- Remove redundant build steps for easypqp-insilico in rust-release workflow
- Specify manifest path for building wheels in python-publish workflow

### 💼 Other

- Readme

### 🚜 Refactor

- Property peptide prediction
- Optimize property prediction
- Remove unused code in runner.rs
- Update dependencies and improve data handling in property prediction and tuning data
- Update redeem-properties dependency to use develop branch and enhance retention time scaling logic
- Simplify retention time scaling logic and enhance protein entry parsing in TSV output
- Remove unused predict_properties field from DLFeatureGenerators and DLFeatureGeneratorSettings

### 📚 Documentation

- Update documentation for read_fasta function
- Add contributing guidelines for commit message conventions

### ⚙️ Miscellaneous Tasks

- Add .gitignore and Cargo.toml files for project setup
- Add util module for reading fasta and json files and optimize property prediction
- Update dependencies and improve parameter handling for generating in-silico libraries
- Add tests and lint code
- Add GitHub Actions workflows for Rust and Python publishing
- Rename binary to "easypqp-insilico" in Cargo.toml
- Update dependencies and add chrono crate for time handling
- Update binary names and build configurations for different platforms
- Update dependencies and add logging for error handling in python wrapped method
- Update value hint for file path in main.rs
- Update CHANGELOG.md
- Update CHANGELOG.md
- Update CHANGELOG.md
- Update CHANGELOG.md

<!-- generated by git-cliff -->
