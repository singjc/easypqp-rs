# easypqp-py: Python Bindings for EasyPQP

Python bindings for [EasyPQP rust library](https://github.com/justinsing/easypqp-rs). Currently, the rust library is mainly used for in-silico peptide query parameter generation.

## Prerequisites

### System Requirements

- **Rust**: 1.70 or newer
- **Python**: 3.10 or newer
- **Cargo**: Latest stable version
- **pip**: Python package manager

### Optional (Linux only)
For optimal binary compatibility on Linux, install `patchelf`:
```bash
# Debian/Ubuntu
sudo apt-get install patchelf

# Arch Linux
sudo pacman -S patchelf

# Via pip (alternative)
pip install maturin[patchelf]
```

## Installation

### Option 1: Development Installation (Editable Mode)

```bash
# Navigate to the easypqp-py directory
cd easypqp-py

# Install in development mode
maturin develop

# Or with optimizations enabled
maturin develop --release
```

### Option 2: Build and Install from Source

```bash
cd easypqp-py

# Build the wheel
maturin build

# Install the built wheel
pip install target/wheels/easypqp_rs-*.whl
```

### Option 3: Install via pip (when published)

```bash
pip install easypqp-rs
```

## Development

### Setting Up Development Environment

```bash
# Clone the repository
git clone https://github.com/justinsing/easypqp-rs.git
cd easypqp-rs/easypqp-py

# Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install maturin
pip install maturin

# Set up Python virtual environment (optional but recommended)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in development mode
maturin develop
```

### Common Build Commands

```bash
# Build in debug mode (faster builds)
maturin develop

# Build in release mode (optimized)
maturin develop --release

# Clean and rebuild
maturin develop --clean

# Build with specific features
maturin develop --features parquet

# Skip building dependencies
maturin develop --skip-install
```

### Testing the Installation

```python
import easypqp_rs

# Test basic functionality
print(easypqp_rs.__version__)  # Check if module loads
```