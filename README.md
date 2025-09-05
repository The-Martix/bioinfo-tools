# Bioinfo Tools

A collection of bioinformatics tools including structools, aminotools, seqtools, chemtools, and more.

## Installation

This package can be installed using pip directly from the source directory or from a Python package repository.

### Install from Source (Recommended for Development)

1. Clone or download this repository
2. Navigate to the project directory
3. Install in development mode:

```bash
pip install -e .
```

This installs the package in "editable" mode, meaning changes to the source code are immediately reflected without reinstalling.

### Install from Source (Standard Installation)

```bash
pip install .
```

### Install with Development Dependencies

To install with additional development tools (pytest, black, flake8):

```bash
pip install -e ".[dev]"
```

## Requirements

- Python 3.8 or higher
- Dependencies are automatically installed with pip:
  - requests
  - matplotlib
  - seaborn
  - rdkit
  - biopython
  - numpy
  - scipy

## Usage

After installation, you can import any of the tool modules:

```python
from structools import *
from aminotools import *
from seqtools import *
from chemtools import *
```

## Alternative Installation (PYTHONPATH Method - Legacy)

If you prefer not to use pip installation, you can still use the older PYTHONPATH method:

### 🪟 Windows
```cmd
setx PYTHONPATH "C:\path\to\bioinfo-tools"
```

### 🐧 Linux / 🍎 macOS
Add to your `~/.bashrc`, `~/.zshrc`, or `~/.profile`:
```bash
export PYTHONPATH="/path/to/bioinfo-tools:$PYTHONPATH"
```

Then run:
```bash
source ~/.bashrc  # or your shell config file
```
