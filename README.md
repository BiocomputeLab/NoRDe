# NoRDe – Nonrepetitive RNA Designer

A modular, extensible Python toolkit for the **generation, filtering, grouping, analysis, and visualization** of RNA sequence variants under structural and evolutionary constraints.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command-Line Usage](#command-line-usage)
- [Configuration](#configuration)
- [Variant Generation Modes](#variant-generation-modes)
- [Main Workflow](#main-workflow)
- [Analysis & Visualization](#analysis--visualization)
- [API Reference](#api-reference)
- [Reproducibility](#reproducibility)
- [Testing](#testing)
- [License & Citation](#license--citation)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)

---

## Overview

This toolkit enables **systematic design and analysis of RNA sequence variants** with user-defined constraints on structure, diversity, and evolutionary conservation. It is suitable for computational biology research, synthetic biology, and RNA engineering projects.

---

## Features

- **Inverse folding** and **conservation-guided** variant generation
- **Flexible filtering** by structure, GC content, and homopolymer runs
- **Diversity-based grouping** and selection
- **Mutation tolerance and conservation analysis**
- **Comprehensive visualization** (heatmaps, diversity metrics)
- **Parallel processing** for scalability
- **Caching** for efficient repeated runs
- **Configurable workflow** via a central config file
- **Extensible modular design**
- **Command-line interface** for reproducible, scriptable workflows

---

## Installation

### Requirements

- Python 3.8+
- [ViennaRNA Python bindings](https://www.tbi.univie.ac.at/RNA/)
- numpy, matplotlib, seaborn, scikit-learn

### Quick Installation

```bash
# Clone the repository
git clone https://github.com/kkariotis/NoRDe
cd NoRDe

# Install dependencies
pip install -r requirements.txt

# Install ViennaRNA Python bindings (see their documentation)
# For macOS: brew install viennarna
# For Ubuntu: sudo apt-get install libviennarna-dev
# Then: pip install viennarna
```

### Development Installation

```bash
# Install in development mode
pip install -e .

# Install development dependencies
pip install -e ".[dev]"
```

---

## Quick Start

1. **Configure parameters** in `tool/config.py` as needed.
2. **Run the main workflow:**

```bash
python -m tool.main
```

Outputs (variants, groups, plots) will be saved in the `output/` directory.

---

## Command-Line Usage

A command-line interface (CLI) is provided for flexible, reproducible runs.
You can set all major parameters directly from the command line:

```bash
python -m tool.cli --mode inverse --n_runs 5 --output output/my_variants.fasta --group_size 10 --n_groups 2 --seed 42
```

**Arguments:**
- `--mode` : Variant generation mode (`auto`, `inverse`, `conservation`)
- `--n_runs` : Number of independent runs
- `--output` : Output FASTA file
- `--group_size` : Number of variants per group
- `--n_groups` : Number of groups
- `--seed` : Random seed for reproducibility

---

## Configuration

All workflow parameters are set in [`tool/config.py`](tool/config.py):

- **Structural constraints:**
  `WILD_TYPE_SCAFFOLD`, `LMAX_THRESHOLD`, `GC_CONTENT_RANGE`
- **Variant generation:**
  `GENERATION_METHOD` (`auto`, `inverse`, `conservation`), `AUTO_SWITCH_LENGTH`
- **Mutation parameters:**
  `CONSERVATION_ATTEMPTS`, `MUTATION_WEIGHTS`, `CONSERVATION_BIAS`
- **Parallelization:**
  `NUM_CORES`
- **Reproducibility:**
  `RANDOM_SEED` (set to integer for reproducible runs)
- **Output:**
  `OUTPUT_DIR`

---

## Variant Generation Modes

- **Inverse folding:** Generates variants matching a target structure.
- **Conservation-guided:** Uses mutation tolerance profiles to guide variant generation.
- **Auto:** Selects method based on sequence length (`AUTO_SWITCH_LENGTH`).

Set mode in `tool/config.py`:
```python
GENERATION_METHOD = "auto"  # or "inverse", "conservation"
```

---

## Main Workflow

The main workflow (`tool/main.py`) performs:

1. **Variant generation** (multiple runs for diversity)
2. **Filtering** (structure, GC, homopolymers, confidence)
3. **Lmax filtering** (limits shared subsequence length)
4. **Grouping** (maximizes diversity)
5. **Saving outputs** (FASTA files, plots)
6. **Diversity and conservation analysis**

---

## Analysis & Visualization

- **Difference heatmaps:** Nucleotide differences vs wild-type
- **Hamming/Lmax heatmaps:** Diversity metrics
- **Conservation analysis:** Mutation tolerance profiles
- **Group diversity metrics:** Intra/inter-group distances

All plots are saved in `output/`.

---

## API Reference

### Main Classes

- `Config`: Central configuration
- `VariantGenerator`: Generates RNA variants
- `VariantFilter`: Filters variants by constraints
- `GroupGenerator`: Groups variants for diversity
- `SequenceAnalyzer`: Structure, GC, homopolymer analysis
- `ConservationAnalyzer`: Mutation tolerance profiling
- `PlotGenerator`, `ConservationVisualizer`: Visualization tools

### Example: Generate and Filter Variants

```python
from tool.config import Config
from tool.core.variant_generation import VariantGenerator
from tool.core.filtering import VariantFilter

wt_struct, _ = SequenceAnalyzer.run_rnafold(Config.WILD_TYPE_SCAFFOLD)
variants = VariantGenerator.inverse_fold_variants(wt_struct, 100)
filtered = VariantFilter.filter_variants(variants, Config.WILD_TYPE_SCAFFOLD, wt_struct)
```

---

## Reproducibility

The toolkit provides comprehensive reproducibility features with important limitations:

### ⚠️ Critical Reproducibility Limitation

**ViennaRNA operations are NOT reproducible** due to C library limitations. This affects:
- RNA structure prediction results
- Inverse folding results
- Final variant sequences

**What IS reproducible:**
- Mutation positions and base choices
- Clustering algorithms (K-means)
- Random sampling operations
- NumPy and Python operations

**What is NOT reproducible:**
- ViennaRNA structure predictions
- ViennaRNA inverse folding
- Final variant sequences (due to ViennaRNA dependency)

### Random Seed Configuration

**Via Config File:**
```python
# In tool/config.py
RANDOM_SEED = 42  # Set to any integer for reproducible runs
```

**Via Command Line:**
```bash
python -m tool.cli --seed 42 --mode conservation --n_runs 5
```

**Via Code:**
```python
from tool.config import Config
Config.set_seed(42)  # Sets seed for all random operations
```

### Reproducibility Features

- **Partial reproducibility:** Mutation and clustering operations are seeded
- **Flexible seeding:** Can be set via config, CLI, or programmatically
- **Caching:** Results are cached for repeated analysis
- **Config file:** All parameters are centralized for easy sharing

### Workarounds for Full Reproducibility

1. **Cache ViennaRNA results** for repeated runs
2. **Document ViennaRNA version** and system environment
3. **Use alternative RNA tools** with seeding support
4. **Implement custom RNA folding** (slower but seedable)

---

## Testing

The toolkit includes comprehensive tests to ensure code quality and functionality:

```bash
# Run all tests
python test_publication_ready.py

# Run with pytest (if installed)
pytest tests/

# Run with coverage
pytest --cov=tool tests/
```

### Code Quality

The codebase follows publication-ready standards:
- **Type hints** throughout all functions
- **Comprehensive error handling** with specific exception types
- **Proper logging** instead of print statements
- **Documentation** for all public functions and classes
- **Modular design** with clear separation of concerns

---

## License & Citation

This toolkit is released under the MIT License.

If you use this toolkit in published research, please cite:

```
.....
```

---

## Contact

For questions, bug reports, or feature requests, please contact:

- Konstantinos Kariotis
- Email: kariotis.konst@gmail.com

---

## Acknowledgements

- ViennaRNA Package authors
    - Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. Algorithms for Molecular Biology, 6(1), 1-12.
- Open-source contributors

---

**Enjoy designing and analyzing RNA variants!**

---
