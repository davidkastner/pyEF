![Graphical Summary of README](docs/_static/header.webp)
==============================
# PyEF
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/davidkastner/pyef/workflows/CI/badge.svg)](https://github.com/davidkastner/pyef/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/pyef/badge/?version=latest)](https://pyef.readthedocs.io/en/latest/?badge=latest)

## Table of Contents
1. **Overview**
2. **Quick Start Guide**
    * Running Your First Analysis in 3 Steps
    * Common Analysis Types
    * Key Parameters
    * Python API Quick Start
3. **Detailed Tutorials**
    * Using PyEF Programmatically (Python API)
    * Using the Command-Line Interface (CLI)
    * Example Workflows
4. **Installation**
    * Download the package from GitHub
    * Creating a python environment
    * Developer install of PyEF
    * Supporting installations
5. **What is included?**
    * File structure
    * Command Line Interface
6. **Documentation**
    * Update the ReadTheDocs
    * GitHub refresher


## 1. Overview
The purpose of PyEF is to make electric field and electrostatics calculations more accessible.
PyEF is a python package optimized to run using molden files from QM calculations.

## 2. Quick Start Guide

### Running Your First Analysis in 3 Steps

**Step 1: Create a job list** (`jobs.csv`)
```csv
path/to/structure.pdb, 25, 26
```
Format: `structure_path, metal_atom_index, bonded_atom_index` (0-indexed)

**Step 2: Create a config file** (`config.yaml`)
```yaml
input: jobs.csv
dielectric: 1
ef: true                              # Electric field analysis
esp: false                            # Electrostatic potential
estab: false                          # Electrostatic stabilization
multiwfn_module: multiwfn
multiwfn_path: /path/to/multiwfn
atmrad_path: /path/to/atmrad
```

**Step 3: Run PyEF**
```bash
pyef -c config.yaml
```

**Output:** Results saved to `jobs_Efielddata.csv`

### Common Analysis Types

**Electric Field on Metal-Ligand Bonds:**
```yaml
ef: true
esp: false
estab: false
```

**Electrostatic Stabilization Energy:**
```yaml
ef: false
esp: false
estab: true
substrate_idxs: [1, 2, 3, 4, 5]  # Your substrate atoms (0-indexed)
```

**Complete Electrostatic Analysis:**
```yaml
ef: true
esp: true
estab: true
substrate_idxs: [1, 2, 3, 4, 5]
multipole_order: 2                # 1=monopole, 2=+dipole, 3=+quadrupole
```

### Key Parameters

| Parameter | Description | Common Values |
|-----------|-------------|---------------|
| `dielectric` | Medium dielectric constant | 1 (vacuum), 4 (protein), 78.5 (water) |
| `charge_types` | Charge partitioning scheme | `Hirshfeld_I` (recommended), `CHELPG`, `Hirshfeld` |
| `multipole_order` | Expansion order | 1 (charges), 2 (+dipoles), 3 (+quadrupoles) |
| `decompose_atomwise` | Per-atom contributions | `true` or `false` |

### Python API Quick Start

```python
from pyef.analysis import Electrostatics

# Initialize
es = Electrostatics(['structure/'], [25], '/scr/', dielectric=4.0)
es.prepData()
es.fix_allECPmolden()

# Calculate E-field
df = es.getEfield('Hirshfeld_I', 'output', 'multiwfn',
                  '/path/to/multiwfn', '/path/to/atmrad',
                  input_bond_indices=[(25, 26)])

# Calculate stabilization
estab_df = es.get_Electrostatic_stabilization(
    '/path/to/multiwfn', 'multiwfn', '/path/to/atmrad',
    substrate_idxs=[1,2,3,4,5], multipole_order=2
)
```

**For detailed examples and advanced options, see the Detailed Tutorials section below.**

---

## 3. Detailed Tutorials

### Using PyEF Programmatically (Python API)
The most common way to use PyEF is by importing it directly in your Python scripts:

```python
from pyef.analysis import Electrostatics

# Initialize the Electrostatics object
folders = ['structure1/', 'structure2/']
metal_indices = [25, 30]  # 0-indexed atom positions
path_to_molden = '/scr/'

es = Electrostatics(folders, metal_indices, path_to_molden, dielectric=4.0)

# Prepare data
es.prepData()
es.fix_allECPmolden()

# Calculate ESP
esp_df = es.getESP(
    charge_types=['Hirshfeld_I'],
    ESPdata_filename='esp_output',
    multiwfn_module='multiwfn',
    multiwfn_path='/path/to/multiwfn',
    atmrad_path='/path/to/atmrad',
    use_multipole=True
)

# Calculate E-field
efield_df = es.getEfield(
    charge_types='Hirshfeld_I',
    Efielddata_filename='efield_output',
    multiwfn_module='multiwfn',
    multiwfn_path='/path/to/multiwfn',
    atmrad_path='/path/to/atmrad',
    input_bond_indices=[(25, 26), (30, 31)],
    multipole_bool=True
)

# Calculate Electrostatic Stabilization
estab_df = es.get_Electrostatic_stabilization(
    multiwfn_path='/path/to/multiwfn',
    multiwfn_module='multiwfn',
    atmrad_path='/path/to/atmrad',
    substrate_idxs=[1, 2, 3, 4, 5],  # Substrate atoms
    env_idxs=[6, 7, 8, 9, 10],       # Environment atoms
    charge_type='Hirshfeld_I',
    multipole_order=2
)
```

### Using the Command-Line Interface (CLI)

PyEF provides a CLI for batch processing multiple structures, which is ideal for high-throughput analysis and workflow automation.

#### Step 1: Create a Configuration File

Create a YAML configuration file (e.g., `pyef_config.yaml`) with your analysis parameters:

```yaml
# ==========================================
# PyEF Configuration File
# ==========================================

# Input CSV file containing job information (required)
input: jobs.csv

# Dielectric constant for the medium (default: 1 for vacuum)
# Common values: 1 (vacuum), 4 (protein), 78.5 (water)
dielectric: 1

# ==========================================
# Analysis Types (set to true to enable)
# ==========================================

# Electric field analysis
ef: true

# Electrostatic potential (ESP) analysis
esp: false

# Electrostatic stabilization analysis
estab: false

# Geometry checking
geometry_check: false

# ==========================================
# Multiwfn Configuration (required)
# ==========================================

multiwfn_module: multiwfn
multiwfn_path: /path/to/multiwfn
atmrad_path: /path/to/atmrad

# ==========================================
# Advanced Options
# ==========================================

# Charge partitioning schemes
# Options: Hirshfeld, Hirshfeld_I, Voronoi, Mulliken, Lowdin,
#          Becke, ADCH, CHELPG, MK, AIM, CM5, RESP, etc.
charge_types:
  - Hirshfeld_I

# Use multipole expansion (monopole + dipole + quadrupole)
multipole: true           # For E-field calculations
use_multipole: false      # For ESP calculations

# Atom-wise decomposition (shows contribution from each atom)
decompose_atomwise: false

# Multipole expansion order (1=monopole, 2=+dipole, 3=+quadrupole)
multipole_order: 2

# ==========================================
# Electrostatic Stabilization Options
# (Only needed if estab: true)
# ==========================================

# Substrate atom indices (the molecule being stabilized)
# Use 0-based indexing
substrate_idxs:
  - 1
  - 2
  - 3
  - 4
  - 5

# Environment atom indices (optional - if not specified, uses all non-substrate atoms)
env_idxs:
  - 10
  - 11
  - 12
```

#### Step 2: Create a Job Batch File

Create a CSV file (e.g., `jobs.csv`) listing all structures to analyze:

```csv
# ==========================================
# PyEF Job Batch File
# ==========================================
# Format: job_path, metal_index, bonded_atom_index
# Lines starting with '#' are comments and will be ignored
# Use 0-based indexing for atom indices

structures/complex1.pdb, 25, 26
structures/complex2.pdb, 30, 31
structures/complex3.pdb, 28, 29

# You can add comments anywhere
structures/complex4.pdb, 32, 33  # This is an inline comment
```

**Column descriptions:**
- **Column 1**: Path to the structure file (PDB/molden format)
- **Column 2**: Index of the metal atom (0-based indexing)
- **Column 3**: Index of the bonded atom (0-based indexing)

#### Step 3: Run the CLI

Execute PyEF with your configuration file:

```bash
pyef -c pyef_config.yaml
```

Or using the long form:

```bash
pyef --config pyef_config.yaml
```

#### CLI Output

The CLI will:
1. Display a welcome banner with PyEF information
2. Parse your configuration and job batch file
3. Process each structure sequentially
4. Perform the requested analyses (ESP, E-field, electrostatic stabilization)
5. Save results to CSV files in the working directory

**Output files:**
- `{job_name}_ESPdata.csv` - ESP analysis results (if esp: true)
- `{job_name}_Efielddata.csv` - E-field analysis results (if ef: true)
- `{job_name}_Estab.csv` - Electrostatic stabilization results (if estab: true)

#### Example Workflows

**Example 1: Electric Field Analysis Only**

`config_efield.yaml`:
```yaml
input: protein_complexes.csv
dielectric: 4.0
ef: true
esp: false
estab: false
multiwfn_module: multiwfn
multiwfn_path: /usr/local/bin/multiwfn
atmrad_path: /usr/local/share/atmrad
charge_types:
  - Hirshfeld_I
multipole: true
```

`protein_complexes.csv`:
```csv
# Active site analysis
proteins/1abc_active_site.pdb, 150, 151
proteins/2def_active_site.pdb, 145, 146
```

Run: `pyef -c config_efield.yaml`

**Example 2: Electrostatic Stabilization Analysis**

`config_estab.yaml`:
```yaml
input: substrate_environment.csv
dielectric: 1
ef: false
esp: false
estab: true
multiwfn_module: multiwfn
multiwfn_path: /usr/local/bin/multiwfn
atmrad_path: /usr/local/share/atmrad
charge_types:
  - Hirshfeld_I
multipole_order: 2
substrate_idxs:  # Substrate molecule atoms
  - 1
  - 2
  - 3
  - 4
  - 5
# env_idxs not specified - will use all non-substrate atoms
decompose_atomwise: true  # Get per-atom contributions
```

Run: `pyef -c config_estab.yaml`

**Example 3: Complete Analysis (ESP + E-field + Stabilization)**

`config_complete.yaml`:
```yaml
input: full_analysis.csv
dielectric: 4.0
ef: true
esp: true
estab: true
geometry_check: true
multiwfn_module: multiwfn
multiwfn_path: /usr/local/bin/multiwfn
atmrad_path: /usr/local/share/atmrad
charge_types:
  - Hirshfeld_I
  - Hirshfeld
multipole: true
use_multipole: true
multipole_order: 3  # Include quadrupoles
decompose_atomwise: true
substrate_idxs:
  - 1
  - 2
  - 3
```

Run: `pyef -c config_complete.yaml`

#### Tips for CLI Usage

**General:**
- Use `#` to add comments in both YAML and CSV files
- Use absolute paths or paths relative to your working directory
- The CLI is designed for batch processing of multiple structures
- Invalid lines in the CSV file will be skipped with a warning

**Performance:**
- ESP calculation: ~15 minutes per 400-atom system (single node)
- E-field calculation: ~1 hour per 400-atom system (single node)
- Electrostatic stabilization: varies with multipole order and atom count

**Dielectric Constants:**
- 1.0 = vacuum
- 2-4 = protein interior
- 20-40 = protein-solvent interface
- 78.5 = water

**Charge Schemes:**
- `Hirshfeld_I`: Most accurate, iterative (recommended)
- `Hirshfeld`: Fast, non-iterative
- `CHELPG`: Fitted to ESP, good for charges
- `Mulliken`: Fast but basis-set dependent

**Multipole Orders:**
- Order 1: Monopole only (charge-charge interactions)
- Order 2: + Dipole (includes charge-dipole, dipole-dipole)
- Order 3: + Quadrupole (includes all terms up to Q×Q)

#### When to Use CLI vs Python API

**Use the CLI when:**
- Processing multiple structures in batch mode
- Integrating PyEF into automated workflows
- Running standardized analysis protocols
- Working in HPC environments with job schedulers
- You want a simple configuration-based approach

**Use the Python API when:**
- Performing interactive analysis and exploration
- Customizing calculations beyond CLI options
- Integrating with other Python analysis tools
- Developing new analysis methods
- You need fine-grained control over parameters

## 4. Installation

### Quick Install (Recommended)
Automated one-command installation with testing:

```bash
git clone git@github.com:davidkastner/pyEF.git
cd pyEF
./install.sh
```

This script will:
- Create a conda environment with all dependencies (including openbabel)
- Install the package in development mode
- Run the test suite to verify installation
- Display activation instructions

### Manual Installation
If you prefer to install manually:

#### Download the package from GitHub
```bash
git clone git@github.com:davidkastner/pyEF.git
```

#### Creating python environment
All the dependencies can be loaded together using the prebuilt environment.yml file.
Compatibility is automatically tested for python versions 3.8 and higher.
Installing all dependencies together has shown to produce more robust installations.

```bash
cd pyEF
conda env create -f environment.yml
conda activate pyef
pip install -e .
```

### Developer install of pyEF
```
cd pyEF
python -m pip install -e .
```

## 5. What's included?
PyEF contains a collection of optimized scripts that work together in pre-built workflows.
If a script is not already included for procedure, many of the functions may be useful in building a procedure.

### File structure

```
.
├── docs                # Readthedocs documentation site
├── pyEF                # Directory containing E3Bind modules
│   ├── cli.py          # Entry point for running E3Bind tools and workflows
│   ├── run             # Runs the workflow
│   ├── analysis        # Analyze the data
│   └── geometry        # Check geometries
└── ...
```


## 6. Documentation
Accurate documentation will always be a high priority for the project.
You can find documentation at the project's [ReadtheDocs](https://pyEF.readthedocs.io/).

### Update the ReadTheDocs

```
make clean
make html
```

### GitHub refresher
#### Push new changes

```
git status
git pull
git add -A .
git commit -m "Change a specific functionality"
git push -u origin main
```

#### Making a pull request
```
git checkout main
git pull

# Before you begin making changes, create a new branch
git checkout -b new-feature-branch
git add -A
git commit -m "Detailed commit message describing the changes"
git push -u origin new-feature-branch

# Visit github.com to add description, submit, merge the pull request

# Once finished on github.com, return to local
git checkout main
git pull

# Delete the remote branch
git branch -d new-feature-branch
```

#### Handle merge conflict

```
git stash push --include-untracked
git stash drop
git pull
```

#### Acknowledgements
Authors: Melissa Manetsch and David W. Kastner

[MolSSi template](https://github.com/molssi/cookiecutter-cms) version 1.6.
