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

**NEW FORMAT (Required):** Explicitly specify paths to .molden and .xyz files:
```csv
# Format: analysis_type, path_to_molden, path_to_xyz, [atom_indices]

# E-field analysis with specific bonds
ef, /path/to/job1/optim.molden, /path/to/job1/optim.xyz, (25, 26), (25, 27)

# ESP analysis with metal center
esp, /path/to/job2/final.molden, /path/to/job2/final.xyz, 30

# Electrostatic stabilization
estab, /path/to/job3/structure.molden, /path/to/job3/structure.xyz

# Combined analysis (ef+esp)
ef+esp, /path/to/job4/optim.molden, /path/to/job4/optim.xyz, 35

# Relative paths work too
ef, ./calculations/job5/optim.molden, ./calculations/job5/optim.xyz, (10, 11)
```

**Supported analysis types:**
- `ef` - Electric field analysis
- `esp` - Electrostatic potential analysis
- `estab` - Electrostatic stabilization analysis
- `ef+esp` - Combined analysis (any combination with `+`)

**Format Details:**
- Column 1: Analysis type (required)
- Column 2: Absolute or relative path to .molden file (required)
- Column 3: Absolute or relative path to .xyz file (required)
- Column 4+: Atom indices (optional, depends on analysis type):
  - For `ef`: Bond pairs as tuples, e.g., `(25, 26), (25, 27)`
  - For `esp`: Single metal atom index, e.g., `30`
  - For `estab`: Not required (uses `substrate_idxs` from config)

**Step 2: Create a config file** (`config.yaml`)
```yaml
input: jobs.csv
dielectric: 1
# Note: ef/esp/estab settings are optional if you specify analysis types in jobs.csv
multiwfn_module: multiwfn
multiwfn_path: /path/to/multiwfn
charge_types:
  - Hirshfeld_I
```

**Config file notes:**
- If you specify analysis types in `jobs.csv` (new format), you don't need `ef/esp/estab` settings
- If you use the old `jobs.csv` format (no analysis types), set `ef/esp/estab` in config
- The `ef/esp/estab` settings in config apply to ALL jobs when using old format

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
| `include_ptchgs` | Include QM/MM point charges | `true` or `false` |
| `ptchg_file` | Path to MM point charge file | Path to charge file |
| `dielectric_scale` | MM charge scaling factor | 1.0 (default), adjusts screening |

### Python API Quick Start

```python
from pyef.analysis import Electrostatics

# Initialize with explicit file paths
molden_paths = ['/path/to/job1/optim.molden', '/path/to/job2/optim.molden']
xyz_paths = ['/path/to/job1/optim.xyz', '/path/to/job2/optim.xyz']

es = Electrostatics(molden_paths, xyz_paths, dielectric=4.0)

# Calculate E-field (no metal indices needed when using bond indices)
df = es.getEfield('Hirshfeld_I', 'output',
                  '/path/to/multiwfn',
                  input_bond_indices=[(25, 26)])

# Calculate stabilization
estab_df = es.getElectrostatic_stabilization(
    '/path/to/multiwfn',
    substrate_idxs=[1,2,3,4,5], multipole_order=2
)

# For QM/MM calculations with point charges
es.includePtChgs('/path/to/pointcharges.txt')
es.set_dielec_scale(1.0)
# Then run E-field or ESP calculations as above
```

**For detailed examples and advanced options, see the Detailed Tutorials section below.**

---

## 3. Detailed Tutorials

### Using PyEF Programmatically (Python API)
The most common way to use PyEF is by importing it directly in your Python scripts:

```python
from pyef.analysis import Electrostatics

# Initialize the Electrostatics object with complete folder paths
folders = ['structure1/scr', 'structure2/scr']  # Complete paths to data
metal_indices = [25, 30]  # 0-indexed atom positions (only needed for ESP)

es = Electrostatics(folders, lst_of_tmcm_idx=metal_indices, dielectric=4.0)

# Prepare data
es.prepData()
es.fix_allECPmolden()

# Calculate ESP
esp_df = es.getESP(
    charge_types=['Hirshfeld_I'],
    ESPdata_filename='esp_output',
    multiwfn_module='multiwfn',
    multiwfn_path='/path/to/multiwfn',
    use_multipole=True
)

# Calculate E-field
efield_df = es.getEfield(
    charge_types='Hirshfeld_I',
    Efielddata_filename='efield_output',
    multiwfn_module='multiwfn',
    multiwfn_path='/path/to/multiwfn',
    input_bond_indices=[(25, 26), (30, 31)],
    multipole_bool=True
)

# Calculate Electrostatic Stabilization
estab_df = es.getElectrostatic_stabilization(
    multiwfn_path='/path/to/multiwfn',
    multiwfn_module='multiwfn',
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
# QM/MM Point Charges (Optional)
# ==========================================

# Include external point charges from MM region in QM/MM calculations
include_ptchgs: false

# Path to point charge file (required if include_ptchgs is true)
# File format: whitespace-delimited, 2 header lines, then columns: charge x y z
# Example file format:
#   Number of point charges
#   charge (e) x (Angstrom) y (Angstrom) z (Angstrom)
#   -0.834  10.234  5.678  -2.345
#    0.417  11.123  6.789   1.234
ptchg_file: path/to/pointcharges.txt

# Dielectric scaling factor for MM charges (default: 1.0)
# MM charges are scaled by 1/sqrt(dielectric_scale)
dielectric_scale: 1.0

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

Create a CSV file (e.g., `jobs.csv`) listing all structures to analyze.

**RECOMMENDED: New Per-Job Format**

Specify analysis type for each job directly in the CSV:

```csv
# ==========================================
# PyEF Job Batch File - Per-Job Analysis Format
# ==========================================
# Format: job_path, analysis_type, [additional parameters]
# Lines starting with '#' are comments
#
# Analysis types: ef, esp, estab, or combinations (ef+esp, esp+estab, etc.)
# Use 0-based indexing for all atom indices

# E-field analysis with specific bonds
structures/complex1.pdb, ef, (25, 26), (25, 27)

# ESP analysis (requires metal index)
structures/complex2.pdb, esp, 30

# Electrostatic stabilization (no indices needed)
structures/complex3.pdb, estab

# Combined analysis
structures/complex4.pdb, ef+esp, 35, 36

# You can mix different analysis types in one batch!
structures/complex5.pdb, ef, (40, 41)
structures/complex6.pdb, esp, 45
```

**Benefits of per-job format:**
- No need to set `ef/esp/estab` in config file
- Different jobs can have different analysis types
- Cleaner, more flexible workflow
- All job info in one place

**OLD FORMAT (Still Supported)**

If you don't specify analysis types in CSV, use config file settings:

```csv
# Simplified format (uses config file analysis settings)
structures/complex1.pdb

# With metal index (for ESP or auto-finding bonds)
structures/complex2.pdb, 25

# With specific bonds (old format)
structures/complex3.pdb, 25, 26

# With tuple bonds (new format)
structures/complex4.pdb, (30, 31), (30, 32)
```

**Column descriptions:**

**New per-job format:**
- **Column 1** (required): Path to the structure file (PDB/molden format)
- **Column 2** (required): Analysis type (`ef`, `esp`, `estab`, or combinations like `ef+esp`)
- **Columns 3+** (varies by analysis):
  - For `ef`: Bond tuples `(atom1, atom2)` or metal_index, bonded_atom_index
  - For `esp`: Metal atom index (0-based)
  - For `estab`: No additional columns needed

**Old format (no analysis type in CSV):**
- **Column 1** (required): Path to the structure file
- **Column 2** (optional): Metal atom index (0-based)
- **Column 3** (optional): Bonded atom index (0-based)
- Or use bond tuples: `(atom1, atom2), (atom3, atom4), ...`

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

`protein_complexes.csv` (with specific bonds):
```csv
# Active site analysis - specify exact bonds (old format)
proteins/1abc_active_site.pdb, 150, 151
proteins/2def_active_site.pdb, 145, 146

# Or use new tuple format for multiple bonds
proteins/3ghi_active_site.pdb, (200, 201), (200, 202), (200, 203)
```

Or `protein_complexes.csv` (simplified - auto-detect bonds):
```csv
# Simplified format - let PyEF auto-detect bonds
proteins/1abc_active_site.pdb
proteins/2def_active_site.pdb
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

`substrate_environment.csv`:
```csv
# Simple format - only job paths needed for estab analysis
structures/system1.pdb
structures/system2.pdb
structures/system3.pdb
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

**Example 4: QM/MM Analysis with Point Charges**

`config_qmmm.yaml`:
```yaml
input: qmmm_systems.csv
dielectric: 4.0
ef: true
esp: true
estab: false
multiwfn_module: multiwfn
multiwfn_path: /usr/local/bin/multiwfn
atmrad_path: /usr/local/share/atmrad
charge_types:
  - Hirshfeld_I
multipole: true
use_multipole: true
# QM/MM point charges from MM region
include_ptchgs: true
ptchg_file: mm_region_charges.txt
dielectric_scale: 1.0
```

`mm_region_charges.txt` (point charge file format):
```
1500
charge(e) x(Ang) y(Ang) z(Ang)
-0.8340  10.234  5.678  -2.345
 0.4170  11.123  6.789   1.234
 0.4170  12.456  7.890   0.123
-0.8340  13.789  8.901  -1.456
...
```

`qmmm_systems.csv`:
```csv
# QM/MM systems - QM region structures
qm_active_site1.pdb, 150, 151
qm_active_site2.pdb, 145, 146
```

Run: `pyef -c config_qmmm.yaml`

**Notes on QM/MM calculations:**
- Point charges represent the MM region atoms in QM/MM calculations
- File format: 2 header lines, then whitespace-delimited columns (charge, x, y, z)
- Coordinates should be in Angstroms, charges in elementary charge units
- The `dielectric_scale` parameter scales MM charges by 1/sqrt(dielectric_scale)
- Point charges are included in both ESP and E-field calculations
- Compatible with Amber, CHARMM, and other MD force field point charges

#### Index Requirements by Analysis Type

PyEF has been designed to minimize required input complexity. **Only ESP analysis requires atom indices!**

| Analysis Type | Required CSV Columns | Example |
|--------------|---------------------|---------|
| **Electrostatic Stabilization** | Column 1: job_path only | `structure.pdb` |
| **Electric Field (auto-detect)** | Column 1: job_path only | `structure.pdb` |
| **Electric Field (specific bonds)** | Columns 1-3: job_path, metal_idx, bond_idx | `structure.pdb, 25, 26` |
| **ESP** | Columns 1-2: job_path, metal_idx | `structure.pdb, 25` |

**Key Points:**
- **ESP analysis**: Requires metal atom index (column 2) - this is validated at runtime
- **E-field analysis**: Can work with just job paths (auto-detect bonds) OR with specific bond indices
- **Estab analysis**: Only needs job paths - substrate atoms specified in config file

#### Tips for CLI Usage

**General:**
- Use `#` to add comments in both YAML and CSV files
- Use absolute paths or paths relative to your working directory
- The CLI is designed for batch processing of multiple structures
- Invalid lines in the CSV file will be skipped with a warning
- **Simplified CSV format**: For most analyses, you only need structure paths!

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

**QM/MM Point Charges:**
- Use for hybrid QM/MM calculations where QM region is analyzed with MM environment
- Point charge file format: 2 header lines, then `charge x y z` (whitespace-delimited)
- First header: number of charges; Second header: column labels
- Coordinates in Angstroms, charges in elementary charge units (e)
- Common sources: Amber prmtop/inpcrd, CHARMM PSF/CRD, GROMACS topology
- The `dielectric_scale` parameter adjusts screening: charges scaled by 1/√(dielectric_scale)
- Point charges contribute to both ESP and E-field calculations
- Use `include_ptchgs: true` in config and specify `ptchg_file` path

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
