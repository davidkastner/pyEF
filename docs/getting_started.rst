Getting Started
===============

*Welcome to PyEF!*

Overview
--------
PyEF is a Python package for analyzing electric fields and electrostatics in molecular systems.
The package is optimized to run using molden files from quantum mechanical calculations and provides
both a command-line interface (CLI) for batch processing and a Python API for interactive analysis.

**Key Features:**

- Electric field analysis on bonds and molecular sites
- Electrostatic potential (ESP) calculations
- Electrostatic stabilization energy analysis
- Support for multipole expansions (monopole, dipole, quadrupole)
- Multiple charge partitioning schemes (Hirshfeld, Hirshfeld_I, CHELPG, etc.)
- QM/MM support with external point charges
- Batch processing via CLI
- Flexible Python API for custom workflows

Installation
------------

Quick Install (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Automated one-command installation with testing:

.. code-block:: bash

    git clone git@github.com:davidkastner/pyEF.git
    cd pyEF
    ./install.sh

This script will:

- Create a conda environment with all dependencies (including openbabel)
- Install the package in development mode
- Run the test suite to verify installation
- Display activation instructions

Manual Installation
~~~~~~~~~~~~~~~~~~~
If you prefer to install manually:

**Step 1: Clone the repository**

.. code-block:: bash

    git clone git@github.com:davidkastner/pyEF.git
    cd pyEF

**Step 2: Create conda environment**

All dependencies can be loaded using the prebuilt environment.yml file.
Compatibility is automatically tested for Python versions 3.8 and higher.

.. code-block:: bash

    conda env create -f environment.yml
    conda activate pyef

**Step 3: Install PyEF**

.. code-block:: bash

    pip install -e .

The ``-e`` flag installs in development mode, allowing you to modify the code.

Quick Start Guide
-----------------

Running Your First Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Step 1: Create a job list (jobs.csv)**

PyEF supports two formats for the job list. The recommended format explicitly specifies analysis types:

.. code-block:: text

    # Format: analysis_type, path_to_molden, path_to_xyz, [atom_indices]

    # E-field analysis with specific bonds
    ef, /path/to/job1/optim.molden, /path/to/job1/optim.xyz, (25, 26), (25, 27)

    # ESP analysis with metal center
    esp, /path/to/job2/final.molden, /path/to/job2/final.xyz, 30

    # Electrostatic stabilization
    estab, /path/to/job3/structure.molden, /path/to/job3/structure.xyz

    # Combined analysis
    ef+esp, /path/to/job4/optim.molden, /path/to/job4/optim.xyz, 35

**Supported analysis types:**

- ``ef`` - Electric field analysis
- ``esp`` - Electrostatic potential analysis
- ``estab`` - Electrostatic stabilization analysis
- ``ef+esp`` - Combined analysis (any combination with +)

**Step 2: Create a config file (config.yaml)**

.. code-block:: yaml

    input: jobs.csv
    dielectric: 1
    multiwfn_module: multiwfn
    multiwfn_path: /path/to/multiwfn
    charge_types:
      - Hirshfeld_I

**Step 3: Run PyEF**

.. code-block:: bash

    pyef -c config.yaml

**Output:** Results saved to ``jobs_Efielddata.csv``

Common Analysis Types
~~~~~~~~~~~~~~~~~~~~~

**Electric Field on Metal-Ligand Bonds**

.. code-block:: yaml

    ef: true
    esp: false
    estab: false
    charge_types:
      - Hirshfeld_I
    multipole: true

**Electrostatic Stabilization Energy**

.. code-block:: yaml

    ef: false
    esp: false
    estab: true
    substrate_idxs: [1, 2, 3, 4, 5]  # Your substrate atoms (0-indexed)
    multipole_order: 2

**Complete Electrostatic Analysis**

.. code-block:: yaml

    ef: true
    esp: true
    estab: true
    substrate_idxs: [1, 2, 3, 4, 5]
    multipole_order: 2  # 1=monopole, 2=+dipole, 3=+quadrupole

Key Parameters
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Parameter
     - Description
     - Common Values
   * - ``dielectric``
     - Medium dielectric constant
     - 1 (vacuum), 4 (protein), 78.5 (water)
   * - ``charge_types``
     - Charge partitioning scheme
     - Hirshfeld_I (recommended), CHELPG, Hirshfeld
   * - ``multipole_order``
     - Expansion order
     - 1 (charges), 2 (+dipoles), 3 (+quadrupoles)
   * - ``decompose_atomwise``
     - Per-atom contributions
     - true or false
   * - ``include_ptchgs``
     - Include QM/MM point charges
     - true or false

Using the Command-Line Interface
---------------------------------

The CLI is designed for batch processing multiple structures. It provides a configuration-based workflow
that's ideal for high-throughput analysis and integration with job schedulers.

Creating a Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a YAML configuration file with your analysis parameters:

.. code-block:: yaml

    # ==========================================
    # PyEF Configuration File
    # ==========================================

    # Input CSV file containing job information (required)
    input: jobs.csv

    # Dielectric constant (default: 1 for vacuum)
    dielectric: 1

    # Analysis Types (set to true to enable)
    ef: true
    esp: false
    estab: false

    # Multiwfn Configuration (required)
    multiwfn_module: multiwfn
    multiwfn_path: /path/to/multiwfn

    # Charge partitioning schemes
    charge_types:
      - Hirshfeld_I

    # Multipole expansion settings
    multipole: true
    multipole_order: 2

    # Atom-wise decomposition
    decompose_atomwise: false

CLI Usage Examples
~~~~~~~~~~~~~~~~~~

**Basic E-field Analysis:**

.. code-block:: bash

    pyef -c config_efield.yaml

**Complete Analysis with All Features:**

.. code-block:: bash

    pyef --config config_complete.yaml

**QM/MM Analysis with Point Charges:**

Create a config file with point charge settings:

.. code-block:: yaml

    include_ptchgs: true
    ptchg_file: mm_region_charges.txt
    dielectric_scale: 1.0

Point charge file format:

.. code-block:: text

    1500
    charge(e) x(Ang) y(Ang) z(Ang)
    -0.8340  10.234  5.678  -2.345
     0.4170  11.123  6.789   1.234
     0.4170  12.456  7.890   0.123
    ...

Using the Python API
--------------------

The Python API provides fine-grained control and is ideal for interactive analysis and custom workflows.

Basic Usage
~~~~~~~~~~~

.. code-block:: python

    from pyef.analysis import Electrostatics

    # Initialize with explicit file paths
    molden_paths = ['/path/to/job1/optim.molden', '/path/to/job2/optim.molden']
    xyz_paths = ['/path/to/job1/optim.xyz', '/path/to/job2/optim.xyz']

    es = Electrostatics(molden_paths, xyz_paths, dielectric=4.0)

    # Calculate E-field
    df = es.getEfield('Hirshfeld_I', 'output',
                      '/path/to/multiwfn',
                      input_bond_indices=[(25, 26)])

    # Calculate stabilization
    estab_df = es.getElectrostatic_stabilization(
        '/path/to/multiwfn',
        substrate_idxs=[1,2,3,4,5],
        multipole_order=2
    )

Advanced Features
~~~~~~~~~~~~~~~~~

**Including QM/MM Point Charges:**

.. code-block:: python

    # For QM/MM calculations with point charges
    es.includePtChgs('/path/to/pointcharges.txt')
    es.set_dielec_scale(1.0)
    # Then run E-field or ESP calculations as above

**ESP Analysis:**

.. code-block:: python

    esp_df = es.getESP(
        charge_types=['Hirshfeld_I'],
        ESPdata_filename='esp_output',
        multiwfn_module='multiwfn',
        multiwfn_path='/path/to/multiwfn',
        use_multipole=True
    )

**Multipole Expansion:**

.. code-block:: python

    # E-field with multipole expansion
    efield_df = es.getEfield(
        charge_types='Hirshfeld_I',
        Efielddata_filename='efield_output',
        multiwfn_module='multiwfn',
        multiwfn_path='/path/to/multiwfn',
        input_bond_indices=[(25, 26), (30, 31)],
        multipole_bool=True
    )

Package Structure
-----------------

.. code-block:: text

    pyEF/
    ├── docs/               # ReadtheDocs documentation
    ├── pyef/               # Main package directory
    │   ├── cli.py          # Command-line interface
    │   ├── run.py          # Workflow execution
    │   ├── analysis.py     # Analysis modules (Electrostatics class)
    │   ├── geometry.py     # Geometry checking utilities
    │   ├── utility.py      # Utility functions
    │   ├── manage.py       # File management
    │   ├── multiwfn_interface.py  # Multiwfn integration
    │   ├── constants.py    # Physical constants
    │   └── tests/          # Test suite
    └── ...

Tips and Best Practices
------------------------

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

- ``Hirshfeld_I``: Most accurate, iterative (recommended)
- ``Hirshfeld``: Fast, non-iterative
- ``CHELPG``: Fitted to ESP, good for charges
- ``Mulliken``: Fast but basis-set dependent

**Multipole Orders:**

- Order 1: Monopole only (charge-charge interactions)
- Order 2: + Dipole (includes charge-dipole, dipole-dipole)
- Order 3: + Quadrupole (includes all terms up to Q×Q)

When to Use CLI vs Python API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Use the CLI when:**

- Processing multiple structures in batch mode
- Integrating PyEF into automated workflows
- Running standardized analysis protocols
- Working in HPC environments with job schedulers

**Use the Python API when:**

- Performing interactive analysis and exploration
- Customizing calculations beyond CLI options
- Integrating with other Python analysis tools
- Developing new analysis methods

Developer Guide
---------------

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

To update the ReadtheDocs site locally:

.. code-block:: bash

    cd docs
    make clean
    make html

The compiled docs will be in the ``_build/html`` directory.

GitHub Workflow
~~~~~~~~~~~~~~~

**Push new changes:**

.. code-block:: bash

    git status
    git pull
    git add -A .
    git commit -m "Descriptive commit message"
    git push -u origin main

**Making a pull request:**

.. code-block:: bash

    git checkout main
    git pull

    # Create a new branch
    git checkout -b new-feature-branch
    git add -A
    git commit -m "Detailed commit message"
    git push -u origin new-feature-branch

    # Visit github.com to submit and merge the pull request

    # After merging, return to main
    git checkout main
    git pull
    git branch -d new-feature-branch

**Handle merge conflict:**

.. code-block:: bash

    git stash push --include-untracked
    git stash drop
    git pull

Copyright
---------

Copyright (c) 2025, Melissa Manetsch and David W. Kastner