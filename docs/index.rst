.. PyEF documentation master file

.. image:: _static/logo-white.svg

PyEF: Electric Field Analysis for Molecular Systems
====================================================

.. container:: .large

   PyEF is a Python package for analyzing electric fields, electrostatic potentials,
   and electrostatic stabilization in molecular systems from quantum mechanical calculations.


.. container:: .buttons

   `Getting Started <getting_started.html>`_
   `API Reference <api.html>`_
   `GitHub <https://github.com/davidkastner/pyef>`_

Overview
--------

PyEF provides a comprehensive toolkit for calculating and analyzing electrostatic interactions
in molecular systems. It is optimized for processing molden files from quantum mechanical calculations
and offers both a command-line interface for high-throughput batch processing and a Python API
for interactive analysis.

Key Features
~~~~~~~~~~~~

**Analysis Capabilities:**

- **Electric Field Analysis**: Calculate electric fields at specific bonds and molecular sites
- **Electrostatic Potential (ESP)**: Compute ESP at metal centers and other points of interest
- **Electrostatic Stabilization**: Quantify stabilization energies between molecular fragments
- **Multipole Expansions**: Support for monopole, dipole, and quadrupole contributions
- **QM/MM Integration**: Include external point charges from MM regions

**Technical Features:**

- Multiple charge partitioning schemes (Hirshfeld, Hirshfeld_I, CHELPG, Mulliken, etc.)
- Atom-wise decomposition of electrostatic contributions
- Dielectric screening effects
- Automated geometry validation
- Integration with Multiwfn for charge analysis

**Interfaces:**

- **CLI**: Batch processing with YAML configuration files
- **Python API**: Fine-grained control for interactive analysis and custom workflows
- **HPC-Ready**: Designed for use in high-performance computing environments

Quick Start
-----------

Installation
~~~~~~~~~~~~

.. code-block:: bash

   git clone git@github.com:davidkastner/pyEF.git
   cd pyEF
   ./install.sh

Basic Usage (CLI)
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Create a job list (jobs.csv)
   echo "ef, path/to/optim.molden, path/to/optim.xyz, (25, 26)" > jobs.csv

   # Create a config file (config.yaml)
   cat > config.yaml << EOF
   input: jobs.csv
   dielectric: 1
   multiwfn_path: /path/to/multiwfn
   charge_types:
     - Hirshfeld_I
   EOF

   # Run analysis
   pyef -c config.yaml

Basic Usage (Python API)
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pyef.analysis import Electrostatics

   # Initialize
   es = Electrostatics(['optim.molden'], ['optim.xyz'], dielectric=4.0)

   # Calculate E-field
   df = es.getEfield('Hirshfeld_I', 'output', '/path/to/multiwfn',
                     input_bond_indices=[(25, 26)])

   # Calculate stabilization
   estab_df = es.getElectrostatic_stabilization('/path/to/multiwfn',
                                                 substrate_idxs=[1,2,3,4,5],
                                                 multipole_order=2)

Use Cases
---------

PyEF is designed for researchers studying:

- **Enzyme catalysis**: Analyzing electric fields in active sites
- **Metalloproteins**: Electrostatic effects at metal centers
- **Organometallic chemistry**: Understanding ligand-metal interactions
- **Protein-ligand binding**: Quantifying electrostatic contributions to binding
- **Reaction mechanisms**: Electrostatic stabilization of transition states
- **Environmental effects**: Modeling solvent and protein dielectric screening

Applications
~~~~~~~~~~~~

**Electric Field Analysis:**
   Measure electric fields along specific bonds to understand how the environment
   influences molecular properties and reactivity.

**Electrostatic Stabilization:**
   Quantify how much an enzyme or protein environment electrostatically stabilizes
   a substrate, intermediate, or transition state.

**QM/MM Calculations:**
   Incorporate point charges from molecular mechanics regions to study electrostatics
   in large biomolecular systems.

Performance
-----------

Typical calculation times (single node):

- **ESP**: ~15 minutes per 400-atom system
- **E-field**: ~1 hour per 400-atom system
- **Electrostatic stabilization**: Varies with multipole order and system size

The package is designed for efficient batch processing on HPC clusters.

Documentation Contents
----------------------

.. toctree::
   :caption: User Guide
   :maxdepth: 2

   getting_started
   api

.. toctree::
   :caption: External Links
   :maxdepth: 1

   GitHub Repository <https://github.com/davidkastner/pyef>
   Author Website <https://kastner.io/>

Citation
--------

If you use PyEF in your research, please cite:

.. code-block:: bibtex

   @software{pyef,
     title = {PyEF: Electric Field Analysis for Molecular Systems},
     author = {Manetsch, Melissa and Kastner, David W.},
     year = {2025},
     url = {https://github.com/davidkastner/pyef}
   }

Support and Contributing
------------------------

- **Issues**: Report bugs and request features on `GitHub Issues <https://github.com/davidkastner/pyef/issues>`_
- **Discussions**: Ask questions and share ideas on `GitHub Discussions <https://github.com/davidkastner/pyef/discussions>`_
- **Contributing**: See our `Contributing Guidelines <https://github.com/davidkastner/pyef/blob/main/CONTRIBUTING.md>`_

License
-------

PyEF is released under the MIT License. See the LICENSE file for details.

Acknowledgments
---------------

**Authors:** Melissa Manetsch and David W. Kastner

Built with the `MolSSI Cookiecutter <https://github.com/molssi/cookiecutter-cms>`_ template.

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`