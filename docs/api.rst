API Documentation
=================

This page contains the API reference for PyEF. The package is organized into several modules,
each providing specific functionality for electric field and electrostatics analysis.

Main Modules
------------

.. autosummary::
   :toctree: autosummary/
   :recursive:

   pyef.analysis
   pyef.cli
   pyef.run
   pyef.geometry
   pyef.utility
   pyef.manage
   pyef.multiwfn_interface
   pyef.constants

Module Descriptions
-------------------

**pyef.analysis**
    Core analysis module containing the ``Electrostatics`` class for calculating electric fields,
    electrostatic potentials, and electrostatic stabilization energies. This is the main module
    for performing calculations.

**pyef.cli**
    Command-line interface module that provides the ``pyef`` command for batch processing.
    Handles argument parsing and workflow orchestration.

**pyef.run**
    Workflow execution module that manages the processing pipeline for multiple jobs,
    including data preparation, analysis execution, and results compilation.

**pyef.geometry**
    Geometry checking and validation utilities. Provides functions for verifying molecular
    geometries and detecting structural issues.

**pyef.utility**
    General utility functions for file I/O, data manipulation, coordinate transformations,
    and other common operations used throughout the package.

**pyef.manage**
    File management utilities for handling molden files, XYZ files, and other molecular
    structure formats. Includes functions for data organization and preprocessing.

**pyef.multiwfn_interface**
    Interface module for communicating with the Multiwfn program. Handles input generation,
    process execution, and output parsing for charge analysis and multipole calculations.

**pyef.constants**
    Physical and mathematical constants used in electrostatics calculations, including
    conversion factors and fundamental constants.

Key Classes and Functions
--------------------------

Electrostatics Class
~~~~~~~~~~~~~~~~~~~~

The ``Electrostatics`` class in ``pyef.analysis`` is the primary interface for all calculations:

.. code-block:: python

    from pyef.analysis import Electrostatics

    # Initialize
    es = Electrostatics(molden_paths, xyz_paths, dielectric=4.0)

    # Main methods:
    # - getEfield(): Calculate electric fields
    # - getESP(): Calculate electrostatic potentials
    # - getElectrostatic_stabilization(): Calculate stabilization energies
    # - includePtChgs(): Include QM/MM point charges
    # - set_dielec_scale(): Set dielectric scaling factor

For detailed method signatures and parameters, see the auto-generated API documentation below.

Complete API Reference
----------------------

.. automodule:: pyef.analysis
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pyef.cli
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pyef.run
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pyef.geometry
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pyef.utility
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pyef.manage
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pyef.multiwfn_interface
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: pyef.constants
   :members:
   :undoc-members:
   :show-inheritance:
