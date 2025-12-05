"""Electrostatics Analysis Module for Quantum Chemistry Calculations.

This module provides the Electrostatics class for computing electrostatic
properties (ESP, electric fields, partial charges) from quantum chemistry
calculations. It interfaces with Multiwfn for charge analysis and supports
various charge partitioning schemes (Hirshfeld, Hirshfeld-I, CHELPG, etc.).

Key Features:
    - Multiple charge partitioning schemes
    - ESP and electric field calculations
    - Multipole moment analysis
    - QMMM point charge integration
    - Dielectric environment effects
    - ECP (Effective Core Potential) support

Dependencies:
    - Multiwfn (external program)
    - NumPy, pandas, scipy
    - OpenBabel, BioPandas
"""

import re
import os
import mmap
import glob
import shutil
import logging
import traceback
import subprocess
import math
import time

import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy.optimize import LinearConstraint, minimize
from collections import deque
from importlib import resources
from distutils.dir_util import copy_tree

#import openbabel
#from biopandas.pdb import PandasPdb

from .geometry import Geometry, Visualize
from .utility import MoldenObject
from .multiwfn_interface import MultiwfnInterface
from . import constants
class Electrostatics:
    """Compute electrostatic properties from quantum chemistry calculations.

    This class processes series of quantum chemistry calculations to extract
    electrostatic properties including ESP, electric fields, and partial charges
    using various charge partitioning schemes via Multiwfn.

    Attributes
    ----------
    lst_of_folders : list of str
        Names of folders containing computation outputs.
    lst_of_tmcm_idx : list of int
        Atom indices where ESP will be computed (typically metal centers).
    folder_to_file_path : str
        Path from folder location to .molden file.
    config : dict
        Configuration dictionary containing:
            - hasECP (bool): Whether ECP was used in calculation
            - includePtChgs (bool): Include point charges in ESP calculation
            - ptChgfp (str): Path to point charge file
            - molden_filename (str): Name of molden file
            - xyzfilename (str): Name of xyz coordinate file
            - rerun (bool): Force recalculation of charges
            - maxIHirshBasis (int): Max basis functions for Hirshfeld-I
            - maxIHirshFuzzyBasis (int): Max basis for fuzzy Hirshfeld-I
            - ECP (str): ECP basis set family
            - dielectric (float): Dielectric constant
            - dielectric_scale (float): Scaling factor for dielectric
            - chgprefix (str): Prefix for charge filenames
            - excludeAtomfromEcalc (list): Atoms to exclude from E-field calc
            - changeDielectBoundBool (bool): Use dielectric=1 for bonded atoms
            - visualize_efield (bool): Create PDB files for atom-wise E-fields
            - visualize_charges (bool): Create PDB files for partial charges
            - visualize_per_bond (bool): Create separate PDB files for each bond
    dict_of_calcs : dict
        Mapping of charge scheme names to Multiwfn command codes.
    dict_of_multipole : dict
        Mapping of multipole schemes to Multiwfn command sequences.
    ptchgdf : pd.DataFrame or None
        DataFrame containing point charge data for QMMM calculations.
    periodic_table : dict
        Mapping of element symbols to atomic numbers.
    amassdict : dict
        Atomic properties (mass, atomic number, covalent radius, valence).

    Examples
    --------
    >>> # Recommended: Use complete folder paths
    >>> folders = ['calc1/scr', 'calc2/scr']
    >>> es = Electrostatics(folders, dielectric=2.0)
    >>>
    >>> # Legacy: Use separate folder names and subdirectory
    >>> folders = ['calc1', 'calc2']
    >>> es = Electrostatics(folders, '/scr/', dielectric=2.0)
    >>>
    >>> # For ESP calculations (metal indices required)
    >>> es_with_metals = Electrostatics(folders, lst_of_tmcm_idx=[0, 0], dielectric=2.0)
    >>> data = es_with_metals.getESP(['Hirshfeld', 'Hirshfeld_I'], 'esp_results', ...)
    """

    def __init__(self, lst_of_folders, folder_to_file_path='', lst_of_tmcm_idx=None,
                 **kwargs):
        """Initialize Electrostatics analysis object.

        Parameters
        ----------
        lst_of_folders : list of str
            List of folder paths containing calculation outputs.
            Can be complete paths (e.g., ['job1/scr', 'job2/scr']) or just
            folder names if using folder_to_file_path.
        folder_to_file_path : str, optional
            Optional subdirectory path to append to each folder (default: '').
            Use '' or '/' if lst_of_folders contains complete paths.
            Legacy usage: '/scr/' to append 'scr/' to each folder name.
        lst_of_tmcm_idx : list of int, optional
            List of atom indices for ESP calculation (0-indexed).
            Only required when running ESP analysis. Not needed for E-field
            or electrostatic stabilization if bond indices are provided.
        **kwargs : dict, optional
            Configuration options to override defaults. See class docstring
            for available configuration keys.

        Notes
        -----
        ECP options: "stuttgart_rsc", "def2", "crenbl", "lanl2dz", "lacvps"
        """
        # Initialize configuration with defaults
        self.config = {
            'hasECP': False,              # Use Effective Core Potentials
            'includePtChgs': False,       # Include QMMM point charges
            'ptChgfp': '',                # Point charge file path
            'molden_filename': 'final_optim.molden',  # Molden file name
            'xyzfilename': 'final_optim.xyz',         # XYZ file name
            'rerun': False,               # Force recalculation
            'maxIHirshBasis': 12000,      # Hirshfeld-I basis limit
            'maxIHirshFuzzyBasis': 6000,  # Fuzzy Hirshfeld-I basis limit
            'ECP': "lacvps",              # ECP basis set family
            'dielectric': 1,              # Dielectric constant
            'dielectric_scale': 1,        # Dielectric scaling factor
            'chgprefix': '',              # Charge filename prefix
            'excludeAtomfromEcalc': [],   # Atoms to exclude from E-field
            'changeDielectBoundBool': False,  # Special dielectric for bonds
            'visualize_efield': False,    # Create PDB files for atom-wise E-fields
            'visualize_charges': False,   # Create PDB files for partial charges
            'visualize_per_bond': False   # Create separate PDB files for each bond
        }
        # Override defaults with user-provided kwargs
        self.config.update(kwargs)

        # Store required parameters
        self.lst_of_folders = lst_of_folders
        self.lst_of_tmcm_idx = lst_of_tmcm_idx if lst_of_tmcm_idx is not None else []
        self.folder_to_file_path = folder_to_file_path

        # Initialize Multiwfn interface with current config
        self.multiwfn = MultiwfnInterface(config=self.config)

        # Periodic table: element symbol -> atomic number
        # Source: molsimplify, http://www.webelements.com/
        # Note: Full capitalization included for compatibility with various
        # software that write .xyz files in all caps
        # Use centralized constants from constants module
        self.periodic_table = constants.PERIODIC_TABLE
        self.amassdict = constants.ATOMIC_MASS_DICT

        # Prepare data files (extract final frames, check file existence)
        self.prepData()

        # Fix molden files if ECP was used in calculations
        if self.config['hasECP']:
            self.fix_allECPmolden()

    # =========================================================================
    # Multiwfn Interface Methods
    # =========================================================================

    def _run_multiwfn(self, command, input_commands, output_file=None,
                     description="Multiwfn calculation", capture_output=True):
        """Centralized Multiwfn runner with proper error handling and output display.

        Parameters
        ----------
        command : str
            The shell command to run Multiwfn (e.g., "multiwfn file.molden")
        input_commands : list of str
            List of commands to send to Multiwfn's stdin
        output_file : str, optional
            Path to redirect stdout to a file
        description : str, optional
            Description of the calculation for error messages
        capture_output : bool, optional
            If True, capture and display output. If False, let it stream to terminal

        Returns
        -------
        tuple
            (stdout, stderr, returncode) from the Multiwfn process

        Raises
        ------
        RuntimeError
            If Multiwfn fails with non-zero exit code
        """
        import sys

        # Prepare the full command with output redirection if specified
        full_command = command
        if output_file:
            full_command = f"{command} > {output_file}"

        print(f"\n{'='*60}")
        print(f"Running: {description}")
        print(f"Command: {full_command}")
        print(f"{'='*60}\n")

        try:
            # Run Multiwfn process
            if capture_output:
                proc = subprocess.Popen(
                    full_command,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=True,
                    text=True,
                    bufsize=1
                )

                # Send commands to stdin
                input_str = "\n".join(input_commands)
                stdout, stderr = proc.communicate(input=input_str)

                # Display output in real-time style
                if stdout:
                    print("MULTIWFN OUTPUT:")
                    print("-" * 60)
                    print(stdout)
                    print("-" * 60)

                if stderr:
                    print("MULTIWFN STDERR:")
                    print("-" * 60)
                    print(stderr, file=sys.stderr)
                    print("-" * 60)

            else:
                # Let output stream directly to terminal
                proc = subprocess.Popen(
                    full_command,
                    stdin=subprocess.PIPE,
                    shell=True,
                    text=True
                )
                input_str = "\n".join(input_commands)
                proc.communicate(input=input_str)
                stdout, stderr = "", ""

            returncode = proc.returncode

            # Check for errors
            if returncode != 0:
                error_msg = f"""
{'='*60}
ERROR: {description} FAILED
{'='*60}
Command: {full_command}
Return code: {returncode}

STDOUT:
{stdout if stdout else '(empty)'}

STDERR:
{stderr if stderr else '(empty)'}

Input commands sent:
{chr(10).join(input_commands)}
{'='*60}
"""
                raise RuntimeError(error_msg)

            print(f"\n✓ {description} completed successfully\n")
            return stdout, stderr, returncode

        except Exception as e:
            error_msg = f"""
{'='*60}
EXCEPTION during {description}
{'='*60}
Command: {full_command}
Error: {str(e)}

Traceback:
{traceback.format_exc()}

Input commands that were to be sent:
{chr(10).join(input_commands)}
{'='*60}
"""
            print(error_msg, file=sys.stderr)
            raise RuntimeError(error_msg) from e

    # =========================================================================
    # Configuration Methods
    # =========================================================================

    def updateCalcSettings(self, key, value):
        """Update calculation settings dictionary.

        Parameters
        ----------
        key : str
            Setting key to update.
        value : any
            New value for the setting.
        """
        self.dict_settings[key] = value

    def setRunrunBool(self, rename, setbool=True):
        """Set rerun flag and charge filename prefix.

        Parameters
        ----------
        rename : str
            Prefix for charge output filenames.
        setbool : bool, default=True
            Whether to force recalculation of charges.
        """
        self.config['chgprefix'] = rename
        self.config['rerun'] = setbool

    def setExcludeAtomFromCalc(self, lstExcludeAtoms):
        """Set atoms to exclude from electric field calculations.

        Parameters
        ----------
        lstExcludeAtoms : list of int
            Atom indices (0-indexed) to exclude from E-field calculations.
        """
        self.config['excludeAtomfromEcalc'] = lstExcludeAtoms

    def includePtChgs(self, name_ptch_file):
        """Include point charges in ESP calculation (for QMMM).

        Parameters
        ----------
        name_ptch_file : str
            Path to point charge file.

        Notes
        -----
        Sets includePtChgs config to True and stores file path.
        Point charges are typically from MM region in QM/MM calculations.
        """
        self.config['ptChgfp'] = name_ptch_file
        self.config['includePtChgs'] = True
        print(f'Point charges to be included via {name_ptch_file}')

    def set_dielec_scale(self, dielec):
        """Set dielectric scaling factor for point charges.

        Parameters
        ----------
        dielec : float
            Scaling factor for dielectric in QMMM calculations.
        """
        self.config['dielectric_scale'] = dielec

    def initialize_excludeAtomsFromEfieldCalc(self, atom_to_exclude):
        """Initialize atoms to exclude from electric field calculation.

        Parameters
        ----------
        atom_to_exclude : list of int
            List of atom indices (0-indexed) to exclude.
        """
        self.config['excludeAtomfromEcalc'] = atom_to_exclude

    def minDielecBonds(self, bool_bonds):
        """Set dielectric to 1 for bonded atoms in ESP calculation.

        Parameters
        ----------
        bool_bonds : bool
            If True, use dielectric=1 for atoms bonded to ESP center.
            Helps avoid over-estimating screening from bound atoms.
        """
        self.config['changeDielectBoundBool'] = bool_bonds

    def runlowmemory(self):
        """Enable low-memory mode for large systems (>300 atoms).

        Notes
        -----
        Switches Hirshfeld-I calculation to slower but lower-memory algorithm.
        Useful for large systems that would otherwise exceed memory limits.
        """
        self.multiwfn.dict_of_multipole = {
            'Hirshfeld': ['3', '2'],
            'Hirshfeld_I': ['4', '-2', '1', '2'],
            'Becke': ['1', '2']
        }
        self.multiwfn.dict_of_calcs = {
            'Hirshfeld': '1', 'Voronoi': '2', 'Mulliken': '5', 'Lowdin': '6',
            'SCPA': '7', 'Becke': '10', 'ADCH': '11', 'CHELPG': '12',
            'MK': '13', 'AIM': '14', 'Hirshfeld_I': ['15', '-2'],
            'CM5': '16', 'EEM': '17', 'RESP': '18', 'PEOE': '19'
        }

    def changeDielectric(self, dlc):
        """Change dielectric constant of the environment.

        Parameters
        ----------
        dlc : float
            Dielectric constant of solvent/environment.

        Notes
        -----
        Used to model solvation effects in ESP calculations.
        Default is 1 (vacuum). Common values: water=80, protein=4.
        """
        self.config['dielectric'] = dlc

    def set_molden_filename(self, new_name):
        """Set the expected molden filename.

        Parameters
        ----------
        new_name : str
            Name of the molden file to read.
        """
        self.config['molden_filename'] = new_name
        self.multiwfn.config['molden_filename'] = new_name

    def set_xyzfilename(self, new_name):
        """Set the expected XYZ coordinate filename.

        Parameters
        ----------
        new_name : str
            Name of the XYZ file to read.
        """
        self.config['xyzfilename'] = new_name
        self.multiwfn.config['xyzfilename'] = new_name

    def rePrep(self):
        """Re-run data preparation (file extraction and validation).

        Notes
        -----
        Useful if files have been updated or regenerated.
        """
        self.prepData()

    # =========================================================================
    # File Preparation and Utility Methods
    # =========================================================================

    def fix_allECPmolden(self):
        """Fix molden files that used Effective Core Potentials.

        Notes
        -----
        ECPs can cause artifacts in molden files that make them incompatible
        with Multiwfn. This method reformats the molden files to fix these issues.
        Processes all folders in lst_of_folders.
        """
        owd = os.getcwd()
        basis_set_fam = self.config['ECP']
        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        print('   > Re-formatting .molden files to fix ECP artifacts')

        for f in list_of_folders:
            folder_path = os.path.join(owd, f + folder_to_molden)
            final_optim_molden = os.path.join(
                folder_path, self.config['molden_filename']
            )
            final_optim_xyz = os.path.join(
                folder_path, self.config['xyzfilename']
            )
            molden = MoldenObject(final_optim_xyz, final_optim_molden)
            molden.fix_ECPmolden(owd, basis_set_fam)

    def prepData(self):
        """Prepare calculation data for analysis.

        Extracts final frames from optimization trajectories, validates file
        existence, and sets up consistent naming conventions for .molden and
        .xyz files across all calculation folders.

        Notes
        -----
        - Extracts last frame from optim.xyz to create final_optim.xyz
        - Locates .molden files and updates config if non-standard names found
        - Removes folders from processing list if files cannot be located
        """

        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        owd = os.getcwd()
        print('   > Pre-processing data')
        backup_xyz = ['xyz.xyz']

        for f in list_of_folders:
            try:
                folder_path = os.path.join(owd, f + folder_to_molden)
                print('      > .molden and .xyz file should be located here: ' + folder_path)

                # Processing optim.xyz to create final_optim.xyz
                final_optim_xyz = os.path.join(folder_path, self.config['xyzfilename'])

                # Copying .molden files to final_optim.molden
                final_optim_molden = os.path.join(folder_path, self.config['molden_filename'])

                #this is for the full optimization cycle
                optim_file_path = os.path.join(folder_path, 'optim.xyz')


                if not os.path.exists(final_optim_molden):
                    print(f'Expected .molden with filename: {self.config["molden_filename"]} in {folder_path}. We could not find a .molden filename with the default prefix, you can alter using: set_molden_filename()')
                    print(f"For now searching for all .molden files in directory. If you only have one .molden file in this directory, the defauly prefix will be altered ")
                    files = glob.iglob(os.path.join(folder_path, "*.molden"))
                    for file in files:
                        if os.path.abspath(file) != os.path.abspath(final_optim_molden):
                            #change the defauly path the the molden file!
                            self.set_molden_filename(os.path.basename(file))
                            #set the backup_xyz filename to this prefix here. Generally a good guess!
                            file_prefix, _ = os.path.splitext(os.path.basename(file))
                            backup_xyz = file_prefix + '.xyz'
                            print(f'Default .molden file name is now changed to {self.config["molden_filename"]}')
                            break
                    else:
                        raise Exception("Unable to locate any .molden file in directory: {f}")
                else:
                    print(f"      > {self.config['molden_filename']} sucessflly located in {folder_path}.")



                if not os.path.exists(final_optim_xyz):
                    try:
                        with open(optim_file_path, 'r') as full_traj:
                            num_atoms = int(full_traj.readline())
                            num_lines = num_atoms + 2
                            # If ends on first optimization of cycle, change cutting pattern
                            full_traj.seek(0)
                            head = [next(full_traj) for _ in range(num_lines)]
                            full_traj.seek(0)
                            with open(os.path.join(folder_path, 'initial_' + os.path.basename(optim_file_path)), 'w') as initxyz:
                                initxyz.writelines(head)
                            with open(final_optim_xyz, 'w') as finalxyz:
                                finalxyz.writelines(deque(full_traj, num_lines))
                    except Exception as e:
                        print(f'Expected .xyz with filename: {self.config["xyzfilename"]} in {folder_path}. If your xyz filename does NOT match the default, you can alter using: set_xyzfilename()')
                        print(f'We will try to use the anticipated prefix from the associated molden file: {backup_xyz}')
                        backup_file = os.path.join(folder_path, backup_xyz)
                        if os.path.exists(backup_file):
                            self.config['xyzfilename'] = backup_xyz
                            logging.info(f'Single point data found. Using {backup_xyz} as fallback.')
                        else:
                            raise Exception("Could not locate .molden or corresponding .xyz file in directory {f}")

                else:
                    print(f'      > {self.config["xyzfilename"]} succesfully located in {folder_path}.')            

                # Copying .molden files to final_optim.molden
                final_optim_molden = os.path.join(folder_path, self.config['molden_filename'])
                backup_xyz = 'final_optim'
            except Exception as e:
                print(f'Could not locate molden and corresponding xyz file in {f} please rename these files if they exist')
                print(f'For now I am just taking this folder out of the list')

                list_of_folders.remove(f)
                os.chdir(owd)
        self.lst_of_folders = list_of_folders
        os.chdir(owd)

    # =========================================================================
    # Multipole and Charge Parsing Methods
    # =========================================================================

    @staticmethod
    def getmultipoles(multipole_name):
        """Parse multipole moments from Multiwfn output file.

        Parameters
        ----------
        multipole_name : str
            Path to multipole moment output file from Multiwfn.

        Returns
        -------
        list of dict
            List of dictionaries containing multipole data for each atom:
            - 'Index': atom index (str)
            - 'Element': element symbol (str)
            - 'Atom_Charge': atomic charge (float)
            - 'Dipole_Moment': dipole moment vector [x, y, z] (list of float)
            - 'Quadrupole_Moment': 3x3 quadrupole tensor (np.ndarray)

        Notes
        -----
        Parses Multiwfn output format with sections separated by asterisks.
        Extracts monopole (charge), dipole, and quadrupole moments for each atom.
        """
        # Read the text file
        with open(multipole_name, 'r') as file:
            text = file.read()
            # Split the text into sections based on the delimiter
            sections = re.split(r'\n\s*[*]+\s*', text)

        # Define a pattern to extract the index, element name, and atomic dipole and quadrupole moment values
        #pattern = r'Atomic multipole moments of\s+(\d+)\s*\(([^)]+)\s*\).*?X=\s*([-+]?\d+\.\d+)\s+Y=\s*([-+]?\d+\.\d+)\s+Z=\s*([-+]?\d+\.\d+).*?XX=\s*([-+]?\d+\.\d+)\s+XY=\s*([-+]?\d+\.\d+)\s+XZ=\s*([-+]?\d+\.\d+).*?YX=\s*([-+]?\d+\.\d+)\s+YY=\s*([-+]?\d+\.\d+)\s+YZ=\s*([-+]?\d+\.\d+).*?ZX=\s*([-+]?\d+\.\d+)\s+ZY=\s*([-+]?\d+\.\d+)\s+ZZ=\s*([-+]?\d+\.\d+)'
        pattern = (r'Atomic multipole moments of\s+(\d+)\s*\(([^)]+)\s*\).*?'
    r'Atomic monopole moment \(electron\):\s*([-+]?\d+\.\d+)\s+Atomic charge:\s*([-+]?\d+\.\d+).*?'
    r'X=\s*([-+]?\d+\.\d+)\s+Y=\s*([-+]?\d+\.\d+)\s+Z=\s*([-+]?\d+\.\d+).*?'
    r'XX=\s*([-+]?\d+\.\d+)\s+XY=\s*([-+]?\d+\.\d+)\s+XZ=\s*([-+]?\d+\.\d+).*?'
    r'YX=\s*([-+]?\d+\.\d+)\s+YY=\s*([-+]?\d+\.\d+)\s+YZ=\s*([-+]?\d+\.\d+).*?'
    r'ZX=\s*([-+]?\d+\.\d+)\s+ZY=\s*([-+]?\d+\.\d+)\s+ZZ=\s*([-+]?\d+\.\d+)'
    )
        atomicDicts = []

        # Iterate through sections and extract the information
        for section in sections:
            match = re.search(pattern, section, re.DOTALL)
            if match:
                index = match.group(1)
                element = match.group(2)
                monopole_moment = float(match.group(3))
                atomic_charge = float(match.group(4))
                x_value = match.group(5)
                y_value = match.group(6)
                z_value = match.group(7)
                dipole_moment = [float(x_value), float(y_value), float(z_value)]

                xx_value = match.group(8)
                xy_value = match.group(9)
                xz_value = match.group(10)
                yx_value = match.group(11)
                yy_value = match.group(12)
                yz_value = match.group(13)
                zx_value = match.group(14)
                zy_value = match.group(15)
                zz_value = match.group(16)
                # Create a 3x3 matrix for the quadrupole moment
                quadrupole_moment = np.array([
                [float(xx_value), float(xy_value), float(xz_value)],
                [float(yx_value), float(yy_value), float(yz_value)],
                [float(zx_value), float(zy_value), float(zz_value)]
                ])

                atomDict = {"Index": index, "Element": element, "Atom_Charge": atomic_charge, 'Dipole_Moment': dipole_moment, 'Quadrupole_Moment': quadrupole_moment}
                atomicDicts.append(atomDict)
        return atomicDicts

    def getPtChgs(self, filename_pt):
        """Load point charges from file (for QMMM calculations).

        Parameters
        ----------
        filename_pt : str
            Path to point charge file (typically from Amber/CHARMM).

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: ['Atom', 'charge', 'x', 'y', 'z']
            All atoms labeled as 'pnt' (point charge).

        Notes
        -----
        File format: whitespace-delimited with 2 header lines to skip.
        """
        chg_df = pd.read_table(
            filename_pt, skiprows=2, delim_whitespace=True,
            names=['charge', 'x', 'y', 'z']
        )
        atm_name = ['pnt']
        atoms = atm_name * len(chg_df['charge'])
        chg_df['Atom'] = atoms
        return chg_df

    @staticmethod
    def mapcount(filename):
        """Rapidly count lines in a file using memory mapping.

        Parameters
        ----------
        filename : str
            Path to file to count.

        Returns
        -------
        int
            Number of lines in the file.

        Notes
        -----
        Uses mmap for efficient line counting in large files.
        """
        f = open(filename, "r+")
        buf = mmap.mmap(f.fileno(), 0)
        lines = 0
        readline = buf.readline

        while readline():
            lines += 1
        return lines

    # =========================================================================
    # Electrostatic Potential (ESP) Calculation Methods
    # =========================================================================

    def calcesp(self, path_to_xyz, espatom_idx, charge_range, charge_file):
        """Calculate electrostatic potential at specified atom.

        Parameters
        ----------
        path_to_xyz : str
            Path to XYZ coordinate file.
        espatom_idx : int
            Atom index (0-indexed) where ESP will be calculated.
        charge_range : list of int
            Indices of atoms to include in ESP calculation.
        charge_file : str
            Path to partial charge file from Multiwfn.

        Returns
        -------
        list
            [ESP (float), atom_symbol (str)]
            ESP in Volts at the specified atom.

        Notes
        -----
        - Uses Coulomb's law: V = k * q / r with dielectric screening
        - Point charges from QMMM included if config['includePtChgs'] is True
        - Bonded atoms treated with dielectric=1 if config['changeDielectBoundBool'] is True
        - Units: Volts (N·m/C = J/C)
        """
        dielectric = self.config['dielectric']
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Convert each column to list for quicker indexing
        atoms = list(df['Atom'])
        charges = list(df['charge'])
        xs = list(df['x'])
        ys = list(df['y'])
        zs = list(df['z'])

        #For QMMM calculation, include point charges in ESP calculation
        if self.config['includePtChgs']:
            df_ptchg  = self.getPtChgs(self.config['ptChgfp'])
            xs = xs + list(df_ptchg['x'])
            ys = ys + list(df_ptchg['y'])
            zs = zs + list(df_ptchg['z'])
            charges = charges + list(df_ptchg['charge'])
            atoms = atoms + list(df_ptchg['Atom'])
            #change charge range to include these new partial charges!
            charge_range = range(0, len(xs))
        # Pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        # Determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        total_esp = 0

        bound_atoms = []
        #create list of bound atoms, these are treated with a different dielectric
        if self.config['changeDielectBoundBool']:
            bound_atoms = self.getBondedAtoms(path_to_xyz, idx_atom)
 
        for idx in charge_range:
            if idx == idx_atom:
                continue
            elif idx in bound_atoms:
                #now account for bound atoms
                 r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                 total_esp = total_esp + (charges[idx]/r)
            else:
                # Calculate esp and convert to units (A to m)
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                total_esp = total_esp + (1/dielectric)*(charges[idx]/r)

        final_esp = constants.COULOMB_CONSTANT*total_esp*constants.ELEMENTARY_CHARGE  #N*m^2/(C^2)*(C/m) = N*m/C = J/C = Volt
        return [final_esp, df['Atom'][idx_atom]]


    def calc_firstTermE(espatom_idx, charge_range, charge_file):
        '''
        Input: espatom_idx: integer of atom index
        charge_range: list of integers of atom indices
        charge_file: string of charge filename
        Output: list of E-field vector and atomic symbol
        '''

        # E in units of V/(ansgrom) = N*m/(C*Angstrom)
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Convert each column to list for quicker indexing
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']

        # Pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        # Determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        position_vec = [xo, yo, zo]
        Ex = 0
        Ey = 0
        Ez = 0

        for idx in charge_range:
            if idx == idx_atom:
                continue
            else:
                # Calculate Efield and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                Ex = Ex - constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(charges[idx]/r)*(1/(xs[idx] - xo))
                Ey = Ey - constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(charges[idx]/r)*(1/(ys[idx] - yo))
                Ez = Ez - constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(charges[idx]/r)*(1/(zs[idx] - zo))

        E_vec = [Ex, Ey, Ez]
        return [E_vec, position_vec, df['Atom'][idx_atom]]



    def calc_firstTermE_atom_decomposable(self, espatom_idx, charge_range, charge_file, df_ptchg=None):
        '''
        Input: espatom_idx: integer of atom index
        charge_range: list of integers of atom indices
        charge_file: string of charge filename
        df_ptchg: optional pre-loaded point charge dataframe
        Output: list of E-field vector and atomic symbol
        '''

        inv_eps = 1/(self.config['dielectric'])
        # E in units of V/(ansgrom) = N*m/(C*Angstrom)
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Convert each column to list for quicker indexing
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']

        # Pick the index of the atom at which the efield should be calculated
        idx_atom = espatom_idx

        # Determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        position_vec = constants.ANGSTROM_TO_M*np.array([xo, yo, zo])
        Ex = 0
        Ey = 0
        Ez = 0 
        atom_wise_additions = []
        for idx in range(0, len(xs)):
            if idx == idx_atom:
                atom_wise_additions.append([0,0,0])
                continue
            elif idx not in charge_range:
                atom_wise_additions.append([0,0,0])
                continue
            else:
                # Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                Ex_contrib = -inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(charges[idx])*((xs[idx] - xo)*constants.ANGSTROM_TO_M)/(r**3)
                Ey_contrib = -inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(charges[idx])*((ys[idx] - yo)*constants.ANGSTROM_TO_M)/(r**3)
                Ez_contrib = -inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(charges[idx])*((zs[idx] - zo)*constants.ANGSTROM_TO_M)/(r**3)
                Ex = Ex + Ex_contrib
                Ey = Ey + Ey_contrib
                Ez = Ez + Ez_contrib
                atom_wise_additions.append([Ex_contrib, Ey_contrib,Ez_contrib ])


        E_vec = [Ex, Ey, Ez]
        QM_charges = charges
                #For QMMM calculation, include point charges in E field calc
        if self.config['includePtChgs']:
            qm_charges = charges
            QM_coords = [xs, ys, zs]
            #check if self.ptchf is population
            if df_ptchg is None:
                # If point charges are enabled but not passed, try to load them
                # Note: This may fail if ptChgfp is not a full path
                df_ptchg = self.getPtChgs(self.config['ptChgfp'])
            MM_xs = list(df_ptchg['x'])
            MM_ys = list(df_ptchg['y'])
            MM_zs = list(df_ptchg['z'])
            init_MM_charges = list(df_ptchg['charge'])

            #make new MM_charges by scaling!
            mm_coords = np.column_stack((np.array(MM_xs), np.array(MM_ys), np.array(MM_zs)))
            mm_charges = np.array(init_MM_charges)
            qm_charges = np.array(QM_charges)
            #MM_charges = Electrostatics.correct_mm_charges(QM_coords, qm_charges, mm_coords, mm_charges)

            MM_charges = mm_charges*(1/(math.sqrt(self.config['dielectric_scale'])))
            #change charge range to include these new partial charges!
            charge_range = range(0, len(MM_xs))
            for chg_idx in charge_range:
                r = (((MM_xs[chg_idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((MM_ys[chg_idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((MM_zs[chg_idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                dist_vec = constants.ANGSTROM_TO_M*np.array([(MM_xs[chg_idx] - xo), (MM_ys[chg_idx] - yo), (MM_zs[chg_idx] - zo)])
                E_to_add = -inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*MM_charges[chg_idx]*(1/(r**3))*dist_vec
                E_vec[0] += E_to_add[0]
                E_vec[1] += E_to_add[1]
                E_vec[2] += E_to_add[2]

        return [constants.VM_TO_VA*np.array(E_vec), position_vec, df['Atom'][idx_atom], constants.VM_TO_VA*np.array(atom_wise_additions)]

    #Helper functions for resp-based adjustement of MD point charges

    def compute_esp(q_coords, q_vals, probe_coords):
        """
        Compute ESP at probe_coords due to point charges at q_coords.
        """
        diff = probe_coords[:, np.newaxis, :] - q_coords[np.newaxis, :, :] 
        dists = np.linalg.norm(diff, axis=2)  # shape: (P, Q)
        mask = dists > 1e-8 
        # Avoid division by zero by masking
        inv_dists = np.zeros_like(dists)
        inv_dists[mask] = 1.0 / dists[mask]
        
        esp = np.dot(inv_dists, q_vals) 
        return esp
    def compute_esp_from_qm(qm_coords, qm_charges, mm_coords):
        '''
        Compute the ESP at MM atoms due to QM charges.
        '''
        return Electrostatics.compute_esp(qm_coords, qm_charges, mm_coords)



    def calc_fullE(self, idx_atom, charge_range, xyz_file, atom_multipole_file, df_ptchg=None):
        '''

        Input: idx_atom: integer of atom index
        charge_range: list of integers of atom indices
        xyz_file: string of xyz filename
        atom_multipole_file: string of multipole filename
        df_ptchg: optional pre-loaded point charge dataframe
        Output: list of E-field vector and atomic symbol, Efields are reported in units of Volts/Angstrom
        '''

        df = Geometry(xyz_file).getGeomInfo()
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Convert each column to list for quicker indexing
        xs = df['X']
        ys = df['Y']
        zs = df['Z']

        
        # Following derivation of E-field strength 
        Monopole_E = np.array([0, 0, 0])

        # Unit conversion
        # Bohr to meters (atomic units)
        b_to_A = 0.529177
        Ex = 0
        Ey = 0
        Ez = 0

        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]

        position_vec = constants.ANGSTROM_TO_M*np.array([xo, yo, zo])


        inv_eps = 1/self.config['dielectric']

        if self.config['includePtChgs'] and df_ptchg is None:
            # If point charges are enabled but not passed, try to load them
            # Note: This may fail if ptChgfp is not a full path
            df_ptchg = self.getPtChgs(self.config['ptChgfp'])
        #load multipole moments from processed outputs 
        lst_multipole_dict = MultiwfnInterface.getmultipoles(atom_multipole_file)

        #make a list for each term in multipole expansion
        lst_multipole_idxs = list(range(0, len(lst_multipole_dict)))

       
        QM_charges = [lst_multipole_dict[idx]["Atom_Charge"] for idx in range(0, len(xs))]
        QM_coords = np.column_stack((np.array(xs), np.array(ys), np.array(zs)))
        #This ensure that atoms we would like to exclude from Efield calculation are excluded from every term in the multipole expansion
        multipole_chg_range = [i for i in lst_multipole_idxs if i in charge_range]
        for idx in multipole_chg_range:
            #If the atom idx is outside of the charge range then skip
            atom_dict = lst_multipole_dict[idx] 
            if idx == idx_atom:
                continue
            else:
                # Units of dipole_vec: 
                dipole_vec = constants.BOHR_TO_M*np.array(atom_dict["Dipole_Moment"])
                #convention of vector pointing towards atom of interest, so positive charges exert positive Efield
                dist_vec = constants.ANGSTROM_TO_M*np.array([(xs[idx] - xo), (ys[idx] - yo), (zs[idx] - zo)])
                dist_arr = np.outer(dist_vec, dist_vec)
                quadrupole_arr = constants.BOHR_TO_M*constants.BOHR_TO_M*atom_dict['Quadrupole_Moment']
                
                # Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                Monopole_E = Monopole_E - inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"]*(1/(r**3))*dist_vec  #neg for sign convention so Efield points toward neg charge
                Ex_quad = inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(1/(r**3))*(constants.ANGSTROM_TO_M*(xs[idx] - xo))*(1/r**4)*dist_arr[1:, 1:]*quadrupole_arr[1:,1:]
                Ey_quad = inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(1/(r**3))*(constants.ANGSTROM_TO_M*(ys[idx] - yo))*(1/r**4)*dist_arr[0:2:3, 0:2:3]*quadrupole_arr[0:2:3,0:2:3]
                Ez_quad = inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(1/(r**3))*(constants.ANGSTROM_TO_M*(zs[idx] - zo))*(1/r**4)*dist_arr[0:2, 0:2]*quadrupole_arr[0:2,0:2]

                Ex = Ex + inv_eps*constants.COULOMB_CONSTANT*(1/(r**3))*dist_vec[0]*( -constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"] + constants.ELEMENTARY_CHARGE*(1/r**2)*np.dot(dipole_vec[1:], dist_vec[1:]))-(1/3)*Ex_quad.sum()
                Ey = Ey + inv_eps*constants.COULOMB_CONSTANT*(1/(r**3))*dist_vec[1]*( -constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"] + constants.ELEMENTARY_CHARGE*(1/r**2)*np.dot(dipole_vec[0:2:3], dist_vec[0:2:3])) -(1/3)*Ey_quad.sum()
                Ez = Ez + inv_eps*constants.COULOMB_CONSTANT*(1/(r**3))*dist_vec[2]*( -constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"] + constants.ELEMENTARY_CHARGE*(1/r**2)*np.dot(dipole_vec[0:2], dist_vec[0:2])) -(1/3)*Ez_quad.sum()
        E_vec = [Ex, Ey, Ez]


        #For QMMM calculation, include point charges in E field calc
        if self.config['includePtChgs']:
            MM_xs = list(df_ptchg['x'])
            MM_ys = list(df_ptchg['y'])
            MM_zs = list(df_ptchg['z'])
            init_MM_charges = list(df_ptchg['charge'])
 
            #make new MM_charges by scaling!
            mm_coords = np.column_stack((np.array(MM_xs), np.array(MM_ys), np.array(MM_zs)))
            mm_charges = np.array(init_MM_charges)
            qm_charges = np.array(QM_charges)
            print(f'Size of mm_coords: {np.shape(mm_coords)} and qm coords: {np.shape(QM_coords)}; mm charges: {np.shape(mm_charges)} and qm_charges: {np.shape(qm_charges)}')
            #MM_charges = Electrostatics.correct_mm_charges(QM_coords, qm_charges, mm_coords, mm_charges)

            MM_charges = mm_charges*(1/(math.sqrt(self.config['dielectric_scale'])))
            #change charge range to include these new partial charges!
            charge_range = range(0, len(MM_xs))
            for chg_idx in charge_range:
                r = (((MM_xs[chg_idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((MM_ys[chg_idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((MM_zs[chg_idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                dist_vec = constants.ANGSTROM_TO_M*np.array([(MM_xs[chg_idx] - xo), (MM_ys[chg_idx] - yo), (MM_zs[chg_idx] - zo)])
                E_to_add = -inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*MM_charges[chg_idx]*(1/(r**3))*dist_vec
                Monopole_E = Monopole_E +  E_to_add
                E_vec[0] += E_to_add[0]
                E_vec[1] += E_to_add[1]
                E_vec[2] += E_to_add[2]
            #Add contributions to Monopole E from point charges to total E
        return [constants.VM_TO_VA*np.array(E_vec), position_vec, df['Atom'][idx_atom], constants.VM_TO_VA*np.array(Monopole_E) ]



    def calc_atomwise_ElectricField(self, idx_atom, charge_range, xyz_file, atom_multipole_file, df_ptchg=None):
        '''

        Input: idx_atom: integer of atom index
        charge_range: list of integers of atom indices
        xyz_file: string of xyz filename
        atom_multipole_file: string of multipole filename
        df_ptchg: optional pre-loaded point charge dataframe
        Output: list of E-field vector and atomic symbol, Efields are reported in units of Volts/Angstrom
        '''
        df = Geometry(xyz_file).getGeomInfo()
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Convert each column to list for quicker indexing
        xs = df['X']
        ys = df['Y']
        zs = df['Z']

        
        # Following derivation of E-field strength 
        Monopole_E = np.array([0, 0, 0])

        # Unit conversion
        # Bohr to meters (atomic units)
        b_to_A = 0.529177
        Ex = 0
        Ey = 0
        Ez = 0

        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]

        position_vec = constants.ANGSTROM_TO_M*np.array([xo, yo, zo])

        inv_eps = 1/self.config['dielectric']

        if self.config['includePtChgs'] and df_ptchg is None:
            # If point charges are enabled but not passed, try to load them
            # Note: This may fail if ptChgfp is not a full path
            df_ptchg = self.getPtChgs(self.config['ptChgfp'])

        #load multipole moments from processed outputs
        lst_multipole_dict = MultiwfnInterface.getmultipoles(atom_multipole_file)

        #make a list for each term in multipole expansion
        # Create mapping: atom_index (0-indexed) -> multipole_dict
        # Note: Multiwfn uses 1-indexed atoms, so we subtract 1 for 0-indexed Python arrays
        # However, if the file is already 0-indexed, the Index field will be "0", "1", etc.
        # We'll use enumerate to ensure proper mapping regardless
        multipole_dict_by_idx = {idx: atom for idx, atom in enumerate(lst_multipole_dict)}

        # Create charge list in proper order
        QM_charges = [atom['Atom_Charge'] for atom in lst_multipole_dict]
        QM_coords = np.column_stack((np.array(xs), np.array(ys), np.array(zs)))
        #This ensure that atoms we would like to exclude from Efield calculation are excluded from every term in the multipole expansion
        # Use actual atom indices that exist in both geometry and multipole data
        multipole_chg_range = [i for i in range(len(lst_multipole_dict)) if i in charge_range]
        multipole_Efield_contribution = []
        for idx in multipole_chg_range:
            #If the atom idx is outside of the charge range then skip
            atom_dict = multipole_dict_by_idx[idx]
            if idx == idx_atom:
                # Append zero contribution for self (maintains consistent array size)
                multipole_Efield_contribution.append([0.0, 0.0, 0.0])
                continue
            else:
                # Units of dipole_vec:
                dipole_vec = constants.BOHR_TO_M*np.array(atom_dict["Dipole_Moment"])
                #convention of vector pointing towards atom of interest, so positive charges exert positive Efield
                dist_vec = constants.ANGSTROM_TO_M*np.array([(xs[idx] - xo), (ys[idx] - yo), (zs[idx] - zo)])
                dist_arr = np.outer(dist_vec, dist_vec)
                quadrupole_arr = constants.BOHR_TO_M*constants.BOHR_TO_M*atom_dict['Quadrupole_Moment']

                # Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                Monopole_E = Monopole_E - inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"]*(1/(r**3))*dist_vec  #neg for sign convention so Efield points toward neg charge
                Ex_quad = inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(1/(r**3))*(constants.ANGSTROM_TO_M*(xs[idx] - xo))*(1/r**4)*dist_arr[1:, 1:]*quadrupole_arr[1:,1:]
                Ey_quad = inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(1/(r**3))*(constants.ANGSTROM_TO_M*(ys[idx] - yo))*(1/r**4)*dist_arr[0:2:3, 0:2:3]*quadrupole_arr[0:2:3,0:2:3]
                Ez_quad = inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*(1/(r**3))*(constants.ANGSTROM_TO_M*(zs[idx] - zo))*(1/r**4)*dist_arr[0:2, 0:2]*quadrupole_arr[0:2,0:2]

                E_x_atomic_contrib = inv_eps*constants.COULOMB_CONSTANT*(1/(r**3))*dist_vec[0]*( -constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"] + constants.ELEMENTARY_CHARGE*(1/r**2)*np.dot(dipole_vec[1:], dist_vec[1:]))-(1/3)*Ex_quad.sum()
                E_y_atomic_contrib = inv_eps*constants.COULOMB_CONSTANT*(1/(r**3))*dist_vec[1]*( -constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"] + constants.ELEMENTARY_CHARGE*(1/r**2)*np.dot(dipole_vec[0:2:3], dist_vec[0:2:3])) -(1/3)*Ey_quad.sum()
                E_z_atomic_contrib = inv_eps*constants.COULOMB_CONSTANT*(1/(r**3))*dist_vec[2]*( -constants.ELEMENTARY_CHARGE*atom_dict["Atom_Charge"] + constants.ELEMENTARY_CHARGE*(1/r**2)*np.dot(dipole_vec[0:2], dist_vec[0:2])) -(1/3)*Ez_quad.sum()

                Ex = Ex + E_x_atomic_contrib
                Ey = Ey + E_y_atomic_contrib
                Ez = Ez + E_z_atomic_contrib
                # Append individual atomic contribution (not cumulative)
                multipole_Efield_contribution.append([E_x_atomic_contrib, E_y_atomic_contrib, E_z_atomic_contrib])

        E_vec = [Ex, Ey, Ez]


        #For QMMM calculation, include point charges in E field calc
        if self.config['includePtChgs']:
            MM_xs = list(df_ptchg['x'])
            MM_ys = list(df_ptchg['y'])
            MM_zs = list(df_ptchg['z'])
            init_MM_charges = list(df_ptchg['charge'])
 
            #make new MM_charges by scaling!
            mm_coords = np.column_stack((np.array(MM_xs), np.array(MM_ys), np.array(MM_zs)))
            mm_charges = np.array(init_MM_charges)
            qm_charges = np.array(QM_charges)
            #MM_charges = Electrostatics.correct_mm_charges(QM_coords, qm_charges, mm_coords, mm_charges)

            MM_charges = mm_charges*(1/(math.sqrt(self.config['dielectric_scale'])))
            #change charge range to include these new partial charges!
            charge_range = range(0, len(MM_xs))
            for chg_idx in charge_range:
                r = (((MM_xs[chg_idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((MM_ys[chg_idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((MM_zs[chg_idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                dist_vec = constants.ANGSTROM_TO_M*np.array([(MM_xs[chg_idx] - xo), (MM_ys[chg_idx] - yo), (MM_zs[chg_idx] - zo)])
                E_to_add = -inv_eps*constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*MM_charges[chg_idx]*(1/(r**3))*dist_vec
                Monopole_E = Monopole_E +  E_to_add
                E_vec[0] += E_to_add[0]
                E_vec[1] += E_to_add[1]
                E_vec[2] += E_to_add[2]
                multipole_Efield_contribution.append(E_to_add)
            #Add contributions to Monopole E from point charges to total E
        return [constants.VM_TO_VA*np.array(E_vec), position_vec, df['Atom'][idx_atom], constants.VM_TO_VA*np.array(Monopole_E), constants.VM_TO_VA*np.array(multipole_Efield_contribution)]
    def Efield_multipleAtoms(self, xyzfilepath, df_substrate, df_env, return_decomposition=False):
        """
        Calculate E-field vectors at substrate atoms due to environment atoms.

        Parameters:
        -----------
        xyzfilepath : str
            Path to XYZ file
        df_substrate : DataFrame
            DataFrame with substrate atom information (must have 'Index', 'Atom_Charge' columns)
        df_env : DataFrame
            DataFrame with environment atom information (must have 'Index', 'Atom_Charge' columns)
        return_decomposition : bool, optional
            If True, return atom-wise E-field contributions (default: False)

        Returns:
        --------
        DataFrame or tuple
            If return_decomposition=False: df_substrate with added 'Efield_x', 'Efield_y', 'Efield_z' columns
            If return_decomposition=True: (df_substrate, atomwise_Efield_dict)
                where atomwise_Efield_dict has (sub_idx, env_idx) keys and [Ex, Ey, Ez] vector values
        """
        df = Geometry(xyzfilepath).getGeomInfo()

        # Convert each column to list for quicker indexing
        atoms = df['Atom']
        xs = df['X']
        ys = df['Y']
        zs = df['Z']

        dielectric = self.config['dielectric']
        inv_eps = 1/dielectric

        substrate_idxs = df_substrate['Index']
        env_idxs = df_env['Index']

        # Store total E-fields at each substrate atom
        total_Efields = []  # List of [Ex, Ey, Ez] vectors
        # Store atom-wise E-field contributions: dict with (sub_idx, env_idx) as key and [Ex, Ey, Ez] as value
        atomwise_Efield_contributions = {}

        for sub_idx in substrate_idxs:
            E_vec = np.array([0.0, 0.0, 0.0])
            xo = xs[sub_idx]
            yo = ys[sub_idx]
            zo = zs[sub_idx]

            for idx in env_idxs:
                # Distance vector from environment atom to substrate atom
                dist_vec = constants.ANGSTROM_TO_M * np.array([xs[idx] - xo, ys[idx] - yo, zs[idx] - zo])
                r = np.linalg.norm(dist_vec)

                # E-field from a point charge: E = k * q / r^2 * r_hat
                # Direction: points away from positive charge, toward negative charge
                env_charge = df_env.loc[df_env['Index'] == idx, 'Atom_Charge'].iloc[0]
                E_contribution = -inv_eps * constants.COULOMB_CONSTANT * constants.ELEMENTARY_CHARGE * env_charge * (1/(r**3)) * dist_vec

                E_vec += E_contribution

                # Store individual contribution
                atomwise_Efield_contributions[(sub_idx, idx)] = E_contribution

            total_Efields.append(E_vec)

        df_substrate = df_substrate.copy()
        # Convert to V/Angstrom for consistency with other E-field calculations
        Efields_array = constants.VM_TO_VA * np.array(total_Efields)
        df_substrate['Efield_x'] = Efields_array[:, 0]
        df_substrate['Efield_y'] = Efields_array[:, 1]
        df_substrate['Efield_z'] = Efields_array[:, 2]
        df_substrate['Efield_magnitude'] = np.linalg.norm(Efields_array, axis=1)

        print(f'Substrate E-fields (V/Å): {Efields_array}')

        if return_decomposition:
            # Convert atomwise contributions to V/Angstrom as well
            atomwise_Efield_VA = {key: constants.VM_TO_VA * vec for key, vec in atomwise_Efield_contributions.items()}
            return df_substrate, atomwise_Efield_VA
        return df_substrate

    def ESP_multipleAtoms(self, xyzfilepath, df_substrate, df_env, return_decomposition=False):
        #will return a list of multipole moments where each is
        df = Geometry(xyzfilepath).getGeomInfo()
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Convert each column to list for quicker indexing
        atoms = df['Atom']
        xs = df['X']
        ys = df['Y']
        zs = df['Z']

        #KJ_J = 10**-3
        #faraday = 23.06   #kcal/(mol*V)
        one_mol = 6.02*(10**23)

        dielectric = self.config['dielectric']

        substrate_idxs = df_substrate['Index']
        env_idxs = df_env['Index']

        env_charges = df_env['Atom_Charge']
        #env_dipoles = df_env['Dipole_Moment']
        #env_quadrupole = df_env['Quadrupole_Moment']

        total_ESPs = []
        # Store atom-wise ESP contributions: dict with (sub_idx, env_idx) as key and ESP contribution as value
        atomwise_ESP_contributions = {}

        for sub_idx in substrate_idxs:
            sub_esp = 0
            xo = xs[sub_idx]
            yo = ys[sub_idx]
            zo = zs[sub_idx]
            for idx in env_idxs:
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                env_contribution = (constants.ELEMENTARY_CHARGE*constants.COULOMB_CONSTANT*df_env.loc[df_env['Index'] == idx, 'Atom_Charge'].iloc[0])/r
                sub_esp = sub_esp + env_contribution
                # Store individual contribution
                atomwise_ESP_contributions[(sub_idx, idx)] = env_contribution
            total_ESPs.append(sub_esp)
        df_substrate = df_substrate.copy()
        df_substrate['ESP'] = np.array(total_ESPs) #Units: N*m^2/(C^2)*(C/m) = N*m/C = J/C = Volt
        print(f'Here is the substrate ESP:{total_ESPs} in volts')

        if return_decomposition:
            return df_substrate, atomwise_ESP_contributions
        return df_substrate

    def ESPfromMultipole(self, xyzfilepath, atom_multipole_file, charge_range, idx_atom):
        '''

        Input: idx_atom: integer of atom index
        charge_range: list of integers of atom indices
        xyfilepath: string of xyz filename
        atom_multipole_file: string of multipole filename
        Output: list of ESP and atomic symbol Units of ESP is Volts!
        '''
        df = Geometry(xyzfilepath).getGeomInfo()
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Convert each column to list for quicker indexing
        atoms = df['Atom']
        xs = df['X']
        ys = df['Y']
        zs = df['Z']

        # Determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]

         # Unit conversion
        #KJ_J = 10**-3
        #faraday = 23.06   #kcal/(mol*V)
        one_mol = 6.02*(10**23)
        #cal_J = 4.184
        dielectric = self.config['dielectric']
        lst_multipole_dict = MultiwfnInterface.getmultipoles(atom_multipole_file)
        #create list of bound atoms, these are treated with a different dielectric
        bound_atoms = Geometry(xyzfilepath).getBondedAtoms(idx_atom)
        total_esp= 0
        for idx in charge_range:
            atom_dict = lst_multipole_dict[idx]
            if idx == idx_atom:
                continue
            elif idx in bound_atoms:
                #default it to exclude bound atoms
                 r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                 total_esp = total_esp + (atom_dict["Atom_Charge"]/r)
            else:
                # Calculate esp and convert to units (A to m)
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                total_esp = total_esp + (1/dielectric)*(atom_dict["Atom_Charge"]/r)

        final_esp = constants.COULOMB_CONSTANT*total_esp*((constants.ELEMENTARY_CHARGE))  #Units: N*m^2/(C^2)*(C/m) = N*m/C = J/C = Volt
        return [final_esp, df['Atom'][idx_atom] ]
 

 # Bond_indices is a list of tuples where each tuple contains the zero-indexed values of location of the atoms of interest
    def E_proj_bondIndices(self, bond_indices, xyz_filepath, atom_multipole_file, all_lines):
        '''
        Input: bond_indices: list of tuples of atom indices
        xyz_filepath: string of xyz filename
        atom_multipole_file: string of multipole filename
        all_lines: list of strings of lines in xyz file
        Output: list of E-field vector and atomic symbol in units of V/Angstrom
        '''
        bonded_atoms = []
        E_projected = []
        E_monopole_proj = []
        bonded_positions = []
        # Determine the Efield vector at point of central metal stom
        bond_lens = []

        
        for atomidxA, atomidxB in bond_indices:
            [A_bonded_E, A_bonded_position, A_bonded_atom, A_monopole_E_bonded]  =  self.calc_fullE(atomidxA, all_lines, xyz_filepath, atom_multipole_file)   
            [B_bonded_E, B_bonded_position, B_bonded_atom, B_monopole_E_bonded]  =  self.calc_fullE(atomidxB, all_lines, xyz_filepath, atom_multipole_file)  
            bond_vec_unnorm = np.subtract(np.array(A_bonded_position), np.array(B_bonded_position)) 
            bond_len = np.linalg.norm(bond_vec_unnorm)
            bond_vec = bond_vec_unnorm/(bond_len)
            # Initialized a bond_dipole_vec as the (bond_vec_unnorm )*(sum of the partial charges).. can just use dipole! 
            # Compute E-field projected along this bond!
            E_proj = (1/2)*np.dot((np.array(A_bonded_E) + np.array(B_bonded_E)), bond_vec)
            E_proj_monopole = (1/2)*np.dot((np.array(A_monopole_E_bonded) + np.array(B_monopole_E_bonded)), bond_vec)
            E_projected.append(E_proj)
            E_monopole_proj.append(E_proj_monopole)
            bonded_atoms.append((A_bonded_atom, B_bonded_atom))
            bonded_positions.append((A_bonded_position, B_bonded_position))
            bond_lens.append(bond_len)
        return [E_projected, bonded_atoms, bond_indices, bond_lens, E_monopole_proj]


    def calc_bond_dipoles(self, bond_indices, xyz_filepath, atom_multipole_file, bool_multipole):
        """Calculate bond dipole moments for given bond indices.

        For each bond between atoms A and B, calculates the bond dipole moment vector
        relative to the geometric center (midpoint) of the bond. This is origin-independent
        for both neutral and charged bonds:
        μ_bond = q_A * (r_A - r_center) + q_B * (r_B - r_center)
        which simplifies to: μ_bond = (q_A - q_B) * (r_A - r_B) / 2

        Parameters
        ----------
        bond_indices : list of tuples
            List of (atomA_idx, atomB_idx) tuples specifying bonds
        xyz_filepath : str
            Path to XYZ structure file
        atom_multipole_file : str
            Path to multipole/charge file
        bool_multipole : bool
            If True, use multipole file; if False, use monopole charges

        Returns
        -------
        list
            List of bond dipole information dictionaries containing:
            - 'bond_dipole_vec': [x, y, z] dipole vector in Debye
            - 'bond_dipole_mag': magnitude of dipole in Debye
            - 'bond_indices': (atomA_idx, atomB_idx)
            - 'charges': (charge_A, charge_B)
        """
        bond_dipoles = []

        # Get geometry information
        df_geom = Geometry(xyz_filepath).getGeomInfo()

        # Get charge information
        if bool_multipole:
            multipole_data = MultiwfnInterface.getmultipoles(atom_multipole_file)
            # Use enumerate for proper 0-indexed mapping, regardless of Index field value
            charge_dict = {idx: atom['Atom_Charge'] for idx, atom in enumerate(multipole_data)}
        else:
            df_charges = pd.read_csv(atom_multipole_file, sep='\s+', names=["Atom", 'x', 'y', 'z', "charge"])
            charge_dict = {idx: charge for idx, charge in enumerate(df_charges['charge'])}

        for atomA_idx, atomB_idx in bond_indices:
            # Get positions (in Angstrom from XYZ file)
            posA = np.array([df_geom.loc[atomA_idx, 'X'],
                           df_geom.loc[atomA_idx, 'Y'],
                           df_geom.loc[atomA_idx, 'Z']])
            posB = np.array([df_geom.loc[atomB_idx, 'X'],
                           df_geom.loc[atomB_idx, 'Y'],
                           df_geom.loc[atomB_idx, 'Z']])

            # Get charges (in elementary charge units)
            chargeA = charge_dict[atomA_idx]
            chargeB = charge_dict[atomB_idx]

            # Calculate bond dipole relative to geometric center (midpoint)
            # This makes the dipole origin-independent for both neutral and charged bonds
            # and consistent with line integral formalism
            r_center = (posA + posB) / 2
            dipole_vec_eA = chargeA * (posA - r_center) + chargeB * (posB - r_center)  # in e*Angstrom
            # Simplified: dipole_vec = (chargeA - chargeB) * (posA - posB) / 2

            # Convert to Debye (1 e*Angstrom = 4.80320427 Debye)
            convert_to_Debye = 4.80320427
            dipole_vec = convert_to_Debye * dipole_vec_eA  # in Debye
            dipole_mag = np.linalg.norm(dipole_vec)  # in Debye

            bond_dipole_info = {
                'bond_dipole_vec': dipole_vec.tolist(),
                'bond_dipole_mag': dipole_mag,
                'bond_indices': (atomA_idx, atomB_idx),
                'charges': (chargeA, chargeB)
            }
            bond_dipoles.append(bond_dipole_info)

        return bond_dipoles

    # Bond_indices is a list of tuples where each tuple contains the zero-indexed values of location of the atoms of interest
    def E_proj_bondIndices_atomwise(self, bond_indices, xyz_filepath, atom_multipole_file, all_lines, bool_multipole, df_ptchg=None):
        '''
        Input: bond_indices: list of tuples of atom indices
        xyz_filepath: string of xyz filename
        atom_multipole_file: string of multipole filename
        all_lines: list of strings of lines in xyz file
        df_ptchg: optional pre-loaded point charge dataframe
        Output: list of E-field vector and atomic symbol in units of V/Angstrom
        '''
        bonded_atoms = []
        E_projected = []
        E_monopole_proj = []
        bonded_positions = []
        # Determine the Efield vector at point of central metal stom
        bond_lens = []
        E_proj_atomwise_list = []  # List to store per-bond atom-wise contributions

        if bool_multipole:
            for atomidxA, atomidxB in bond_indices:
                [A_bonded_E, A_bonded_position, A_bonded_atom, A_monopole_E_bonded, A_atom_wise_E]  =  self.calc_atomwise_ElectricField(atomidxA, all_lines, xyz_filepath, atom_multipole_file, df_ptchg=df_ptchg)
                [B_bonded_E, B_bonded_position, B_bonded_atom, B_monopole_E_bonded, B_atom_wise_E]  =  self.calc_atomwise_ElectricField(atomidxB, all_lines, xyz_filepath, atom_multipole_file, df_ptchg=df_ptchg)
                bond_vec_unnorm = np.subtract(np.array(A_bonded_position), np.array(B_bonded_position))
                bond_len = np.linalg.norm(bond_vec_unnorm)
                bond_vec = bond_vec_unnorm/(bond_len)
                # Initialized a bond_dipole_vec as the (bond_vec_unnorm )*(sum of the partial charges).. can just use dipole!
                # Compute E-field projected along this bond!
                E_proj = (1/2)*np.dot((np.array(A_bonded_E) + np.array(B_bonded_E)), bond_vec)
                E_proj_monopole = (1/2)*np.dot((np.array(A_monopole_E_bonded) + np.array(B_monopole_E_bonded)), bond_vec)
                E_projected.append(E_proj)
                E_monopole_proj.append(E_proj_monopole)
                bonded_atoms.append((A_bonded_atom, B_bonded_atom))
                bonded_positions.append((A_bonded_position, B_bonded_position))
                bond_lens.append(bond_len)
                E_proj_atomwise = (1/2)*(np.array(A_atom_wise_E) + np.array(B_atom_wise_E))@ bond_vec.T

                E_proj_atomwise_list.append(E_proj_atomwise)

        else:
            for atomidxA, atomidxB in bond_indices:
                [A_bonded_E, A_bonded_position, A_bonded_atom,  A_atom_wise_E]  =  self.calc_firstTermE_atom_decomposable(atomidxA, all_lines, atom_multipole_file, df_ptchg=df_ptchg)
                [B_bonded_E, B_bonded_position, B_bonded_atom, B_atom_wise_E]  =  self.calc_firstTermE_atom_decomposable(atomidxB, all_lines, atom_multipole_file, df_ptchg=df_ptchg)

                bond_vec_unnorm = np.subtract(np.array(A_bonded_position), np.array(B_bonded_position))
                bond_len = np.linalg.norm(bond_vec_unnorm)
                bond_vec = bond_vec_unnorm/(bond_len)
                # Initialized a bond_dipole_vec as the (bond_vec_unnorm )*(sum of the partial charges).. can just use dipole!
                # Compute E-field projected along this bond!
                E_proj = (1/2)*np.dot((np.array(A_bonded_E) + np.array(B_bonded_E)), bond_vec)
                E_projected.append(E_proj)
                bonded_atoms.append((A_bonded_atom, B_bonded_atom))
                bonded_positions.append((A_bonded_position, B_bonded_position))
                bond_lens.append(bond_len)

                E_proj_atomwise = (1/2)*((np.array(A_atom_wise_E) + np.array(B_atom_wise_E))@ bond_vec.T)

                E_proj_atomwise_list.append(E_proj_atomwise)

        # For backwards compatibility, return the last bond's atomwise data as the 5th element
        # and return the full list as the 6th element
        E_proj_atomwise = E_proj_atomwise_list[-1] if E_proj_atomwise_list else None
        return [E_projected, bonded_atoms, bond_indices, bond_lens, E_proj_atomwise, E_proj_atomwise_list]

    

    def E_proj_first_coord(self, metal_idx, xyz_file_path, atom_multipole_file, all_lines):
        '''Function to calculate Efield projection accounting only for electrostatic effects of directly bound atoms 
        Input: metal_idx: integer of atom index
        xyz_file_path: string of xyz filename
        atom_multipole_file: string of multipole filename
        all_lines: list of strings of lines in xyz file
        Output: list of E-field vector and atomic symbol
        '''
        bonded_atoms = []
        E_projected = []
        E_monopole_proj = []
        bonded_positions = []
        # Determine the Efield vector at point of central metal stom
        [center_E, center_position, center_atom, Monopole_E_center]  =  self.calc_fullE(metal_idx, all_lines, xyz_file_path, atom_multipole_file)
        lst_bonded_atoms = self.getBondedAtoms(xyz_file_path, metal_idx) 
        bond_lens = []
        for bonded_atom_idx in lst_bonded_atoms:
            [bonded_E, bonded_position, bonded_atom, monopole_E_bonded]  =  self.calc_fullE(bonded_atom_idx, all_lines, xyz_file_path, atom_multipole_file)    
            bond_vec_unnorm = np.subtract(np.array(center_position), np.array(bonded_position)) 
            bond_len = np.linalg.norm(bond_vec_unnorm)
            bond_vec = bond_vec_unnorm/(bond_len)

            # Initialized a bond_dipole_vec as the (bond_vec_unnorm )*(sum of the partial charges).. can just use dipole! 
            # Compute E-field projected along this bond!
            E_proj = (1/2)*np.dot((np.array(bonded_E) + np.array(center_E)), bond_vec)
            E_proj_monopole= (1/2)*np.dot((Monopole_E_center + monopole_E_bonded), bond_vec)
            E_projected.append(E_proj)
            E_monopole_proj.append(E_proj_monopole)
            bonded_atoms.append(bonded_atom)
            bonded_positions.append(bonded_position)
            bond_lens.append(bond_len)
        return [E_projected, bonded_atoms, lst_bonded_atoms, bond_lens, E_monopole_proj]
            
    def esp_coord_shell(self, metal_idx, charge_file, path_to_xyz, n_shells=1):
        """Calculate ESP accounting for electrostatic contributions from n coordination shells.

        This unified method replaces esp_first_coord() and esp_second_coord() by accepting
        a parameter to specify how many coordination shells to include.

        Parameters
        ----------
        metal_idx : int
            Atom index of the central metal
        charge_file : str
            Path to charge filename
        path_to_xyz : str
            Path to xyz structure file
        n_shells : int, optional
            Number of coordination shells to include (default=1)
            n_shells=1: first coordination shell only
            n_shells=2: first and second coordination shells

        Returns
        -------
        float
            ESP value for the specified coordination shell(s)
        """
        print(f'The index of the metal atom is: {metal_idx}')

        # Build list of atoms to include based on n_shells
        atoms_to_include = []

        if n_shells >= 1:
            # First coordination shell: atoms directly bonded to metal
            first_shell_atoms = self.getBondedAtoms(path_to_xyz, metal_idx)
            atoms_to_include.extend(first_shell_atoms)

        if n_shells >= 2:
            # Second coordination shell: atoms bonded to first shell atoms
            for first_shell_atom in first_shell_atoms:
                second_shell_atoms = self.getBondedAtoms(path_to_xyz, first_shell_atom)
                atoms_to_include.extend(second_shell_atoms)

        if n_shells >= 3:
            # Third and higher coordination shells can be added here if needed
            raise NotImplementedError(f"Coordination shells beyond n=2 are not yet implemented (requested n={n_shells})")

        # Remove duplicates while preserving order
        unique_atoms = list(dict.fromkeys(atoms_to_include))

        # Calculate ESP
        [coord_shell_ESP, atom_type] = self.calcesp(path_to_xyz, metal_idx, unique_atoms, charge_file)
        return coord_shell_ESP


    def charge_atom(filename, atom_idx):
        '''Compute charges and return total charge and partial charge of atom
        Input: filename: string of charge filename
        atom_idx: integer of atom index
        Output: list of total charge and partial charge of atom'''
        df = pd.read_csv(filename, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']
        total_charge = np.sum(df['charge'])
        partial_charge_atom = charges[atom_idx]
        return [total_charge, partial_charge_atom]

    #This just gets the charge if the charge is in monopole setup
    def charge_atoms(self, chg_filename, xyzfilename):
        '''Input: filename: string of charge filename or multipole filename... if the xyzfilename is included must be multipole! '''
        if len(xyzfilename) > 1:
            dict_multipole = MultiwfnInterface.getmultipoles(chg_filename)  #Index, Element, Atom_Charge, Dipole_Moment....
            df_multipole = pd.DataFrame(dict_multipole)
            df_geom = Geometry(xyzfilename).getGeomInfo()
            merged_df = pd.merge(df_geom, df_multipole, left_index=True, right_index=True)
            merged_df['charge'] = merged_df['Atom_Charge']
            merged_df['x'] = merged_df['X']
            merged_df['y'] = merged_df['Y']
            merged_df['z'] = merged_df['Z']
            df = merged_df[["Atom",'x', 'y', 'z', "charge"]]
        else:
            df = pd.read_csv(chg_filename, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        return df


    def getAtomInfo(filename, atom_idx):
        '''
        Input: filename: string of charge filename
        atom_idx: integer of atom index
        Output: list of total charge and partial charge of atom'''
        df = pd.read_csv(filename, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']
        total_charge = np.sum(df['charge'])
        partial_charge_atom = charges[atom_idx]
        atom_type = atoms[atom_idx]
        return [atom_type, partial_charge_atom]


    def getAtomsInfo(filename, atom_indices):
        ''' Input: filename: string of charge filename
        atom_indices: list of integers of atom indices
        Output: list of total charge and partial charge of atom'''
        df = pd.read_csv(filename, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']
        total_charge = np.sum(df['charge'])
        res_charge = 0
        for atm_idx in atom_indices:
            print(f'Atom index: {atm_idx} and charge: {charges[atm_idx]}')
            res_charge = res_charge + charges[atm_idx]
        return [res_charge]


    def esp_bydistance(self, path_to_xyz, espatom_idx,  charge_file):
        '''
        Input: path_to_xyz: string of xyz filename
        espatom_idx: integer of atom index
        charge_file: string of charge filename
        Output: list of ESP and atomic symbol
        '''
        dielectric = self.config['dielectric']
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])

        # Unit conversion
        KJ_J = 10**-3
        faraday = 23.06   #kcal/(mol*V)
        one_mol = 6.02*(10**23)
        cal_J = 4.184

        # Convert each column to list for quicker indexing
        atoms = list(df['Atom'])
        charges = list(df['charge'])
        xs = list(df['x'])
        ys = list(df['y'])
        zs = list(df['z'])

        # Pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        # Determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        chargeo = charges[idx_atom]

        #For QMMM calculation, include point charges in ESP calculation
        if self.config['includePtChgs']:
            ptchg_filename = self.config['ptChgfp']
            if not ptchg_filename or ptchg_filename.strip() == '':
                raise ValueError(
                    "Point charge file path is not set. Please call set_ptChgfile() "
                    "with the filename (e.g., 'ptchg.dat') before using includePtChgs=True"
                )
            init_file_path = path_to_xyz[0:-len(self.folder_to_file_path + self.config['xyzfilename'])]
            full_ptchg_fp = init_file_path + ptchg_filename
            df_ptchg = self.getPtChgs(full_ptchg_fp)
            xs = xs + list(df_ptchg['x'])
            ys = ys + list(df_ptchg['y'])
            zs = zs + list(df_ptchg['z'])
            charges = charges + list(df_ptchg['charge'])
            atoms = atoms + list(df_ptchg['Atom'])

        total_esp = 0
        # Create an ordering of the atoms based on distance from the central atom
        total_atoms = len(xs)
        distances = []
        esps = []
        bound_atoms = []
        if self.config['changeDielectBoundBool']:
            bound_atoms = self.getBondedAtoms(path_to_xyz, idx_atom)
        for idx in range(0, total_atoms):
            if idx == idx_atom:
                continue
            elif idx in bound_atoms:
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                esps.append(constants.COULOMB_CONSTANT*constants.ELEMENTARY_CHARGE*charges[idx]/r) #units of N*m/C =Volt
                distances.append(r)
            else:
                r = (((xs[idx] - xo)*constants.ANGSTROM_TO_M)**2 + ((ys[idx] - yo)*constants.ANGSTROM_TO_M)**2 + ((zs[idx] - zo)*constants.ANGSTROM_TO_M)**2)**(0.5)
                distances.append(r)
                esps.append(constants.COULOMB_CONSTANT*(1/dielectric)*constants.ELEMENTARY_CHARGE*charges[idx]/r) #untis of N*m/C = Volt
        # Now we sort the distance list, and use sorted indices to sort the
        atm_lst = list(atoms)
        atm_lst.pop(idx_atom)

        chg_lst = list(charges)
        chg_lst.pop(idx_atom)
        chg_arr = np.array(chg_lst)
        
        dist_arr = np.array(distances)
        init_sorted_idx = np.argsort(dist_arr)
        lst_sorted_idx = list(init_sorted_idx)
        lst_sorted_idx.append(len(lst_sorted_idx))
        lst_sorted_idx.remove(idx_atom)
        sorted_idx = np.array(lst_sorted_idx)

        esp_arr = np.array(esps)
        sorted_esps = esp_arr.take(init_sorted_idx)
        cumulative_esps = np.cumsum(sorted_esps)

        sorted_dist = dist_arr.take(init_sorted_idx)
        sorted_partial_charges = chg_arr.take(init_sorted_idx)
        sorted_atomTypes = [atm_lst[i] for i in init_sorted_idx]
       
        return [sorted_dist, sorted_esps, cumulative_esps, sorted_idx, sorted_partial_charges, sorted_atomTypes]
  

        # list_of_folders = the list of the folders that contain the desired files
        # new_dir: the [post-folder path to the scr folder that contains the .molden and optim.xyz file themselfs
        # dict of calcs, calculations to be performed by multiwavefunction with the corresponding keys
        # newfilanme: desired name of the .csv fiole that will be createcd in getData cotnaining all of the ESP/other data extracted un the file

    # ============================================================================
    # CONSOLIDATED ESP AND E-FIELD METHODS
    # These unified methods replace the previous separate implementations
    # ============================================================================

    def getESP(self, charge_types, ESPdata_filename, multiwfn_module, multiwfn_path, atmrad_path,
               use_multipole=False, include_decay=False, include_coord_shells=False, dielectric=1):
        """Unified ESP calculation method supporting multiple modes.

        This consolidated method replaces getESPData(), getESPMultipole(), and getESPDecay()
        to eliminate code duplication while preserving all functionality.

        Parameters
        ----------
        charge_types : list of str or str
            Charge partitioning scheme(s). Options: 'Hirshfeld', 'Voronoi', 'Mulliken',
            'Lowdin', 'SCPA', 'Becke', 'ADCH', 'CHELPG', 'MK', 'AIM', 'Hirshfeld_I',
            'CM5', 'EEM', 'RESP', 'PEOE'. Can be a single string (for multipole mode)
            or list of strings.
        ESPdata_filename : str
            Output CSV filename (without extension).
        multiwfn_module : str
            Module name containing Multiwfn executable.
        multiwfn_path : str
            Path to Multiwfn executable.
        atmrad_path : str
            Path to atmrad executable.
        use_multipole : bool, optional
            If True, use multipole expansion; if False, use monopole charges (default: False).
        include_decay : bool, optional
            If True, include distance-sorted ESP decay analysis (default: False).
        include_coord_shells : bool, optional
            If True, include first and second coordination shell ESP (default: False).
        dielectric : float, optional
            Dielectric constant (default: 1).

        Returns
        -------
        pd.DataFrame
            DataFrame with ESP results saved to CSV.

        Examples
        --------
        >>> # Monopole mode (replaces getESPData)
        >>> es.getESP(['Hirshfeld', 'Hirshfeld_I'], 'esp_mono', ...)

        >>> # Multipole mode (replaces getESPMultipole)
        >>> es.getESP('Hirshfeld_I', 'esp_multi', ..., use_multipole=True)

        >>> # Decay analysis mode (replaces getESPDecay)
        >>> es.getESP(['Hirshfeld'], 'esp_decay', ..., include_decay=True, include_coord_shells=True)
        """
        # Ensure charge_types is a list
        if isinstance(charge_types, str):
            charge_types = [charge_types]

        # Validate that metal indices are provided for ESP calculation
        if not self.lst_of_tmcm_idx:
            raise ValueError(
                "ESP calculation requires metal atom indices. "
                "Please provide lst_of_tmcm_idx when initializing the Electrostatics object."
            )

        self.config['dielectric'] = dielectric
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders
        final_structure_file = self.config['xyzfilename']

        owd = os.getcwd()
        allspeciesdict = []
        counter = 0

        for f in list_of_file:
            print(f'-----------------{f}------------------')
            atom_idx = metal_idxs[counter]
            counter += 1
            results_dict = {}
            results_dict['Name'] = f
            file_path_xyz = f"{owd}/{f + folder_to_molden}{final_structure_file}"

            # Calculate total lines for monopole calculations
            if not use_multipole or include_decay:
                total_lines = Electrostatics.mapcount(file_path_xyz)
                init_all_lines = range(0, total_lines - 2)
                all_lines = [x for x in init_all_lines if x not in self.config['excludeAtomfromEcalc']]

            for charge_type in charge_types:
                print(f'Partial Charge Scheme: {charge_type}')
                try:
                    # Use centralized partitionCharge() to get/compute charges
                    comp_cost = self.multiwfn.partitionCharge(
                        multipole_bool=use_multipole,
                        f=f,
                        folder_to_molden=folder_to_molden,
                        multiwfn_path=multiwfn_path,
                        atmrad_path=atmrad_path,
                        charge_type=charge_type,
                        owd=owd
                    )

                    if comp_cost == -1:
                        print(f"WARNING: Charge calculation failed for {charge_type} in {f}")
                        continue

                    file_path_multipole = f"{owd}/{f + folder_to_molden}Multipole{charge_type}.txt"
                    file_path_monopole = f"{owd}/{f + folder_to_molden}Charges{charge_type}.txt"
                    path_to_xyz = f"{owd}/{f + folder_to_molden}{final_structure_file}"

                    # Calculate ESP using appropriate method
                    if use_multipole:
                        [final_esp, atom_name] = self.ESPfromMultipole(
                            file_path_xyz, file_path_multipole, all_lines, atom_idx
                        )
                        results_dict['Total ESP'] = final_esp
                        results_dict['Atom'] = atom_name
                    else:
                        [ESP_all, atom_type] = self.calcesp(
                            file_path_xyz, atom_idx, all_lines, file_path_monopole
                        )
                        [total_charge, partial_charge_atom] = Electrostatics.charge_atom(
                            file_path_monopole, atom_idx
                        )
                        results_dict['Atoms'] = atom_type
                        results_dict['Total Charge'] = total_charge
                        results_dict[f'Partial Charge {charge_type}'] = partial_charge_atom
                        results_dict[f'ESP {charge_type}'] = ESP_all

                    # Add decay analysis if requested
                    if include_decay:
                        try:
                            [sorted_distances, sorted_esps, cum_esps, sorted_cum_idx,
                             sorted_cum_chg, sorted_atomTypes] = self.esp_bydistance(
                                path_to_xyz, atom_idx, file_path_monopole
                            )
                            results_dict['Sorted Distances'] = sorted_distances
                            results_dict[f'Sorted ESP {charge_type}'] = sorted_esps
                            results_dict[f'Cumulative ESP {charge_type}'] = cum_esps
                            results_dict[f'Dist Sorted Idxs {charge_type}'] = sorted_cum_idx
                            results_dict[f'Dist Sorted Partial Charges {charge_type}'] = sorted_cum_chg
                            results_dict[f'Dist Sorted Atom Types {charge_type}'] = sorted_atomTypes
                        except Exception as e:
                            print(f'Warning: Decay analysis failed for {charge_type}: {e}')

                    # Add coordination shell analysis if requested
                    if include_coord_shells:
                        try:
                            ESP_fcoord = self.esp_first_coord(atom_idx, file_path_monopole, path_to_xyz)
                            ESP_scoord = self.esp_second_coord(atom_idx, file_path_monopole, path_to_xyz)
                            results_dict[f'{charge_type} ESP First Coor Shell (kcal/mol)'] = ESP_fcoord
                            results_dict[f'{charge_type} ESP Second Coor Shell (kcal/mol)'] = ESP_scoord
                        except Exception as e:
                            print(f'Warning: Coordination shell analysis failed for {charge_type}: {e}')

                except Exception as e:
                    logging.exception(f'Exception during ESP calculation for {charge_type}')
                    continue

            allspeciesdict.append(results_dict)

        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(f"{ESPdata_filename}.csv")
        return df

    def getEfield(self, charge_types, Efielddata_filename, multiwfn_module, multiwfn_path, atmrad_path,
                  multipole_bool=False, input_bond_indices=[], auto_find_bonds=False,
                  decompose_atomwise=False, visualize=None, dielectric=1):
        """Unified E-field calculation method supporting multiple modes.

        This consolidated method replaces getEFieldMultipole(), getEfield_acrossBond(),
        and getEfield_decomposable() to eliminate code duplication while preserving all functionality.

        Parameters
        ----------
        charge_types : str or list of str
            Charge partitioning scheme. Options: 'Hirshfeld', 'Becke', 'Hirshfeld_I', etc.
        Efielddata_filename : str
            Output CSV filename (without extension).
        multiwfn_module : str
            Module name containing Multiwfn executable.
        multiwfn_path : str
            Path to Multiwfn executable.
        atmrad_path : str
            Path to atmrad executable.
        multipole_bool : bool, optional
            If True, use multipole expansion; if False, use monopole charges (default: True).
        input_bond_indices : list, optional
            Bond indices as list of tuples [(atomA, atomB), ...].
            If empty and auto_find_bonds=False, uses metal_idx from lst_of_tmcm_idx.
            If empty and auto_find_bonds=True, automatically finds bonded atoms.
        auto_find_bonds : bool, optional
            If True and input_bond_indices is empty, automatically find bonds to adjacent atoms
            using getBondedAtoms() for each metal center (default: False).
        decompose_atomwise : bool, optional
            If True, compute atom-wise E-field decomposition (default: False).
        visualize : bool or None, optional
            Override config visualization settings. None uses config defaults (default: None).
        dielectric : float, optional
            Dielectric constant (default: 1).

        Returns
        -------
        pd.DataFrame
            DataFrame with E-field results saved to CSV.

        Examples
        --------
        >>> # Basic multipole mode with specified bonds
        >>> es.getEfield('Hirshfeld_I', 'efield', ..., input_bond_indices=[(0,1), (0,2)])

        >>> # Auto-find bonds to adjacent atoms
        >>> es.getEfield('Hirshfeld_I', 'efield', ..., auto_find_bonds=True)

        >>> # Monopole mode with atom-wise decomposition
        >>> es.getEfield('Hirshfeld', 'efield', ..., multipole_bool=False, decompose_atomwise=True)
        """
        # Ensure charge_types is a string (E-field works with single scheme)
        if isinstance(charge_types, list):
            if len(charge_types) > 1:
                print(f"Warning: Multiple charge types provided. Using first: {charge_types[0]}")
            charge_type = charge_types[0]
        else:
            charge_type = charge_types

        self.config['dielectric'] = dielectric
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders
        final_structure_file = self.config['xyzfilename']

        # Handle visualization settings
        if visualize is None:
            viz_efield = self.config.get('visualize_efield', False)
            viz_charges = self.config.get('visualize_charges', False)
            viz_per_bond = self.config.get('visualize_per_bond', False)
        else:
            viz_efield = visualize
            viz_charges = visualize
            viz_per_bond = False

        owd = os.getcwd()
        allspeciesdict = []
        counter = 0

        for f in list_of_file:
            try:
                results_dict = {}
                file_path_xyz = f"{owd}/{f + folder_to_molden}{final_structure_file}"
                total_lines = Electrostatics.mapcount(file_path_xyz)
                init_all_lines = range(0, total_lines - 2)

                # Handle excludeAtomfromEcalc - support both flat list and nested list
                exclude_atoms = self.config.get('excludeAtomfromEcalc', [])
                if exclude_atoms and isinstance(exclude_atoms[0], (list, tuple, np.ndarray)):
                    # Nested structure: different exclusions per folder (list of lists/arrays)
                    exclude_list = list(exclude_atoms[counter])
                else:
                    # Flat list/array: same exclusions for all folders
                    exclude_list = list(exclude_atoms) if hasattr(exclude_atoms, '__iter__') else []
                all_lines = [x for x in init_all_lines if x not in exclude_list]

                # Determine bond indices
                if auto_find_bonds or (not input_bond_indices):
                    # Auto-find bonded atoms (requires metal indices)
                    if not self.lst_of_tmcm_idx:
                        raise ValueError(
                            "Auto-finding bonds requires metal atom indices. "
                            "Please provide lst_of_tmcm_idx when initializing the Electrostatics object, "
                            "or specify bond indices explicitly via input_bond_indices parameter."
                        )
                    atom_idx = metal_idxs[counter]
                    bonded_atoms_list = self.getBondedAtoms(file_path_xyz, atom_idx)
                    bond_indices_to_use = [(atom_idx, bonded_idx) for bonded_idx in bonded_atoms_list]
                    print(f"Auto-found {len(bond_indices_to_use)} bonds for atom {atom_idx}: {bond_indices_to_use}")
                elif counter < len(input_bond_indices):
                    bond_indices_to_use = input_bond_indices[counter] if isinstance(input_bond_indices[counter], list) else [input_bond_indices[counter]]
                else:
                    bond_indices_to_use = []

                if not bond_indices_to_use:
                    print(f"Warning: No bond indices for {f}, skipping")
                    counter += 1
                    continue

                # Partition charges
                comp_cost = self.multiwfn.partitionCharge(
                    multipole_bool, f, folder_to_molden, multiwfn_path, atmrad_path, charge_type, owd
                )

                if comp_cost == -1:
                    print(f"Warning: Charge calculation failed for {f}")
                    counter += 1
                    continue

                file_path_multipole = f"{owd}/{f + folder_to_molden}Multipole{charge_type}.txt"
                file_path_charges = f"{owd}/{f + folder_to_molden}Charges{charge_type}.txt"

                # Handle point charges
                df_ptchg = None
                if self.config.get('includePtChgs', False):
                    ptchg_filename = self.config.get('ptChgfp', '')
                    if not ptchg_filename or ptchg_filename.strip() == '':
                        raise ValueError(
                            "Point charge file path is not set. Please call set_ptChgfile() "
                            "with the filename before using includePtChgs=True"
                        )
                    full_ptchg_fp = os.path.join(owd, f, ptchg_filename)
                    df_ptchg = self.getPtChgs(full_ptchg_fp)

                # Select appropriate file path
                path_to_pol = file_path_multipole if multipole_bool else file_path_charges
                temp_xyz = file_path_xyz if multipole_bool else ''

                # Calculate E-field with atom-wise decomposition
                [proj_Efields, bondedAs, bonded_idx, bond_lens, E_proj_atomwise,
                 E_proj_atomwise_list] = self.E_proj_bondIndices_atomwise(
                    bond_indices_to_use, file_path_xyz, path_to_pol, all_lines,
                    multipole_bool, df_ptchg=df_ptchg
                )

                # Calculate bond dipoles for each bond
                bond_dipoles = self.calc_bond_dipoles(
                    bond_indices_to_use, file_path_xyz, path_to_pol, multipole_bool
                )

                # Get charge data for visualization
                df = Electrostatics.charge_atoms(self, path_to_pol, temp_xyz)

                # Visualization: E-field atom-wise contributions
                if viz_efield and E_proj_atomwise is not None:
                    # Create full-length array including zeros for excluded atoms
                    total_atoms = total_lines - 2
                    # E_atomwise contains QM contributions followed by MM contributions
                    # Only use QM contributions (length = len(all_lines)) for visualization
                    n_qm_atoms = len(all_lines)

                    # Always create separate PDB files for each bond when multiple bonds exist
                    if len(E_proj_atomwise_list) > 1:
                        # Multiple bonds: create one PDB per bond
                        for bond_idx, (E_atomwise, bond_pair) in enumerate(zip(E_proj_atomwise_list, bonded_idx)):
                            # Extract only QM contributions (first n_qm_atoms elements)
                            E_atomwise_qm = E_atomwise[:n_qm_atoms]
                            # Expand E_atomwise to include all atoms (zeros for excluded)
                            E_atomwise_full = np.zeros(total_atoms)
                            E_atomwise_full[all_lines] = E_atomwise_qm
                            pdbName = f'Efield_bond{bond_idx}_{bond_pair[0]}-{bond_pair[1]}_{f.rstrip("/").replace("/", "_")}_{charge_type}_.pdb'
                            Visualize(file_path_xyz).makePDB(path_to_pol, E_atomwise_full, pdbName)
                    elif len(E_proj_atomwise_list) == 1:
                        # Single bond: create one PDB with bond info in name
                        E_atomwise_qm = E_proj_atomwise_list[0][:n_qm_atoms]
                        E_atomwise_full = np.zeros(total_atoms)
                        E_atomwise_full[all_lines] = E_atomwise_qm
                        bond_pair = bonded_idx[0]
                        pdbName = f'Efield_bond0_{bond_pair[0]}-{bond_pair[1]}_{f.rstrip("/").replace("/", "_")}_{charge_type}_.pdb'
                        Visualize(file_path_xyz).makePDB(path_to_pol, E_atomwise_full, pdbName)
                    else:
                        # Fallback: use the combined E_proj_atomwise
                        E_proj_atomwise_qm = E_proj_atomwise[:n_qm_atoms]
                        E_proj_atomwise_full = np.zeros(total_atoms)
                        E_proj_atomwise_full[all_lines] = E_proj_atomwise_qm
                        pdbName = f'Efield_cont_{f.rstrip("/").replace("/", "_")}_{charge_type}_.pdb'
                        Visualize(file_path_xyz).makePDB(path_to_pol, E_proj_atomwise_full, pdbName)

                # Visualization: partial charges
                if viz_charges:
                    charge_b_col = df['charge']
                    pdbName = f'Charges_{f.rstrip("/").replace("/", "_")}_{charge_type}_.pdb'
                    Visualize(file_path_xyz).makePDB(path_to_pol, charge_b_col, pdbName)

                # Store results
                results_dict['Max Eproj'] = max(abs(np.array(proj_Efields)))
                results_dict['Projected_Efields V/Angstrom'] = proj_Efields
                results_dict['Bonded Atoms'] = bondedAs
                results_dict['Bonded Indices'] = bonded_idx
                results_dict['Bond Lengths'] = bond_lens
                results_dict['Bond Dipole Vectors (Debye)'] = [bd['bond_dipole_vec'] for bd in bond_dipoles]
                results_dict['Bond Dipole Magnitudes (Debye)'] = [bd['bond_dipole_mag'] for bd in bond_dipoles]
                results_dict['Bond Charges'] = [bd['charges'] for bd in bond_dipoles]
                results_dict['Comp Cost'] = comp_cost
                results_dict['Folder'] = f

                allspeciesdict.append(results_dict)
                counter += 1

            except Exception as e:
                print(f"Error processing {f}: {e}")
                logging.exception(f'Exception during E-field calculation')
                counter += 1
                continue

        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(f"{Efielddata_filename}.csv")
        return df



    def get_residues(self, auto_gen, solute_solvent_dict):
        res_dict = {}
        if auto_gen:
            print(f'I am automatically determining the residues and then assigning to a dict')
        else:
            solute_idxs = solute_solvent_dict['Solute Indices']
            res_dict['Solute'] = solute_idxs
            solvent_idxs = solute_solvent_dict['Solvent Indices']
            num_atoms_solvent = solute_solvent_dict['Num Atoms in Solvent']
            num_solvent_res = int(len(solvent_idxs)/num_atoms_solvent)
            for solv_num in range(0, num_solvent_res):
                res_dict[f'Solvent_{solv_num}'] = solvent_idxs[solv_num*num_atoms_solvent:(solv_num+1)*(num_atoms_solvent)]
        return res_dict



    def partitionCharge(self, multipole_bool, f, folder_to_molden, multiwfn_path, atmrad_path, charge_type,  owd):
        '''
        Partition electron density using Multiwfn to generate partial charges or multipole moments.

        This is the centralized interface for all charge partitioning schemes (Hirshfeld, CHELPG, etc.).
        Checks if calculation was previously completed and reuses results unless rerun=True.

        Parameters:
        -----------
        multipole_bool : bool
            True to compute multipole moments, False for monopoles only
        f : str
            Folder name containing the calculation
        folder_to_molden : str
            Path from folder to molden file location
        multiwfn_path : str
            Path to Multiwfn executable
        atmrad_path : str
            Path to atmrad directory (for Hirshfeld-I calculations)
        charge_type : str
            Partitioning scheme (e.g., 'Hirshfeld', 'CHELPG', 'Hirshfeld_I')
        owd : str
            Original working directory

        Returns:
        --------
        comp_cost : float
            Computation time in seconds, 0 if calculation was previously completed (skipped),
            or -1 if calculation failed
        '''
        molden_filename = self.config['molden_filename']
        final_structure_file = self.config['xyzfilename']
        print(f"Do we need to run calcs?: {self.config['rerun']}")
        comp_cost = 0  # Default to 0 for "already exists" case
        num_atoms = 0
        need_to_run_calculation = True
        os.chdir(owd)
        os.chdir(f + folder_to_molden)
        #subprocess.call(multiwfn_module, shell=True)
        file_path_multipole = f"{os.getcwd()}/{self.config['chgprefix']}Multipole{charge_type}.txt"
        file_path_monopole = f"{os.getcwd()}/{self.config['chgprefix']}Charges{charge_type}.txt"
        file_path_xyz = f"{os.getcwd()}/{final_structure_file}"
        #file_path_xyz = f"{os.getcwd()}/{f + folder_to_molden}{final_structure_file}"

        #check if previous calculations fully converged for desired multipole/charge scheme 
        if multipole_bool:
            if os.path.exists(file_path_multipole):
                with open(file_path_multipole, 'r') as file:
                    contents = file.read()
                    if "Calculation took up" in contents:
                        need_to_run_calculation = False
        else:
            if os.path.exists(file_path_monopole):
                need_to_run_calculation = False
                #Dynamically get path to package settings.ini file
                #path_to_ini_file = str(ini_path)
            

        #If you need to run the calculations, urn either the multipole or the monopole calculation!
        if need_to_run_calculation or self.config['rerun']:
            print(f'Running Calculation')
            try: 
                start = time.time()
                if multipole_bool:
                    chg_prefix, _ = os.path.splitext(molden_filename)
                    #Dynamically get path to package settings.ini file
                    with resources.path('pyef.resources', 'settings.ini') as ini_path:
                        path_to_init_file = str(ini_path)
                        command = f"{multiwfn_path} {molden_filename} -set {path_to_init_file}"

                    multiwfn_commands = ['15', '-1'] + self.multiwfn.dict_of_multipole[charge_type] + ['0', 'q']
                    num_atoms = Electrostatics.mapcount(final_structure_file) - 2
                    if charge_type == 'Hirshfeld_I':
                        #get the number of basis functions
                        num_basis = MoldenObject(file_path_xyz, molden_filename).countBasis()
                        atmrad_src = atmrad_path
                        copy_tree(atmrad_src, os.getcwd() + '/atmrad/')
                        print(f'Current num of basis is: {num_basis}')
                        print(f'The current max num is: {self.config["maxIHirshFuzzyBasis"]}')
                        if  num_basis > self.config['maxIHirshFuzzyBasis']:
                            print(f'Number of basis functions: {num_basis}')
                            multiwfn_commands = ['15', '-1'] + ['4', '-2', '1', '2'] + ['0', 'q']
                            print(f'I-Hirshfeld command should be low memory and slow to accomodate large system')

                    # Use centralized Multiwfn runner
                    self.multiwfn.run_multiwfn(
                        command=command,
                        input_commands=multiwfn_commands,
                        output_file=file_path_multipole,
                        description=f"Multipole {charge_type} calculation for {f}"
                    )

                else:
                    command = f"{multiwfn_path} {molden_filename}"
                    chg_prefix,  _ = os.path.splitext(molden_filename)
                    calc_command = self.multiwfn.dict_of_calcs[charge_type]
                    commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                    if charge_type == 'CHELPG':
                        commands = ['7', calc_command, '1','\n', 'y', '0', 'q']
                    elif charge_type == 'Hirshfeld_I':
                        num_basis = MoldenObject(file_path_xyz, molden_filename).countBasis()
                        print(f'Number of basis functions: {num_basis}')
                        if num_basis > self.config['maxIHirshBasis']:
                            commands = ['7', '15', '-2', '1', '\n', 'y', '0', 'q']
                        atmrad_src = atmrad_path
                        copy_tree(atmrad_src, os.getcwd() + '/atmrad/')
                        #if too many atoms will need to change calc to run with reasonable memory

                    # Use centralized Multiwfn runner
                    self.multiwfn.run_multiwfn(
                        command=command,
                        input_commands=commands,
                        description=f"Monopole {charge_type} calculation for {f}"
                    )
                    os.rename(f'{chg_prefix}.chg', file_path_monopole)
                end = time.time()
                comp_cost = end - start
            except Exception as e:
                print(f"\n{'='*60}")
                print(f"ERROR in partitionCharge for {f}")
                print(f"Charge type: {charge_type}, Multipole: {multipole_bool}")
                print(f"{'='*60}")
                print(f"Error: {str(e)}")
                print(f"\nFull traceback:")
                print(traceback.format_exc())
                print(f"{'='*60}\n")
                # Issue could be from lost memory or Multiwfn failure
                os.chdir(owd)
                return -1  # Return -1 to indicate failure instead of comp_cost
        os.chdir(owd)
        return comp_cost

    def get_residueDipoles(self, charge_type,  multiwfn_module, multiwfn_path, atmrad_path, multipole_bool, solute_indices = [], num_atoms_solvent=0):
        '''
        Function computes a partial charges of all atoms in one system... if the desired file already exists the just return it!!
        Inputs:
        ------- 
        charge_type: string, possible choices include: 'Hirshfeld', 'Voronoi', 'Mulliken',  'Lowdin', 'SCPA', 'Becke', 'ADCH', 'CHELPG', 'MK', 'AIM', 'Hirshfeld_I', 'CM5', 'EEM', 'RESP', 'PEOE'}
                    Also accepts multipole options including: 'Hirshfeld': ['3', '2'], 'Hirshfeld_I': ['4', '1', '2'],   'Becke': ['1', '2']
        '''
        lst_dicts = []
        need_to_run_calculation = True


        # Access Class Variables
        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        owd = os.getcwd() # Old working directory
        molden_filename = self.config['molden_filename']
        final_structure_file = self.config['xyzfilename']
        frame_num_lst = []


        for f in list_of_folders:
            try:
                comp_cost = self.multiwfn.partitionCharge(multipole_bool, f, folder_to_molden, multiwfn_path, atmrad_path, charge_type, owd)
                #If the calculation is not successful, continue
                if comp_cost == -1:
                    print(f"Warning: Charge calculation failed for {f}, skipping")
                    continue
                file_path_multipole = f"{os.getcwd()}/{f + folder_to_molden}Multipole{charge_type}.txt"
                file_path_monopole = f"{os.getcwd()}/{f + folder_to_molden}Charges{charge_type}.txt"
                file_path_xyz = f"{os.getcwd()}/{f + folder_to_molden}{final_structure_file}"
                    
                if multipole_bool:
                    xyz_fp =  file_path_xyz
                    chg_fp = file_path_multipole
                    df_charge_atoms = self.charge_atoms(chg_fp, xyz_fp)
                else:
                    chg_fp = file_path_monopole
                    xyz_fp = ''
                    df_charge_atoms = self.charge_atoms(chg_fp, xyz_fp)
                solute_solvent_dict = {}
                solute_solvent_dict['Solute Indices'] = solute_indices
                solute_solvent_dict['Num Atoms in Solvent'] = num_atoms_solvent
                total_atoms = Electrostatics.mapcount(file_path_xyz) -2
                all_indices = list(np.arange(0, total_atoms))
                solute_solvent_dict['Solvent Indices'] = [x for x in all_indices if x not in solute_indices]
                res_dict = self.get_residues(False, solute_solvent_dict)
                df_dip = self.getdipole_residues(res_dict, df_charge_atoms) 
                #re-order res by closest to the solute
                all_centroids = np.vstack(df_dip['centroid'].values)
                solute_name = 'Solute'
                ref_coord = df_dip.loc[df_dip['Label'] == solute_name, 'centroid']
                ref_coord_arr = ref_coord.values[0]
                # Compute distances to reference
                distances = np.linalg.norm(all_centroids - ref_coord_arr, axis=1)
                # Get sort order (indices of sorted distances)
                df_dip['Distance from solute'] = distances
                init_dipole = df_dip['Dipole']
            
                df_dip.sort_values(by='Distance from solute', ascending=True)
                df_dip.to_csv(f'frm{f}_chg{charge_type}_DistDependentDipoles.csv')
                sorted_dips = df_dip['Dipole']
                avg_dips = np.average(sorted_dips[1:])
            
                first_5A = df_dip.loc[(df_dip['Distance from solute'] < 5), 'Dipole'][1:]
                v5to7A = df_dip.loc[(df_dip['Distance from solute'] < 7) & (df_dip['Distance from solute'] > 5), 'Dipole']
                v7to11A = df_dip.loc[(df_dip['Distance from solute'] < 11) & (df_dip['Distance from solute'] > 7), 'Dipole']
                above11A = df_dip.loc[(df_dip['Distance from solute'] > 11), 'Dipole']

                #return a list of the solvents

                dip_dict = {'DipoleSolute': sorted_dips[0], 'AvgDipSolv': avg_dips, 'Avg5Ang': np.average(first_5A), 'Avg5to7Ang': np.average(v5to7A) , 'Avg7to11Ang': np.average(v7to11A), 'Avgabove11Ang': np.average(above11A), 'CompCost': comp_cost, 'Num5Avg':len(first_5A), 'Num5to7':len(v5to7A), 'Num7to11A': len(v7to11A), 'Numabove11A': len(above11A), 'Total_atoms': total_atoms, 'FrameNum': f }
                lst_dicts.append(dip_dict)
            except Exception as e:
                print(e)
                print(f'Error Traceback: {traceback.print_exc()}')
                continue

        all_file_df = pd.DataFrame(lst_dicts)
        all_file_df.to_csv('Dipoles.csv')
        return all_file_df






    def getpartialchgs(self, charge_types, lst_atom_idxs, partial_chg_filename, multiwfn_path, multiwfn_module, atmrad_path):
        '''
        Function computes partial charges on a select set of atoms using the charge scheme specified in charge types. Note atom indices will be carried over between csvs
        
    Inputs:
        --
        charge_types: list(str)
            list of strings
        lst_of_atom_idxs: list(int)
            list of integers denoting atom indices (0 indexed!)
        partial_chg_filename: string
            Name of the output file name    
        multiwfn_path: string
            Path to the multiwfn executable
        multiwfn_module: string
            Name of the module that contains the multiwfn executable
        atmrad_path: string
            Path to the atmrad executable
        
       Outputs:
        --
        df: Pandas Dataframe with partial charge info
        
        
        Notes
        -----
        Will Create a csv file entitled partial_chg_filename.csv with partial charge info
        '''
       # Access Class Variables
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # Old working directory
        allspeciesdict = []
        counter = 0  # Iterator to account for atomic indices of interest
        for f in list_of_file:
            print('-----------------' + str(f) + '------------------')
            counter = counter + 1
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            subprocess.call(multiwfn_module, shell=True)
            command_A = f"{multiwfn_path} {self.config['molden_filename']}"
            results_dir = os.getcwd() + '/'

            results_dict = {}
            results_dict['Name'] = f
           

            for key in charge_types:
                print('Partial Charge Scheme:' + str(key))
                try:
                    full_file_path = os.getcwd() +'/final_optim_' +key+'.txt'
                    path_to_xyz = os.getcwd() + '/' + self.config['xyzfilename']
                    if key == "Hirshfeld_I":
                        atmrad_src = atmrad_path
                        copy_tree(atmrad_src, results_dir + 'atmrad/')
                    try:
                        for atom_idx in lst_atom_idxs:
                            [atom_type, partial_charge_atom] = Electrostatics.getAtomInfo(full_file_path, atom_idx)
                            results_dict[f'{key} Charge {atom_idx} {atom_type}'] = partial_charge_atom
                    except Exception as e:
                        print('The Exception is: ' + str(e))
                        print(traceback.format_exc())
                        print('Error when trying to access electrostatic information: Attemtping to re-compute partial charges of type: ' + str(key))

                        # Re-run multiwfn computation of partial charge
                        proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                        calc_command = self.multiwfn.dict_of_calcs[key]
                        commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                        if key == 'CHELPG':
                            commands = ['7', calc_command, '1','\n', 'y', '0', 'q']
                        output = proc.communicate("\n".join(commands).encode())
                        new_name = 'final_optim_' +key+'.txt'
                        os.rename('final_optim.chg', new_name)


                        for atom_idx in lst_atom_idxs:
                            [atom_type, partial_charge_atom] = Electrostatics.getAtomInfo(full_file_path, atom_idx)
                            results_dict[f'{key} Charge {atom_idx} {atom_type}'] = partial_charge_atom

                except Exception as e:
                    logging.exception('An Exception was thrown')
                    continue
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(partial_chg_filename +'.csv')
        return df

    def compute_dipole(res_coords, res_chgs, nuc_chgs):
        '''
           res_coords in units of angstroms, res_chgs are in a.u. units
           returns dipole in units of Debye
        '''
        convert_Debye = 4.80320427

        # Center coordinates at geometric center (simple average)
        # This makes the dipole moment origin-independent and consistent with line integral formalism
        center = np.average(res_coords, axis=0)
        coords_centered = res_coords - center
        dipole_vector = convert_Debye*np.sum(coords_centered*res_chgs[:, np.newaxis], axis=0)  #shape of (3,)
        dipole_magnitude = np.linalg.norm(dipole_vector)
        return dipole_vector, dipole_magnitude


    def getdipole_residues(self, res_dict, df):
        ''' res_dict: dictionary with strings mapped to list of lists(int)   
          strings denote name of residues which are mapped to the associated atom indices in the lists(0 indexed!)
        '''
        dipole_vectors = []
        dipole_magnitudes = []
        res_labels = []
        centroid_lst = []

        arr_coords = df[['x', 'y', 'z']].to_numpy()
        atoms = df['Atom']
        arr_chgs = np.array(df['charge'])
        for res_name in res_dict.keys():
            res_labels.append(res_name)
            res_indices = res_dict[res_name]
            res_coords = arr_coords[res_indices]
            res_centroid = res_coords.mean(axis=0)
            res_chgs = arr_chgs[res_indices]
            atom_vals = [atoms[i] for i in res_indices]
            nuc_chgs = [self.amassdict[atom_val][1] for atom_val in atom_vals]
            dip_vec, dip_mag = Electrostatics.compute_dipole(res_coords, res_chgs, nuc_chgs)
            dipole_vectors.append(dip_vec)
            dipole_magnitudes.append(dip_mag)
            centroid_lst.append(res_centroid)
        
        df = pd.DataFrame({
           'Label': res_labels,
           'Dipole_vec': dipole_vectors,
           'Dipole' : dipole_magnitudes,
           'centroid' : centroid_lst
        })
        return df



    def getcharge_residues(self, charge_types, res_dict, partial_chg_filename, multiwfn_path, multiwfn_module, atmrad_path, multipole_mode=True, polarization_scheme='Hirshfeld_I'):
        '''Function computes partial charges on a select set of atoms using the charge scheme specified in charge types. Note atom indices will be carried over between csvs

        Inputs:
        -------
        charge_types: list(str)
            list of strings
        res_dict: dictionary with strings mapped to list of lists(int)   
            strings denote name of residues which are mapped to the associated atom indices in the lists(0 indexed!)
        partial_chg_filename: string
            Name of the output file name
        multiwfn_path: string
            Path to the multiwfn executable
        multiwfn_module: string
            Name of the module that contains the multiwfn executable
        atmrad_path: string
            Path to the atmrad executable
        multipole_mode: boolean
            If True, will use multipole mode to compute partial charges, if False will use the standard method

         Outputs:
        -------
        df: Pandas Dataframe with partial charge info
        
        Notes
        -----
        Will Create a csv file entitled partial_chg_filename.csv with partial charge info
        '''
       # Access Class Variables
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # Old working directory
        allspeciesdict = []
        counter = 0  # Iterator to account for atomic indices of interest

        if multipole_mode:
            final_structure_file = self.config['xyzfilename']
            polarization_file = "Multipole" + polarization_scheme + ".txt"
            molden_filename = self.config['molden_filename']

            file_idx = 0
            for f in list_of_file:
                print('-----------------' + str(f) + '------------------')
                counter = counter + 1
                os.chdir(owd)
                os.chdir(f + folder_to_molden)
                results_dir = os.getcwd() + '/'
                # Check if the atomic polarizations have been computed
                path_to_pol = os.path.join(os.getcwd(), polarization_file)
                backup_path_to_pol = os.path.join(os.getcwd(), "final_optim_polarization.txt")
                #Initially we assume that we do need to run the calculation
                need_to_run_calculation = True
                if os.path.exists(path_to_pol):
                    with open(path_to_pol, 'r') as file:
                        contents = file.read()
                        if "Calculation took up" in contents:
                            print(f"   > Atomic Multipoles already calculated here: {f}!!")
                            need_to_run_calculation = False
                elif os.path.exists(backup_path_to_pol):    
                    path_to_pol = backup_path_to_pol
                    with open(path_to_pol, 'r') as file:
                        contents = file.read()
                        if "Calculation took up" in contents:
                            print(f"   > Atomic Multipoles already calculated here: {f}!!")
                            need_to_run_calculation = False
                if need_to_run_calculation:
                    # Dynamically get path to package settings.ini file
                    with resources.path('pyef.resources', 'settings.ini') as ini_path:
                        path_to_ini_file = str(ini_path)
                        Command_Polarization = f"{multiwfn_path} {molden_filename} -set {path_to_ini_file} > {polarization_file}"
                    subprocess.call(multiwfn_module, shell=True)
                    command_A = f"{multiwfn_path} final_optim.molden"
                    print('Atomic Multipole Calculation initialized')
                    # Now Run the calculation for atomic dipole and quadrupole moment
                    atmrad_src = atmrad_path
                    copy_tree(atmrad_src, os.getcwd() + '/atmrad/')
                    print(f"   > Submitting Multiwfn job using: {Command_Polarization}")
                    proc = subprocess.Popen(Command_Polarization, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                    polarization_commands = ['15', '-1'] + self.multiwfn.dict_of_multipole[polarization_scheme] + ['0', 'q']
                    proc.communicate("\n".join(polarization_commands).encode())
              


                results_dict = {}
                results_dict['Name'] = f


                lst_multipole_dict = MultiwfnInterface.getmultipoles(path_to_pol)
                for res_name in res_dict.keys():
                   all_res_idxs = res_dict[res_name]
                   res_idxs = all_res_idxs[file_idx]
                   res_chg = 0
                   for idx in res_idxs:
                       atom_dict = lst_multipole_dict[idx]
                       chg = atom_dict["Atom_Charge"]
                       res_chg += chg
                   results_dict[f'{res_name}_{polarization_scheme}_Charge'] = res_chg
                   print(f'Current charge is: {res_chg}')
                file_idx += 1
                allspeciesdict.append(results_dict)
         
        else:
            for f in list_of_file:
                print('-----------------' + str(f) + '------------------')
                counter = counter + 1
                os.chdir(owd)
                os.chdir(f + folder_to_molden)
                #working in key mode!
                for key in charge_types:
                    try:
                        full_file_path = os.getcwd() +'/final_optim_' +key+'.txt'
                        if key == "Hirshfeld_I":
                            atmrad_src = atmrad_path
                            copy_tree(atmrad_src, results_dir + 'atmrad/')
                        try:
                            for res_name in res_dict.keys():
                                res_indices = res_dict[res_name]
                                [res_charge] = Electrostatics.getAtomsInfo(full_file_path, res_indices)
                                results_dict[f'{res_name} Charge'] = res_charge

                        except Exception as e:
                            print('The Exception is: ' + str(e))
                            print(traceback.format_exc())
                            print('Error when trying to access electrostatic information: Attemtping to re-compute partial charges of type: ' + str(key))

                            # Re-run multiwfn computation of partial charge
                            proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                            calc_command = self.multiwfn.dict_of_calcs[key]
                            commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                            if key == 'CHELPG':
                                commands = ['7', calc_command, '1','\n', 'y', '0', 'q']
                            output = proc.communicate("\n".join(commands).encode())
                            new_name = 'final_optim_' +key+'.txt'
                            os.rename('final_optim.chg', new_name)


                            for res_name in res_dict.keys():
                                all_res_idxs = res_dict[res_name]
                                res_indices= all_res_idxs[file_idx]
                                [res_charge] = Electrostatics.getAtomsInfo(full_file_path, res_indices)
                                results_dict[f'{res_name}_{key}_Charge'] = res_charge
                                file_idx += 1
                    except Exception as e:
                        logging.exception('An Exception was thrown')
                        continue 
                allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(partial_chg_filename +'.csv')
        return df


    def analyze_atomwise_contributions(self, df_atomwise, top_n=10):
        """
        Analyze and summarize atom-wise contributions to electrostatic stabilization.

        Parameters:
        -----------
        df_atomwise : DataFrame
            DataFrame containing atom-wise contribution data (output from getElectrostatic_stabilization
            with decompose_atomwise=True)
        top_n : int, default=10
            Number of top/bottom contributors to display

        Returns:
        --------
        dict containing:
            - 'summary': DataFrame with total contributions per environment atom (sorted)
            - 'top_stabilizing': DataFrame with top stabilizing atoms (most negative contribution)
            - 'top_destabilizing': DataFrame with top destabilizing atoms (most positive contribution)
            - 'statistics': dict with overall statistics
        """
        # Get unique environment atoms and their total contributions
        env_summary = df_atomwise.groupby(['Env_Atom_Index', 'Env_Element']).agg({
            'Total_Energy_Contribution_kcal_mol': 'first',  # All rows for same env atom have same total
            'Energy_Contribution_kcal_mol': 'sum'  # Verify by summing individual contributions
        }).reset_index()

        env_summary = env_summary.sort_values('Total_Energy_Contribution_kcal_mol', ascending=True)

        # Top stabilizing atoms (most negative values = most stabilizing)
        top_stabilizing = env_summary.head(top_n).copy()
        top_stabilizing['Rank'] = range(1, len(top_stabilizing) + 1)

        # Top destabilizing atoms (most positive values = most destabilizing)
        top_destabilizing = env_summary.tail(top_n).iloc[::-1].copy()
        top_destabilizing['Rank'] = range(1, len(top_destabilizing) + 1)

        # Calculate statistics
        statistics = {
            'total_stabilization_energy': env_summary['Total_Energy_Contribution_kcal_mol'].sum(),
            'mean_contribution': env_summary['Total_Energy_Contribution_kcal_mol'].mean(),
            'std_contribution': env_summary['Total_Energy_Contribution_kcal_mol'].std(),
            'num_stabilizing_atoms': (env_summary['Total_Energy_Contribution_kcal_mol'] < 0).sum(),
            'num_destabilizing_atoms': (env_summary['Total_Energy_Contribution_kcal_mol'] > 0).sum(),
            'most_stabilizing_atom': env_summary.iloc[0]['Env_Atom_Index'] if len(env_summary) > 0 else None,
            'most_destabilizing_atom': env_summary.iloc[-1]['Env_Atom_Index'] if len(env_summary) > 0 else None,
        }

        return {
            'summary': env_summary,
            'top_stabilizing': top_stabilizing,
            'top_destabilizing': top_destabilizing,
            'statistics': statistics
        }

    def _compute_atomwise_dipole_contributions(self, atomwise_Efield, df_substrate, df_env, df_all_atoms, one_mol, filename):
        """
        Compute atom-wise energy contributions from environment E-fields interacting with substrate dipoles.

        The interaction energy is: E = -μ·E (negative dot product of dipole moment with electric field)

        Parameters:
        -----------
        atomwise_Efield : dict
            Dictionary with (substrate_idx, env_idx) as keys and E-field vector [Ex, Ey, Ez] (in V/Å) as values
        df_substrate : DataFrame
            DataFrame containing substrate atom information
        df_env : DataFrame
            DataFrame containing environment atom information
        df_all_atoms : DataFrame
            DataFrame containing all atom information (for element names)
        one_mol : float
            Avogadro's number
        filename : str
            Name of the current file being processed

        Returns:
        --------
        DataFrame with columns: FileName, Env_Atom_Index, Env_Element, Total_Energy_Contribution_kcal_mol,
                                Substrate_Atom_Index, Substrate_Element, Energy_Contribution_kcal_mol
        """
        atomwise_data = []

        # Get substrate dipole mapping (dipole moments are in Bohr from multipoles)
        substrate_dipoles = {}
        substrate_elements = {}
        for _, row in df_substrate.iterrows():
            substrate_dipoles[row['Index']] = np.array(row['Dipole_Moment'])  # Expected to be [μx, μy, μz] in Bohr
            if 'Element' in df_all_atoms.columns:
                substrate_elements[row['Index']] = df_all_atoms.loc[df_all_atoms['Index'] == row['Index'], 'Element'].iloc[0]
            else:
                substrate_elements[row['Index']] = 'Unknown'

        # Get environment element mapping
        env_elements = {}
        for _, row in df_env.iterrows():
            if 'Element' in df_all_atoms.columns:
                env_elements[row['Index']] = df_all_atoms.loc[df_all_atoms['Index'] == row['Index'], 'Element'].iloc[0]
            else:
                env_elements[row['Index']] = 'Unknown'

        # Compute energy contribution for each (substrate, environment) pair
        # Energy = -μ·E where μ is dipole moment and E is electric field
        for (sub_idx, env_idx), efield_contribution in atomwise_Efield.items():
            # efield_contribution is [Ex, Ey, Ez] in V/Å
            # substrate dipole is in Bohr, need to convert to appropriate units
            dipole_bohr = substrate_dipoles[sub_idx]

            # Convert dipole from Bohr to C·m
            dipole_Cm = constants.BOHR_TO_M * constants.ELEMENTARY_CHARGE * dipole_bohr

            # Convert E-field from V/Å to V/m
            efield_Vm = efield_contribution / constants.ANGSTROM_TO_M

            # Energy = -μ·E (in Joules)
            energy_contribution_J = -np.dot(dipole_Cm, efield_Vm)

            # Convert to kcal/mol
            energy_contribution_kcal_mol = energy_contribution_J * one_mol / 4184

            atomwise_data.append({
                'FileName': filename,
                'Env_Atom_Index': env_idx,
                'Env_Element': env_elements.get(env_idx, 'Unknown'),
                'Substrate_Atom_Index': sub_idx,
                'Substrate_Element': substrate_elements.get(sub_idx, 'Unknown'),
                'Energy_Contribution_kcal_mol': energy_contribution_kcal_mol,
                'Efield_x_V_per_A': efield_contribution[0],
                'Efield_y_V_per_A': efield_contribution[1],
                'Efield_z_V_per_A': efield_contribution[2]
            })

        df_atomwise = pd.DataFrame(atomwise_data)

        # Add total contribution from each environment atom (summed over all substrate atoms)
        env_totals = df_atomwise.groupby('Env_Atom_Index')['Energy_Contribution_kcal_mol'].sum().reset_index()
        env_totals.columns = ['Env_Atom_Index', 'Total_Energy_Contribution_kcal_mol']

        # Merge back to get total contributions
        df_atomwise = df_atomwise.merge(env_totals, on='Env_Atom_Index', how='left')

        # Sort by total contribution to identify most important atoms
        df_atomwise = df_atomwise.sort_values('Total_Energy_Contribution_kcal_mol', ascending=False)

        return df_atomwise

    def _visualize_atomwise_contributions(self, df_atomwise, file_path_xyz, charge_file, filename, charge_type):
        """
        Create PDB visualization of atom-wise energy contributions.

        Parameters:
        -----------
        df_atomwise : DataFrame
            DataFrame with atom-wise contribution data from _compute_atomwise_contributions
        file_path_xyz : str
            Path to XYZ file
        charge_file : str
            Path to charge/multipole file
        filename : str
            Folder/file name for output
        charge_type : str
            Type of charge calculation used
        """
        # Get total number of atoms from XYZ file
        total_atoms = Electrostatics.mapcount(file_path_xyz) - 2

        # Initialize contribution array (all atoms start at 0)
        b_factor_contributions = np.zeros(total_atoms)

        # Get unique environment atoms and their total contributions
        env_contributions = df_atomwise.groupby('Env_Atom_Index').agg({
            'Total_Energy_Contribution_kcal_mol': 'first'
        }).reset_index()

        # Fill in the environment atom contributions
        for _, row in env_contributions.iterrows():
            atom_idx = int(row['Env_Atom_Index'])
            contribution = row['Total_Energy_Contribution_kcal_mol']
            b_factor_contributions[atom_idx] = contribution

        # Create PDB file with contributions in B-factor column
        pdb_name = f'Estabilization_atomwise_{filename}{charge_type}.pdb'
        Visualize(file_path_xyz).makePDB(charge_file, b_factor_contributions, pdb_name)
        print(f"Created visualization PDB: {pdb_name}")

    def _compute_atomwise_contributions(self, atomwise_ESP, df_substrate, df_env, df_all_atoms, one_mol, filename):
        """
        Compute atom-wise energy contributions from environment atoms to substrate stabilization.

        Parameters:
        -----------
        atomwise_ESP : dict
            Dictionary with (substrate_idx, env_idx) as keys and ESP contribution (in Volts) as values
        df_substrate : DataFrame
            DataFrame containing substrate atom information
        df_env : DataFrame
            DataFrame containing environment atom information
        df_all_atoms : DataFrame
            DataFrame containing all atom information (for element names)
        one_mol : float
            Avogadro's number
        filename : str
            Name of the current file being processed

        Returns:
        --------
        DataFrame with columns: FileName, Env_Atom_Index, Env_Element, Total_Energy_Contribution_kcal_mol,
                                Substrate_Atom_Index, Substrate_Element, Energy_Contribution_kcal_mol
        """
        atomwise_data = []

        # Get substrate charge mapping
        substrate_charges = {}
        substrate_elements = {}
        for _, row in df_substrate.iterrows():
            substrate_charges[row['Index']] = row['Atom_Charge']
            if 'Element' in df_all_atoms.columns:
                substrate_elements[row['Index']] = df_all_atoms.loc[df_all_atoms['Index'] == row['Index'], 'Element'].iloc[0]
            else:
                substrate_elements[row['Index']] = 'Unknown'

        # Get environment element mapping
        env_elements = {}
        for _, row in df_env.iterrows():
            if 'Element' in df_all_atoms.columns:
                env_elements[row['Index']] = df_all_atoms.loc[df_all_atoms['Index'] == row['Index'], 'Element'].iloc[0]
            else:
                env_elements[row['Index']] = 'Unknown'

        # Compute energy contribution for each (substrate, environment) pair
        # Energy = ESP * charge_substrate (in atomic units, then convert to kcal/mol)
        for (sub_idx, env_idx), esp_contribution in atomwise_ESP.items():
            # ESP is in Volts, charge is in elementary charges
            # Energy = ESP (V) * charge (e) = ESP (J/C) * charge (e) * e (C/e) = Energy (J)
            # Convert to kcal/mol: multiply by Avogadro's number and divide by 4184 J/kcal
            energy_contribution_J = esp_contribution * substrate_charges[sub_idx] * constants.ELEMENTARY_CHARGE
            energy_contribution_kcal_mol = energy_contribution_J * one_mol / 4184

            atomwise_data.append({
                'FileName': filename,
                'Env_Atom_Index': env_idx,
                'Env_Element': env_elements.get(env_idx, 'Unknown'),
                'Substrate_Atom_Index': sub_idx,
                'Substrate_Element': substrate_elements.get(sub_idx, 'Unknown'),
                'Energy_Contribution_kcal_mol': energy_contribution_kcal_mol,
                'ESP_Contribution_V': esp_contribution
            })

        df_atomwise = pd.DataFrame(atomwise_data)

        # Add total contribution from each environment atom (summed over all substrate atoms)
        env_totals = df_atomwise.groupby('Env_Atom_Index')['Energy_Contribution_kcal_mol'].sum().reset_index()
        env_totals.columns = ['Env_Atom_Index', 'Total_Energy_Contribution_kcal_mol']

        # Merge back to get total contributions
        df_atomwise = df_atomwise.merge(env_totals, on='Env_Atom_Index', how='left')

        # Sort by total contribution to identify most important atoms
        df_atomwise = df_atomwise.sort_values('Total_Energy_Contribution_kcal_mol', ascending=False)

        return df_atomwise
    

    def _build_interaction_tensor(self, r_i, r_j, multipole_order, dielectric):
        """
        Build the multipole interaction tensor T_ij using AMOEBA formalism.

        This tensor contains derivatives of 1/r and enables computation of all
        multipole-multipole interactions via: E = M_i^T · T_ij · M_j

        Parameters:
        -----------
        r_i, r_j : array-like
            Position vectors of atoms i and j in Angstroms
        multipole_order : int
            Maximum multipole order (1=monopole, 2=+dipole, 3=+quadrupole)
        dielectric : float
            Dielectric constant

        Returns:
        --------
        T : ndarray
            Interaction tensor with dimensions based on multipole_order
            - Order 1: 1×1 (monopole-monopole only)
            - Order 2: 4×4 (monopole + dipole)
            - Order 3: 9×9 (monopole + dipole + quadrupole, traceless)

        Notes:
        ------
        Tensor elements represent spatial derivatives of Coulomb potential:
        - T[0,0]: 1/r (monopole-monopole)
        - T[0,1:4]: ∂(1/r)/∂r_j (monopole-dipole)
        - T[1:4,1:4]: ∂²(1/r)/∂r_i∂r_j (dipole-dipole)
        - Higher orders follow similar pattern
        """
        # Convert to meters and compute distance
        r_vec = (r_j - r_i) * constants.ANGSTROM_TO_M
        r = np.linalg.norm(r_vec)
        r_hat = r_vec / r  # Unit vector from i to j

        # Coulomb prefactor with dielectric
        # Note: charges are dimensionless (in units of e), so we need e² for energy
        k = constants.COULOMB_CONSTANT * (constants.ELEMENTARY_CHARGE ** 2) / dielectric

        if multipole_order == 1:
            # Monopole-monopole only: T = k/r
            T = np.array([[k / r]])

        elif multipole_order == 2:
            # Monopole + dipole interactions (4×4 tensor)
            T = np.zeros((4, 4))

            # T[0,0]: monopole-monopole interaction (q_i × q_j / r)
            T[0, 0] = k / r

            # T[0, 1:4]: monopole-dipole interaction (∂/∂r_j of 1/r)
            # Derivative: ∂(1/r)/∂r_j = -r_hat/r²
            T[0, 1:4] = -k * r_hat / (r**2)

            # T[1:4, 0]: dipole-monopole interaction (∂/∂r_i of 1/r)
            # By symmetry with opposite sign: ∂(1/r)/∂r_i = +r_hat/r²
            T[1:4, 0] = k * r_hat / (r**2)

            # T[1:4, 1:4]: dipole-dipole interaction (∂²/∂r_i∂r_j of 1/r)
            # Formula: (3·r_hat⊗r_hat - I)/r³
            I = np.eye(3)
            T[1:4, 1:4] = k * (3 * np.outer(r_hat, r_hat) - I) / (r**3)

        elif multipole_order == 3:
            # Monopole + dipole + quadrupole interactions (9×9 tensor)
            T = np.zeros((9, 9))

            # Monopole-monopole (same as order 2)
            T[0, 0] = k / r

            # Monopole-dipole (same as order 2)
            T[0, 1:4] = -k * r_hat / (r**2)
            T[1:4, 0] = k * r_hat / (r**2)

            # Dipole-dipole (same as order 2)
            I = np.eye(3)
            T[1:4, 1:4] = k * (3 * np.outer(r_hat, r_hat) - I) / (r**3)

            # Monopole-quadrupole: T[0, 4:9] and T[4:9, 0]
            # Uses traceless quadrupole representation (5 independent components)
            # ∂²(1/r)/∂r_α∂r_β derivatives
            for idx, (alpha, beta) in enumerate([(0,0), (0,1), (0,2), (1,1), (1,2)]):
                # Monopole-quadrupole term
                delta_ab = 1.0 if alpha == beta else 0.0
                T_mq = k * (3*r_hat[alpha]*r_hat[beta] - delta_ab) / (r**3)
                T[0, 4+idx] = T_mq
                T[4+idx, 0] = T_mq

            # Dipole-quadrupole: T[1:4, 4:9] and T[4:9, 1:4]
            # Third derivatives: ∂³(1/r)/∂r_i∂r_α∂r_β
            for i in range(3):
                for idx, (alpha, beta) in enumerate([(0,0), (0,1), (0,2), (1,1), (1,2)]):
                    delta_ia = 1.0 if i == alpha else 0.0
                    delta_ib = 1.0 if i == beta else 0.0
                    delta_ab = 1.0 if alpha == beta else 0.0

                    T_dq = k * (-15*r_hat[i]*r_hat[alpha]*r_hat[beta]
                               + 3*(delta_ia*r_hat[beta] + delta_ib*r_hat[alpha] + delta_ab*r_hat[i])) / (r**4)
                    T[1+i, 4+idx] = T_dq
                    T[4+idx, 1+i] = T_dq

            # Quadrupole-quadrupole: T[4:9, 4:9]
            # Fourth derivatives: ∂⁴(1/r)/∂r_α∂r_β∂r_γ∂r_δ
            for idx1, (alpha, beta) in enumerate([(0,0), (0,1), (0,2), (1,1), (1,2)]):
                for idx2, (gamma, delta) in enumerate([(0,0), (0,1), (0,2), (1,1), (1,2)]):
                    delta_ab = 1.0 if alpha == beta else 0.0
                    delta_ag = 1.0 if alpha == gamma else 0.0
                    delta_ad = 1.0 if alpha == delta else 0.0
                    delta_bg = 1.0 if beta == gamma else 0.0
                    delta_bd = 1.0 if beta == delta else 0.0
                    delta_gd = 1.0 if gamma == delta else 0.0

                    T_qq = k * (105*r_hat[alpha]*r_hat[beta]*r_hat[gamma]*r_hat[delta]
                               - 15*(r_hat[alpha]*r_hat[beta]*delta_gd + r_hat[gamma]*r_hat[delta]*delta_ab
                                    + r_hat[alpha]*r_hat[gamma]*delta_bd + r_hat[alpha]*r_hat[delta]*delta_bg
                                    + r_hat[beta]*r_hat[gamma]*delta_ad + r_hat[beta]*r_hat[delta]*delta_ag)
                               + 3*(delta_ab*delta_gd + delta_ag*delta_bd + delta_ad*delta_bg)) / (r**5)
                    T[4+idx1, 4+idx2] = T_qq
        else:
            raise ValueError(f"multipole_order must be 1, 2, or 3. Got {multipole_order}")

        return T

    def _build_multipole_vector(self, atom_data, multipole_order):
        """
        Build multipole moment vector M for an atom using AMOEBA formalism.

        Parameters:
        -----------
        atom_data : dict or pd.Series
            Atom data containing 'Atom_Charge', 'Dipole_Moment', 'Quadrupole_Moment'
        multipole_order : int
            Maximum multipole order (1=monopole, 2=+dipole, 3=+quadrupole)

        Returns:
        --------
        M : ndarray
            Multipole vector [q, μ_x, μ_y, μ_z, Q_xx, Q_xy, ...] in SI-based units
            - Order 1: [q] (1 element)
            - Order 2: [q, μ_x, μ_y, μ_z] (4 elements)
            - Order 3: [q, μ_x, μ_y, μ_z, Q_xx, Q_xy, Q_xz, Q_yy, Q_yz] (9 elements, traceless)

        Notes:
        ------
        - Charges are dimensionless (in units of elementary charge)
        - Dipoles are converted from Bohr to meters
        - Quadrupoles are converted from Bohr² to meters²
        - Uses traceless quadrupole representation (Q_zz = -Q_xx - Q_yy)
        """
        M = [atom_data['Atom_Charge']]  # Monopole (dimensionless, in units of e)

        if multipole_order >= 2:
            # Add dipole components (convert from atomic units to SI)
            # Multiwfn outputs dipole in Bohr (atomic units)
            # In a.u., dipole has units of e·a₀, so we only convert length
            dipole_bohr = np.array(atom_data['Dipole_Moment'])
            dipole_m = dipole_bohr * constants.BOHR_TO_M  # Now in units of e·meters
            M.extend(dipole_m)

        if multipole_order >= 3:
            # Add quadrupole components (convert from Bohr² to m²)
            # Use traceless representation: 5 independent components
            Q = atom_data['Quadrupole_Moment']
            Q_m2 = Q * (constants.BOHR_TO_M**2)

            # Extract 5 independent components (Q is symmetric and traceless)
            # Components: Q_xx, Q_xy, Q_xz, Q_yy, Q_yz
            # (Q_zz is determined by traceless condition: Q_zz = -Q_xx - Q_yy)
            M.extend([Q_m2[0, 0], Q_m2[0, 1], Q_m2[0, 2], Q_m2[1, 1], Q_m2[1, 2]])

        return np.array(M)

    def getElectrostatic_stabilization(self, multiwfn_path, multiwfn_module, atmrad_path,
                                              substrate_idxs, charge_type='Hirshfeld_I',
                                              name_dataStorage='estaticFile_tensor', env_idxs=None,
                                              decompose_atomwise=False, visualize=None,
                                              multipole_order=2, substrate_multipole_order=None,
                                              env_multipole_order=None):
        """
        Compute electrostatic stabilization using AMOEBA-style multipole tensor formalism.

        This function implements the complete multipole expansion for intermolecular
        electrostatic interactions using the tensor formalism:

            E_total = Σᵢ∈substrate Σⱼ∈environment M_i^T · T_ij · M_j

        where M is the multipole vector and T_ij is the interaction tensor containing
        all derivatives of the Coulomb potential. This naturally includes ALL
        multipole-multipole interaction terms:

        - Order 1 (monopole): charge-charge (q×q)
        - Order 2 (+ dipole): + charge-dipole (q×μ) + dipole-dipole (μ×μ)
        - Order 3 (+ quadrupole): + charge-quadrupole (q×Q) + dipole-quadrupole (μ×Q) + quadrupole-quadrupole (Q×Q)

        This is the formalism used in the AMOEBA polarizable force field and provides
        a complete treatment of distributed multipole interactions.

        Parameters:
        -----------
        multiwfn_path : str
            Path to Multiwfn executable
        multiwfn_module : str
            Module name containing Multiwfn executable
        atmrad_path : str
            Path to atmrad executable
        substrate_idxs : list
            List of substrate atom indices (the molecule being stabilized)
        charge_type : str, optional
            Charge partitioning scheme (default: 'Hirshfeld_I')
        name_dataStorage : str, optional
            Output filename prefix (default: 'estaticFile_tensor')
        env_idxs : list or None, optional
            List of environment atom indices. If None, uses all non-substrate atoms (default: None)
        decompose_atomwise : bool, optional
            If True, compute atom-wise decomposition showing each substrate atom's contribution (default: False)
        visualize : bool or None, optional
            If True, create PDB files with energy contributions in B-factor column.
            None uses config defaults (default: None)
        multipole_order : int, optional
            Order of multipole expansion:
            - 1: monopole only (charge-charge)
            - 2: monopole + dipole (includes charge-charge, charge-dipole, dipole-dipole)
            - 3: monopole + dipole + quadrupole (includes all terms up to Q×Q)
            (default: 2)
        substrate_multipole_order : int or None, optional
            Multipole order for substrate atoms. If None, uses multipole_order.
            Set to 1 for QM/MM where substrate is MM (charges only). (default: None)
        env_multipole_order : int or None, optional
            Multipole order for environment atoms. If None, uses multipole_order.
            Set to 1 for QM/MM where environment is MM (charges only). (default: None)

        Returns:
        --------
        pd.DataFrame or tuple of pd.DataFrame
            If decompose_atomwise=False: returns DataFrame with total stabilization energies
            If decompose_atomwise=True: returns (total_df, atomwise_df) tuple

        Notes:
        ------
        **Comparison with other methods:**
        - `getElectrostatic_stabilization()`: Computes q·V (substrate charge in environment potential)
        - `get_Electrostatic_stabilization_dipole()`: Computes -μ·E (substrate dipole in environment field)
        - This function: Computes DIRECT pairwise multipole-multipole interactions

        **Advantages of tensor formalism:**
        - Complete: automatically includes all multipole-multipole terms
        - Symmetric: properly treats both substrate and environment multipoles
        - Efficient: single matrix multiplication per atom pair
        - Extensible: easy to add higher-order terms

        **Important:**
        - Requires Multiwfn multipole analysis for orders ≥ 2
        - **Automatic fallback**: If multipole files don't exist, automatically uses charges-only files
        - **QM/MM support**: Substrate and environment can have different multipole orders (e.g., QM substrate with multipoles, MM environment with charges only)
        - Results are the electrostatic STABILIZATION energy (environment → substrate)
        - Energy is in kcal/mol
        - Follows AMOEBA force field formalism (Ren & Ponder, J. Phys. Chem. B 2003)

        Examples:
        ---------
        >>> # Compute monopole + dipole + dipole-dipole interactions
        >>> df = estat.get_Electrostatic_stabilization_tensor(
        ...     multiwfn_path, multiwfn_module, atmrad_path,
        ...     substrate_idxs=[[0, 1, 2]],  # 3-atom substrate
        ...     multipole_order=2  # Include dipole-dipole terms
        ... )

        >>> # Include quadrupole terms with atom-wise decomposition
        >>> df, df_atomwise = estat.get_Electrostatic_stabilization_tensor(
        ...     multiwfn_path, multiwfn_module, atmrad_path,
        ...     substrate_idxs=[[0, 1, 2]],
        ...     multipole_order=3,
        ...     decompose_atomwise=True
        ... )

        >>> # QM/MM: QM substrate (with dipoles) interacting with MM environment (charges only)
        >>> df = estat.get_Electrostatic_stabilization_tensor(
        ...     multiwfn_path, multiwfn_module, atmrad_path,
        ...     substrate_idxs=[[0, 1, 2]],  # QM region
        ...     substrate_multipole_order=2,  # QM: charges + dipoles
        ...     env_multipole_order=1         # MM: charges only
        ... )
        """
        dielectric = self.config['dielectric']
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders
        final_structure_file = self.config['xyzfilename']

        # Handle visualization settings
        if visualize is None:
            viz_contributions = self.config.get('visualize_estabilization', False)
        else:
            viz_contributions = visualize

        owd = os.getcwd()
        allspeciesdict = []
        all_atomwise_data = []
        counter = 0
        one_mol = 6.02e23

        # Set default multipole orders for substrate and environment
        if substrate_multipole_order is None:
            substrate_multipole_order = multipole_order
        if env_multipole_order is None:
            env_multipole_order = multipole_order

        # Determine if we need multipole analysis (orders 2 and 3 need dipoles/quadrupoles)
        # We need multipoles if EITHER substrate or environment needs them
        need_multipoles = (substrate_multipole_order >= 2 or env_multipole_order >= 2)

        for f in list_of_file:
            substrate_idx = substrate_idxs[counter]
            results_dict = {}
            file_path_xyz = f"{os.getcwd()}/{f + folder_to_molden}{final_structure_file}"
            total_lines = Electrostatics.mapcount(file_path_xyz)
            init_all_lines = np.arange(0, total_lines - 2)

            if not env_idxs:
                env_idx = [x for x in init_all_lines if x not in substrate_idx]
            else:
                env_idx = env_idxs[counter]

            # Partition charges (with multipole analysis if needed)
            comp_cost = self.multiwfn.partitionCharge(need_multipoles, f, folder_to_molden,
                                            multiwfn_path, atmrad_path, charge_type, owd)

            if comp_cost == -1:
                print(f"Warning: Charge calculation failed for {f}, skipping")
                counter += 1
                continue

            # Get geometry information
            geom = Geometry(file_path_xyz)
            df_geom = geom.getGeomInfo()

            # Load charge/multipole data with automatic fallback
            multipole_name = f"{os.getcwd()}/{f + folder_to_molden}Multipole{charge_type}.txt"
            monopole_name = f"{os.getcwd()}/{f + folder_to_molden}Charges{charge_type}.txt"

            # Try to load multipole file first, fall back to charges if not available
            multipole_available = os.path.exists(multipole_name)

            if need_multipoles and not multipole_available:
                print(f"Warning: Multipole file not found for {f}, falling back to charges-only")
                print(f"  Requested multipole_order={multipole_order}, but only charges available")
                print(f"  Setting both substrate and environment to monopole-only")
                substrate_multipole_order = 1
                env_multipole_order = 1
                need_multipoles = False

            if multipole_available and need_multipoles:
                # Load multipole data
                atomicDict = MultiwfnInterface.getmultipoles(multipole_name)
                df_all = pd.DataFrame(atomicDict)
                df_all['Index'] = df_all['Index'].astype(int) - 1
                multipole_dict = {int(row['Index']): row for _, row in df_all.iterrows()}
            else:
                # Load charges-only data
                df_all = pd.read_csv(monopole_name, sep='\s+',
                                    names=["Element", 'x', 'y', 'z', "Atom_Charge"])
                df_all["Index"] = range(0, len(df_all))
                df_all['Dipole_Moment'] = [[0.0, 0.0, 0.0]] * len(df_all)
                df_all['Quadrupole_Moment'] = [np.zeros((3, 3))] * len(df_all)
                multipole_dict = {row['Index']: row for _, row in df_all.iterrows()}

            # Now handle substrate and environment separately for QM/MM support
            df_substrate = df_all[df_all["Index"].isin(substrate_idx)].copy()
            df_env = df_all[df_all["Index"].isin(env_idx)].copy()

            # If substrate is monopole-only but we loaded multipoles, zero out higher terms
            if substrate_multipole_order == 1 and multipole_available:
                df_substrate['Dipole_Moment'] = [[0.0, 0.0, 0.0]] * len(df_substrate)
                df_substrate['Quadrupole_Moment'] = [np.zeros((3, 3))] * len(df_substrate)

            # If environment is monopole-only but we loaded multipoles, zero out higher terms
            if env_multipole_order == 1 and multipole_available:
                df_env['Dipole_Moment'] = [[0.0, 0.0, 0.0]] * len(df_env)
                df_env['Quadrupole_Moment'] = [np.zeros((3, 3))] * len(df_env)

            # Compute multipole tensor interactions
            total_energy = 0.0
            atomwise_contributions = []

            # Loop over all substrate-environment atom pairs
            for sub_idx in substrate_idx:
                sub_atom_energy = 0.0

                # Get substrate atom data using dictionary lookup for multipoles or direct dataframe for monopoles
                if need_multipoles:
                    if sub_idx not in multipole_dict:
                        print(f"Warning: Substrate atom {sub_idx} not found in multipole data, skipping")
                        continue
                    sub_row = multipole_dict[sub_idx]
                else:
                    sub_row = df_substrate[df_substrate['Index'] == sub_idx].iloc[0]

                # Get substrate atom position
                r_i = np.array([df_geom['X'][sub_idx],
                               df_geom['Y'][sub_idx],
                               df_geom['Z'][sub_idx]])

                # Build substrate multipole vector with substrate-specific order
                M_i = self._build_multipole_vector(sub_row, substrate_multipole_order)

                for env_idx_atom in env_idx:
                    # Get environment atom data using dictionary lookup for multipoles or direct dataframe for monopoles
                    if need_multipoles:
                        if env_idx_atom not in multipole_dict:
                            print(f"Warning: Environment atom {env_idx_atom} not found in multipole data, skipping")
                            continue
                        env_row = multipole_dict[env_idx_atom]
                    else:
                        env_row = df_env[df_env['Index'] == env_idx_atom].iloc[0]

                    # Get environment atom position
                    r_j = np.array([df_geom['X'][env_idx_atom],
                                   df_geom['Y'][env_idx_atom],
                                   df_geom['Z'][env_idx_atom]])

                    # Build environment multipole vector with environment-specific order
                    M_j = self._build_multipole_vector(env_row, env_multipole_order)

                    # Build interaction tensor with max order (tensor must accommodate both)
                    tensor_order = max(substrate_multipole_order, env_multipole_order)
                    T_ij = self._build_interaction_tensor(r_i, r_j, tensor_order, dielectric)

                    # Pad multipole vectors to match tensor dimensions if needed
                    # This handles QM/MM cases where substrate and environment have different orders
                    if len(M_i) < len(M_j):
                        M_i_padded = np.pad(M_i, (0, len(M_j) - len(M_i)))
                        M_j_padded = M_j
                    elif len(M_j) < len(M_i):
                        M_i_padded = M_i
                        M_j_padded = np.pad(M_j, (0, len(M_i) - len(M_j)))
                    else:
                        M_i_padded = M_i
                        M_j_padded = M_j

                    # Compute pairwise interaction: E = M_i^T · T_ij · M_j
                    E_pair_J = M_i_padded.T @ T_ij @ M_j_padded  # Energy in Joules
                    E_pair_kcal_mol = E_pair_J * one_mol / 4184.0  # Convert to kcal/mol

                    sub_atom_energy += E_pair_kcal_mol
                    total_energy += E_pair_kcal_mol

                # Store atom-wise decomposition if requested
                if decompose_atomwise:
                    atom_data = {
                        'FileName': f,
                        'Atom_Index': sub_idx,
                        'Atom_Symbol': df_geom.loc[df_geom.index == sub_idx, 'Atom'].iloc[0],
                        'Charge': sub_row['Atom_Charge'],
                        'Energy_Contribution_kcal_mol': sub_atom_energy
                    }

                    # Add multipole information if substrate has multipoles
                    if substrate_multipole_order >= 2:
                        dipole = np.array(sub_row['Dipole_Moment'])
                        atom_data.update({
                            'Dipole_x_Bohr': dipole[0],
                            'Dipole_y_Bohr': dipole[1],
                            'Dipole_z_Bohr': dipole[2],
                            'Dipole_magnitude_Bohr': np.linalg.norm(dipole)
                        })

                    if substrate_multipole_order >= 3:
                        Q = sub_row['Quadrupole_Moment']
                        atom_data.update({
                            'Quadrupole_xx_Bohr2': Q[0, 0],
                            'Quadrupole_yy_Bohr2': Q[1, 1],
                            'Quadrupole_zz_Bohr2': Q[2, 2],
                            'Quadrupole_xy_Bohr2': Q[0, 1]
                        })

                    atomwise_contributions.append(atom_data)

            # Store results for this file
            results_dict['Total_Energy_kcal_mol'] = total_energy
            results_dict['Substrate_Multipole_Order'] = substrate_multipole_order
            results_dict['Environment_Multipole_Order'] = env_multipole_order
            results_dict['Num_Substrate_Atoms'] = len(substrate_idx)
            results_dict['Num_Environment_Atoms'] = len(env_idx)
            results_dict['FileName'] = f
            allspeciesdict.append(results_dict)

            if decompose_atomwise:
                all_atomwise_data.extend(atomwise_contributions)

                # Visualize atom-wise contributions if requested
                if viz_contributions:
                    df_atomwise_temp = pd.DataFrame(atomwise_contributions)
                    # Create PDB visualization with substrate atom contributions in B-factor
                    try:
                        total_atoms = len(df_geom)
                        b_factors = np.zeros(total_atoms)

                        # Fill in substrate atom contributions
                        for _, row in df_atomwise_temp.iterrows():
                            atom_idx = int(row['Atom_Index'])
                            if 0 <= atom_idx < total_atoms:
                                b_factors[atom_idx] = row['Energy_Contribution_kcal_mol']

                        # Create PDB file
                        pdb_name = f'Estabilization_tensor_{f}_sub{substrate_multipole_order}_env{env_multipole_order}.pdb'
                        if multipole_available and (substrate_multipole_order >= 2 or env_multipole_order >= 2):
                            Visualize(file_path_xyz).makePDB(multipole_name, b_factors, pdb_name)
                        else:
                            Visualize(file_path_xyz).makePDB(monopole_name, b_factors, pdb_name)
                        print(f"    Created visualization PDB: {pdb_name}")
                    except Exception as e:
                        print(f"    Warning: Could not create visualization PDB: {e}")

            print(f"File {f}: Total energy (substrate order {substrate_multipole_order}, env order {env_multipole_order}): {total_energy:.4f} kcal/mol")

            os.chdir(owd)
            counter += 1

        # Create output DataFrames
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(name_dataStorage + '.csv', index=False)

        if decompose_atomwise:
            df_atomwise = pd.DataFrame(all_atomwise_data)
            df_atomwise.to_csv(name_dataStorage + '_atomwise.csv', index=False)
            return df, df_atomwise

        return df