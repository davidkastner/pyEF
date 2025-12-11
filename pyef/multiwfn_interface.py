"""
Multiwfn Interface Module

This module provides a clean interface to the Multiwfn external quantum chemistry
analysis program for charge partitioning and multipole moment calculations.

Classes
-------
MultiwfnInterface
    Main interface class for running Multiwfn calculations

Functions
---------
getmultipoles
    Static function to parse multipole moment output files
"""

import os
import re
import sys
import time
import subprocess
import traceback
import numpy as np
from importlib import resources
from distutils.dir_util import copy_tree


class MoldenObject:
    """Helper class for parsing molden files to count basis functions."""

    def __init__(self, file_path_xyz, molden_filename):
        self.file_path_xyz = file_path_xyz
        self.molden_filename = molden_filename

    def countBasis(self):
        """Count the number of basis functions in the molden file."""
        # This is a simplified version - you may need to import from the actual module
        basis_count = 0
        with open(self.molden_filename, 'r') as f:
            for line in f:
                if '[MO]' in line:
                    break
                if 'Ene=' in line or 'Occup=' in line:
                    basis_count += 1
        return basis_count // 2 if basis_count > 0 else 0


class MultiwfnInterface:
    """
    Interface class for Multiwfn charge partitioning and multipole calculations.

    This class encapsulates all interactions with the Multiwfn external program,
    providing methods for charge partitioning, multipole moment calculations,
    and parsing of Multiwfn output files.

    Attributes
    ----------
    dict_of_calcs : dict
        Mapping of charge scheme names to Multiwfn command codes
    dict_of_multipole : dict
        Mapping of multipole schemes to Multiwfn command sequences
    config : dict
        Configuration dictionary with Multiwfn settings
    """

    # Multiwfn command codes for charge calculation methods
    CHARGE_SCHEMES = {
        'Hirshfeld': '1', 'Voronoi': '2', 'Mulliken': '5',
        'Lowdin': '6', 'SCPA': '7', 'Becke': '10', 'ADCH': '11',
        'CHELPG': '12', 'MK': '13', 'AIM': '14', 'Hirshfeld_I': '15',
        'CM5': '16', 'EEM': '17', 'RESP': '18', 'PEOE': '19'
    }

    # Multiwfn command codes for multipole moment calculations
    MULTIPOLE_SCHEMES = {
        'Hirshfeld': ['3', '2'],
        'Hirshfeld_I': ['4', '1', '2'],
        'Becke': ['1', '2']
    }

    def __init__(self, config=None):
        """
        Initialize the Multiwfn interface.

        Parameters
        ----------
        config : dict, optional
            Configuration dictionary with keys:
            - molden_filename: Name of molden file
            - xyzfilename: Name of XYZ geometry file
            - rerun: Whether to rerun existing calculations
            - hasECP: Whether calculation used ECP
            - maxIHirshFuzzyBasis: Max basis functions for Hirshfeld-I fuzzy
            - maxIHirshBasis: Max basis functions for Hirshfeld-I
        """
        self.config = config if config is not None else {
            'molden_filename': 'final_optim.molden',
            'xyzfilename': 'final_optim.xyz',
            'rerun': False,
            'hasECP': False,
            'maxIHirshFuzzyBasis': 1500,
            'maxIHirshBasis': 2000
        }

        self.dict_of_calcs = self.CHARGE_SCHEMES
        self.dict_of_multipole = self.MULTIPOLE_SCHEMES

    def run_multiwfn(self, command, input_commands, output_file=None,
                     description="Multiwfn calculation", capture_output=True):
        """
        Centralized Multiwfn runner with proper error handling and output display.

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
                    command,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=True,
                    text=True
                )

                # Send input commands
                input_str = '\n'.join(input_commands) + '\n'
                stdout, stderr = proc.communicate(input=input_str)

                # Save to output file if specified
                if output_file:
                    with open(output_file, 'w') as f:
                        f.write(stdout)

                # Display output
                if stdout:
                    print("STDOUT:")
                    print(stdout[:500])  # Print first 500 chars
                    if len(stdout) > 500:
                        print(f"... (truncated {len(stdout) - 500} characters)")

                if stderr:
                    print("\nSTDERR:")
                    print(stderr)

                returncode = proc.returncode
            else:
                # Let output stream to terminal
                proc = subprocess.Popen(
                    command,
                    stdin=subprocess.PIPE,
                    shell=True,
                    text=True
                )
                input_str = '\n'.join(input_commands) + '\n'
                proc.communicate(input=input_str)
                stdout, stderr, returncode = "", "", proc.returncode

            if returncode != 0:
                raise RuntimeError(
                    f"{description} failed with exit code {returncode}\n"
                    f"Command: {full_command}"
                )

            print(f"\n{'='*60}")
            print(f"Completed: {description}")
            print(f"{'='*60}\n")

            return stdout, stderr, returncode

        except Exception as e:
            print(f"\n{'='*60}")
            print(f"ERROR running Multiwfn")
            print(f"Description: {description}")
            print(f"Command: {full_command}")
            print(f"{'='*60}")
            print(f"Error: {str(e)}")
            print(f"{'='*60}\n")
            raise

    def partitionCharge(self, multipole_bool, f,
                       multiwfn_path, atmrad_path, charge_type, owd,
                       molden_filename=None, xyz_filename=None):
        """
        Partition electron density using Multiwfn to generate partial charges or multipole moments.

        This is the centralized interface for all charge partitioning schemes (Hirshfeld, CHELPG, etc.).
        Checks if calculation was previously completed and reuses results unless rerun=True.

        Parameters
        ----------
        multipole_bool : bool
            True to compute multipole moments, False for monopoles only
        f : str
            Complete path to folder containing the calculation
        multiwfn_path : str
            Path to Multiwfn executable
        atmrad_path : str
            Path to atmrad directory (for Hirshfeld-I calculations)
        charge_type : str
            Partitioning scheme (e.g., 'Hirshfeld', 'CHELPG', 'Hirshfeld_I')
        owd : str
            Original working directory
        molden_filename : str, optional
            Name of the molden file (not full path, just filename).
            If None, uses self.config['molden_filename']
        xyz_filename : str, optional
            Name of the xyz file (not full path, just filename).
            If None, uses self.config['xyzfilename']

        Returns
        -------
        float
            Computation time in seconds, 0 if calculation was previously completed (skipped),
            or -1 if calculation failed
        """
        # Handle f being a list (take first element)
        if isinstance(f, list):
            f = f[0]

        # Use provided filenames or fall back to config defaults
        if molden_filename is None:
            molden_filename = self.config['molden_filename']
        if xyz_filename is None:
            final_structure_file = self.config['xyzfilename']
            # Handle xyzfilename being a list (take first element)
            if isinstance(final_structure_file, list):
                final_structure_file = final_structure_file[0]
        else:
            final_structure_file = xyz_filename
        print(f"Do we need to run calcs?: {self.config['rerun']}")
        comp_cost = 0  # Default to 0 for "already exists" case
        num_atoms = 0
        need_to_run_calculation = True

        os.chdir(owd)
        os.chdir(f)

        file_path_multipole = f"{os.getcwd()}/Multipole{charge_type}.txt"
        file_path_monopole = f"{os.getcwd()}/Charges{charge_type}.txt"
        file_path_xyz = f"{os.getcwd()}/{final_structure_file}"

        # Check if previous calculations fully converged for desired multipole/charge scheme
        if multipole_bool:
            if os.path.exists(file_path_multipole):
                with open(file_path_multipole, 'r') as file:
                    contents = file.read()
                    if "Calculation took up" in contents:
                        need_to_run_calculation = False
        else:
            if os.path.exists(file_path_monopole):
                need_to_run_calculation = False

        # If you need to run the calculations, run either the multipole or the monopole calculation!
        if need_to_run_calculation or self.config['rerun']:
            print(f'Running Calculation')
            try:
                start = time.time()
                if multipole_bool:
                    chg_prefix, _ = os.path.splitext(molden_filename)
                    # Dynamically get path to package settings.ini file
                    with resources.path('pyef.resources', 'settings.ini') as ini_path:
                        path_to_init_file = str(ini_path)
                        command = f"{multiwfn_path} {molden_filename} -set {path_to_init_file}"

                    multiwfn_commands = ['15', '-1'] + self.dict_of_multipole[charge_type] + ['0', 'q']
                    num_atoms = self.count_atoms(final_structure_file)

                    if charge_type == 'Hirshfeld_I':
                        # Get the number of basis functions
                        num_basis = MoldenObject(file_path_xyz, molden_filename).countBasis()
                        atmrad_src = atmrad_path
                        copy_tree(atmrad_src, os.getcwd() + '/atmrad/')
                        print(f'Current num of basis is: {num_basis}')
                        print(f'The current max num is: {self.config["maxIHirshFuzzyBasis"]}')
                        if num_basis > self.config['maxIHirshFuzzyBasis']:
                            print(f'Number of basis functions: {num_basis}')
                            multiwfn_commands = ['15', '-1'] + ['4', '-2', '1', '2'] + ['0', 'q']
                            print(f'I-Hirshfeld command should be low memory and slow to accommodate large system')

                    # Use centralized Multiwfn runner
                    self.run_multiwfn(
                        command=command,
                        input_commands=multiwfn_commands,
                        output_file=file_path_multipole,
                        description=f"Multipole {charge_type} calculation for {f}"
                    )

                else:
                    command = f"{multiwfn_path} {molden_filename}"
                    chg_prefix, _ = os.path.splitext(molden_filename)
                    calc_command = self.dict_of_calcs[charge_type]
                    commands = ['7', calc_command, '1', 'y', '0', 'q']  # for atomic charge type

                    if charge_type == 'CHELPG':
                        commands = ['7', calc_command, '1', '\n', 'y', '0', 'q']
                    elif charge_type == 'Hirshfeld_I':
                        num_basis = MoldenObject(file_path_xyz, molden_filename).countBasis()
                        print(f'Number of basis functions: {num_basis}')
                        if num_basis > self.config['maxIHirshBasis']:
                            commands = ['7', '15', '-2', '1', '\n', 'y', '0', 'q']
                        atmrad_src = atmrad_path
                        copy_tree(atmrad_src, os.getcwd() + '/atmrad/')

                    # Use centralized Multiwfn runner
                    self.run_multiwfn(
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
                return -1  # Return -1 to indicate failure

        os.chdir(owd)
        return comp_cost

    @staticmethod
    def count_atoms(xyz_filename):
        """Count number of atoms in XYZ file (mapcount replacement)."""
        with open(xyz_filename, 'r') as f:
            for i, line in enumerate(f):
                pass
        return i + 1  # Total lines including first two header lines

    @staticmethod
    def getmultipoles(multipole_name):
        """
        Parse multipole moments from Multiwfn output file.

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

                atomDict = {
                    "Index": index,
                    "Element": element,
                    "Atom_Charge": atomic_charge,
                    'Dipole_Moment': dipole_moment,
                    'Quadrupole_Moment': quadrupole_moment
                }
                atomicDicts.append(atomDict)

        return atomicDicts
