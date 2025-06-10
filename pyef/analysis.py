import re
import os
import mmap
import glob
import shutil
import logging
import traceback
import subprocess
import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
from collections import deque
from importlib import resources
from distutils.dir_util import copy_tree
import openbabel
from biopandas.pdb import PandasPdb
from .geometry import Geometry
import math
import time
class Electrostatics:
    '''
    A class to compute Electrostatic properties for series of computations housed in folders that contain .molden output files

    ...

    Attributes
    ----------
    lst_of_folders: list of strings
        list contains the name of all folders in which information for computation is contained
    lst_of_tmcm_idx: list of integers
        list contains integer indices of the metal atom at which ESP will be computed (typically center metal of TMC)
    folder_to_file_path: string
        string that indicates the filepath starting from the folder location (in lst_of_folders) to the .molden file
    hasECP: boolean
        indicates if an effective core potential was used... in this case, molden file will need to be re-formatted to be compatible multiwfn!
    includePtChgs: boolean
        indicates if point charges should be included in ESP calculation

    '''
    def __init__(self, lst_of_folders, lst_of_tmcm_idx, folder_to_file_path, hasECP=False, includePtChgs=False):
        self.lst_of_folders = lst_of_folders
        self.lst_of_tmcm_idx = lst_of_tmcm_idx
        self.folder_to_file_path = folder_to_file_path
        self.dict_of_calcs =  {'Hirshfeld': '1', 'Voronoi':'2', 'Mulliken': '5', 'Lowdin': '6', 'SCPA': '7', 'Becke': '10', 'ADCH': '11', 'CHELPG': '12', 'MK':'13', 'AIM': '14', 'Hirshfeld_I': '15', 'CM5':'16', 'EEM': '17', 'RESP': '18', 'PEOE': '19'}
        self.dict_of_multipole = {'Hirshfeld': ['3', '2'], 'Hirshfeld_I': ['4', '1', '2'],   'Becke': ['1', '2'] }
        self.dielectric_scale = 1
        self.dielectric = 1
        self.ptChgs = includePtChgs
        self.ptChgfp = ''
        self.ptchgdf = None
        self.molden_filename = 'final_optim.molden'
        self.xyzfilename = 'final_optim.xyz'
        #default setting does not generate PDB files
        self.excludeAtomfromEcalc = []
        #To avoid over-estimating screening from bound atoms, set dielectric to 1 for primary bound atoms in ESP calv
        self.changeDielectBoundBool = False
        # Dictionary is originally from molsimplify, # Data from http://www.webelements.com/ (last accessed May 13th 2015)
        # Palladium covalent radius seemed to be under-estimated in original implementation, so changed to 1.39 per https://webelements.com/palladium/atom_sizes.html
        # Dictionary from molsimplify, https://molsimplify.readthedocs.io/en/latest/_modules/molSimplify/Classes/globalvars.html
        # Some covalent radii updated with 2008 updated values from webelements.com accessed 1/18/24 
        #Note that full capitalization also included since some software write ou .xyz with all caps
        self.periodic_table = {"H": 1, "He": 2,
            "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
            "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
            "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
            "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
            "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
            "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
            "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
            "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
            "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
            "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
            "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118}
        self.amassdict = {'X': (1.0, 0, 0.77, 0), 'H': (1.0079, 1, 0.37, 1),
             'D': (2.0141, 1, 0.37, 1), 'He': (4.002602, 2, 0.46, 2),
             'Li': (6.94, 3, 1.33, 1), 'Be': (9.0121831, 4, 1.02, 2), 'B': (10.83, 5, 0.85, 3),
             'C': (12.0107, 6, 0.77, 4), 'N': (14.0067, 7, 0.75, 5), 'O': (15.9994, 8, 0.73, 6),
             'F': (18.9984, 9, 0.71, 7), 'Ne': (20.1797, 10, 0.67, 8), 'Na': (22.99, 11, 1.55, 1),
             'Mg': (24.30, 12, 1.39, 2), 'Al': (26.98, 13, 1.26, 3), 'Si': (28.08, 14, 1.16, 4),
             'P': (30.9738, 15, 1.06, 5), 'S': (32.065, 16, 1.02, 6), 'Cl': (35.453, 17, 0.99, 7),
             'Ar': (39.948, 18, 0.96, 8), 'K': (39.10, 19, 1.96, 1), 'Ca': (40.08, 20, 1.71, 2),
             'Sc': (44.96, 21, 1.7, 3), 'Ti': (47.867, 22, 1.36, 4), 'V': (50.94, 23, 1.34, 5),
             'Cr': (51.9961, 24, 1.27, 6), 'Mn': (54.938, 25, 1.39, 7), 'Fe': (55.84526, 26, 1.32, 8),
             'Co': (58.9332, 27, 1.26, 9), 'Ni': (58.4934, 28, 1.24, 10), 'Cu': (63.546, 29, 1.38, 11),
             'Zn': (65.39, 30, 1.31, 12), 'Ga': (69.72, 31, 1.26, 3), 'Ge': (72.63, 32, 1.22, 4),
             'As': (74.92, 33, 1.21, 5), 'Se': (78.96, 34, 1.16, 6), 'Br': (79.904, 35, 1.14, 7),
             'Kr': (83.798, 36, 1.17, 8), 'Rb': (85.47, 37, 2.10, 1), 'Sr': (87.62, 38, 1.85, 2),
             'Y': (88.91, 39, 1.63, 3), 'Zr': (91.22, 40, 1.48, 4), 'Nb': (92.91, 41, 1.47, 5),
             'Mo': (95.96, 42, 1.45, 6), 'Tc': (98.9, 43, 1.56, 7), 'Ru': (101.1, 44, 1.25, 8),
             'Rh': (102.9, 45, 1.35, 9), 'Pd': (106.4, 46, 1.38, 10), 'Ag': (107.9, 47, 1.45, 11),
             'Cd': (112.4, 48, 1.48, 12), 'In': (111.818, 49, 1.42, 3), 'Sn': (118.710, 50, 1.41, 4),
             'Sb': (121.760, 51, 1.40, 5), 'Te': (127.60, 52, 1.99, 6), 'I': (126.90447, 53, 1.40, 7),
             'Xe': (131.293, 54, 1.31, 8), 'Cs': (132.9055, 55, 2.44, 1), 'Ba': (137.327, 56, 1.96, 2),
             'La': (138.9, 57, 1.69, 3), 'Ce': (140.116, 58, 1.63, 4), 'Pr': (140.90766, 59, 1.76, 5),
             'Nd': (144.242, 60, 1.74, 6), 'Pm': (145, 61, 1.73, 7), 'Sm': (150.36, 62, 1.72, 8),
             'Eu': (151.964, 63, 1.68, 9), 'Gd': (157.25, 64, 1.69, 10), 'Tb': (158.92535, 65, 1.68, 11),
             'Dy': (162.500, 66, 1.67, 12), 'Ho': (164.93033, 67, 1.66, 13), 'Er': (167.259, 68, 1.65, 14),
             'Tm': (168.93422, 69, 1.64, 15), 'Yb': (173.045, 70, 1.70, 16), 'Lu': (174.9668, 71, 1.62, 3),
             'Hf': (178.5, 72, 1.50, 8), 'Ta': (180.9, 73, 1.38, 5), 'W': (183.8, 74, 1.46, 6),
             'Re': (186.2, 75, 1.59, 7), 'Os': (190.2, 76, 1.28, 8), 'Ir': (192.2, 77, 1.37, 9),
             'Pt': (195.1, 78, 1.36, 10), 'Au': (197.0, 79, 1.44, 11), 'Hg': (200.6, 80, 1.49, 2),
             'Tl': (204.38, 81, 1.44, 3), 'Pb': (207.2, 82, 1.47, 4), 'Bi': (208.9804, 83, 1.51, 5),
             'Po': (208.98, 84, 1.90, 6), 'At': (209.99, 85, 2.00, 7), 'Rn': (222.6, 86, 142, 4),
             'Fr': (223.02, 87, 3.48, 8), 'Ra': (226.03, 88, 2.01, 2), 'Ac': (277, 89, 1.86, 3),
             'Th': (232.0377, 90, 1.75, 4), 'Pa': (231.04, 91, 2.00, 5), 'U': (238.02891, 92, 1.70, 6),
             'Np': (237.05, 93, 1.90, 7), 'Pu': (244.06, 94, 1.75, 8), 'Am': (243.06, 95, 1.80, 9),
             'Cm': (247.07, 96, 1.69, 10), 'Bk': (247.07, 97, 1.68, 11), 'Cf': (251.08, 98, 1.68, 12)}        
        self.prepData()

    def includePtChgs(self, name_ptch_file):
        ''' Function to include point charges in ESP calculation
        Input: name_ptch_file: string of point charge filename
        '''
        self.ptChgfp = name_ptch_file
        self.ptChgs = True
        print(f'Point charges to be included via {name_ptch_file}')
    def set_dielec_scale(self, dielec):
        self.dielectric_scale = dielec
    def excludeAtomsFromEfieldCalc(self, atom_to_exclude):
        ''' Function to exclude atoms from Efield calculation
        Input: atom_to_exclude: list of integers of atom indices
        '''
        self.excludeAtomfromEcalc = atom_to_exclude

    def minDielecBonds(self, bool_bonds):
        ''' Function to change dielectric of bound atoms to 1
        Input: bool_bonds: boolean to change dielectric of bound atoms to 1
        '''
        self.changeDielectBoundBool = bool_bonds

    def runlowmemory(self):
        ''' Run the lower memory, longer time simulations. Useful when using Hirshfeld-I for simulations with more than 300 atoms
        '''
        self.dict_of_multipole = {'Hirshfeld': ['3', '2'], 'Hirshfeld_I': ['4', '-2', '1', '2'],   'Becke': ['1', '2'] }
        self.dict_of_calcs = {'Hirshfeld': '1', 'Voronoi':'2', 'Mulliken': '5', 'Lowdin': '6', 'SCPA': '7', 'Becke': '10', 'ADCH': '11', 'CHELPG': '12', 'MK':'13', 'AIM': '14', 'Hirshfeld_I': ['15','-2'], 'CM5':'16', 'EEM': '17', 'RESP': '18', 'PEOE': '19'}

    def changeDielectric(self, dlc):
        ''' Function to change dielectric of solvent
        Input: dlc: float   dielectric constant of solvent
        '''
        self.dielectric = dlc

    def set_molden_filename(self, new_name):
        self.molden_filename = new_name
        self.prepData()
    def set_xyzfilename(self, new_name):
        self.xyzfilename = new_name
        self.prepData()

    def fix_ECPmolden(self):
        """Prepares output terachem data for analysis, mainly isolating final .xyz frame and naming .molden file appropriotely"""

        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        owd = os.getcwd()
        print('   > Re-formatting .molden file to fix ECP artifacts')
        for f in list_of_folders:
            os.chdir(owd)
            print('   > Changing directory: ' + str(f + folder_to_molden))
            os.chdir(f + folder_to_molden)
            with open(self.molden_filename, 'r') as file:
                content = file.read()
            pattern_au = re.compile(r'(Au\s+\d+\s+)(\d+)')
            content = pattern_au.sub(r'\g<1>19', content)

            pattern_fe = re.compile(r'(Fe\s+\d+\s+)(\d+)')
            content = pattern_fe.sub(r'\g<1>8', content)

            pattern_i = re.compile(r'(I\s+\d+\s+)(\d+)')
            content = pattern_i.sub(r'\g<1>7', content)

            with open(self.molden_filename, 'w') as file:
                file.write(content)
            print("      > Molden file is fixed\n")
        os.chdir(owd)

    def prepData(self):
        """Prepares output terachem data for analysis, mainly isolating final .xyz frame and naming .molden file appropriotely"""

        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        owd = os.getcwd()
        print('   > Pre-processing data')
        backup_xyz = 'xyz.xyz'

        for f in list_of_folders:
            folder_path = os.path.join(owd, f + folder_to_molden)
            print('      > .molden and .xyz file should be located here: ' + folder_path)

            # Processing optim.xyz to create final_optim.xyz
            final_optim_xyz = os.path.join(folder_path, self.xyzfilename)

            # Copying .molden files to final_optim.molden
            final_optim_molden = os.path.join(folder_path, self.molden_filename)

            #this is for the full optimization cycle
            optim_file_path = os.path.join(folder_path, 'optim.xyz')


            if not os.path.exists(final_optim_molden):
                print(f'Expected .molden with filename: {self.molden_filename} in {folder_path}. We could not find a .molden filename with the default prefix, you can alter using: set_molden_filename()')
                print(f"For now searching for all .molden files in directory. If you only have one .molden file in this directory, the defauly prefix will be altered ")
                try:
                    files = glob.iglob(os.path.join(folder_path, "*.molden"))
                    for file in files:
                        if os.path.abspath(file) != os.path.abspath(final_optim_molden):
                            #change the defauly path the the molden file!
                            self.set_molden_filename(os.path.basename(file))
                            #set the backup_xyz filename to this prefix here. Generally a good guess!
                            file_prefix, _ = os.path.splitext(os.path.basename(file))
                            backup_xyz = file_prefix + '.xyz'
                            print(f'Default .molden file name is now changed to {self.molden_filename}')
                            break
                except Exception as e:
                    logging.exception('An Exception was thrown while copying molden files.')
            else:
                print(f"      > {self.molden_filename} sucessflly located in {folder_path}.")



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
                    print(f'Expected .xyz with filename: {self.xyzfilename} in {folder_path}. If your xyz filename does NOT match the default, you can alter using: set_xyzfilename()')
                    print(f'We will try to use the anticipated prefix from the associated molden file: {backup_xyz}')
                    backup_file = os.path.join(folder_path, backup_xyz)
                    if os.path.exists(backup_file):
                        self.xyzfilename = backup_xyz
                        logging.info(f'Single point data found. Using {backup_xyz} as fallback.')

                    else:
                        logging.exception(f'An unexpected error occurred while processing optim.xyz: {e}')
            else:
                print(f'      > {self.xyzfilename} succesfully located in {folder_path}.')            

            # Copying .molden files to final_optim.molden
            final_optim_molden = os.path.join(folder_path, self.molden_filename)
            backup_xyz = 'final_optim'

        os.chdir(owd)


    def getmultipoles(multipole_name):
        '''
        Input: multipole_name: string of multipole filename
        Output: list of dictionaries of multipole moments
        '''
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
    #Function will process a point charge file (as generated by amber) and return a dataframe with partial charges and coordinates for use in ESP calculation
    def getPtChgs(self, filename_pt):
        '''
        Input: filename_pt: string of point charge filename
        Output: dataframe with partial charges and coordinates
        '''
        chg_df = pd.read_table(filename_pt, skiprows=2, delim_whitespace=True, names=['charge', 'x', 'y', 'z'])
        atm_name = ['pnt']
        atoms = atm_name*len(chg_df['charge'])
        chg_df['Atom'] = atoms
        self.ptchgdf  = chg_df
        return chg_df

    # Define the functions to calculate the ESP:
    def mapcount(filename):
        """Function to rapidly count the number of lines in a file
        Input: filename: string of filename
        Output: integer of number of lines in file
        """

        f = open(filename, "r+")
        buf = mmap.mmap(f.fileno(), 0)
        lines = 0
        readline = buf.readline

        while readline():
            lines += 1
        return lines


    def calcesp(self, path_to_xyz, espatom_idx, charge_range, charge_file):
        """
        Calculate the esp in units of Volts
        Input: path_to_xyz: string of xyz filename
        espatom_idx: integer of atom index
        charge_range: list of integers of atom indices
        charge_file: string of charge filename
        Output: list of ESP and atomic symbol

        Notes
        -----
        Espatom_idx should be the index of the atom at which the esp should be calculated
        Charge_range should provide an array with the index of each atom to be considered in calculation of the ESP
        Charge_file is the filepath for the .txt file containing the hirshfeld charges (generated by multiwfn)

        """
        dielectric = self.dielectric
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2) aka N*m^2/C^2

        # Convert each column to list for quicker indexing
        atoms = list(df['Atom'])
        charges = list(df['charge'])
        xs = list(df['x'])
        ys = list(df['y'])
        zs = list(df['z'])

        #For QMMM calculation, include point charges in ESP calculation
        if self.ptChgs:
            df_ptchg  = self.ptchgdf
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

        # Unit conversion
        A_to_m = 10**(-10)
        #faraday = 23.06   #kcal/(mol*V)
        C_e = 1.6023*(10**-19)
        #cal_J = 4.184

        bound_atoms = []
        #create list of bound atoms, these are treated with a different dielectric
        if self.changeDielectBoundBool:
            bound_atoms = self.getBondedAtoms(path_to_xyz, idx_atom)
 
        for idx in charge_range:
            if idx == idx_atom:
                continue
            elif idx in bound_atoms:
                #now account for bound atoms
                 r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                 total_esp = total_esp + (charges[idx]/r)
            else:
                # Calculate esp and convert to units (A to m)
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                total_esp = total_esp + (1/dielectric)*(charges[idx]/r)

        final_esp = k*total_esp*((C_e))  #N*m^2/(C^2)*(C/m) = N*m/C = J/C = Volt
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
        k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

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

        # Unit conversion
        A_to_m = 10**(-10)
        C_e = 1.6023*(10**-19)
        bohr_to_m = 0.52918*10**(-10) 
        for idx in charge_range:
            if idx == idx_atom:
                continue
            else:
                # Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                Ex = Ex - k*C_e*(charges[idx]/r)*(1/(xs[idx] - xo))   
                Ey = Ey - k*C_e*(charges[idx]/r)*(1/(ys[idx] - yo))
                Ez = Ez - k*C_e*(charges[idx]/r)*(1/(zs[idx] - zo))

        E_vec = [Ex, Ey, Ez]
        return [E_vec, position_vec, df['Atom'][idx_atom]]

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

    def update_mm_charges_based_on_esp(esp, mm_charges, alpha=0.01):
        '''
        Update the MM charges based on the ESP
        '''
        delta_q = alpha * esp  # Simple linear response
        updated_mm_charges = mm_charges + delta_q
        return updated_mm_charges

    def update_mm_charges_drude(esp, mm_charges, polarization, beta=0.01):
        delta_p = beta * esp  # Induced dipole change
        polarization += delta_p
        updated_mm_charges = mm_charges + polarization
        return updated_mm_charges, polarization


    def resp_correction_objective(delta_q, mm_coords, q_mm_orig, qm_coords, esp_target, restraint=0.1):
        """
        Objective function: least-squares ESP error + restraint term.
        """
        q_new = q_mm_orig + delta_q
        print(f'delta q: {delta_q}')
        esp_fit = Electrostatics.compute_esp(mm_coords, q_new, qm_coords)
        esp_err = np.sum((esp_fit - esp_target) ** 2)
        penalty = restraint * np.sum(q_mm_orig ** 2)
        restraint_term = restraint * np.sum(delta_q ** 2)
        return esp_err + restraint_term

    def correct_mm_charges( qm_coords, qm_charges, mm_coords, mm_charges):
        """
        Apply minimal correction to existing MM charges to match QM ESP.
        """

        #esp_target = Electrostatics.compute_esp(qm_coords, qm_charges, mm_coords)
        #delta_q0 = np.zeros_like(q_mm_orig)

        esp_from_QM = Electrostatics.compute_esp_from_qm(qm_coords, qm_charges, mm_coords)
        alpha = 1/2.97
        new_chgs = Electrostatics.update_mm_charges_based_on_esp(esp_from_QM, mm_charges, alpha)

        return list(new_chgs)



    def calc_fullE(self, idx_atom, charge_range, xyz_file, atom_multipole_file):
        '''
        
        Input: idx_atom: integer of atom index
        charge_range: list of integers of atom indices
        xyz_file: string of xyz filename
        atom_multipole_file: string of multipole filename
        Output: list of E-field vector and atomic symbol, Efields are reported in units of Volts/Angstrom
        '''

        df = self.getGeomInfo(xyz_file)
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2 =  N*m^2/C^2)

        # Convert each column to list for quicker indexing
        xs = df['X']
        ys = df['Y']
        zs = df['Z']

        
        # Following derivation of E-field strength 
        Monopole_E = np.array([0, 0, 0])

        # Unit conversion
        A_to_m = 10**(-10)
        # Bohr to meters (atomic units)
        b_to_m = 5.291772109*(10**-11)
        b_to_A = 0.529177
        C_e = 1.6023*(10**-19)
        Vm_to_VA = 10**(-10)
        Ex = 0
        Ey = 0
        Ez = 0

        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]

        position_vec = A_to_m*np.array([xo, yo, zo])


        inv_eps = 1/self.dielectric

        if self.ptChgs:
            df_ptchg  = self.ptchgdf

        #load multipole moments from processed outputs 
        lst_multipole_dict = Electrostatics.getmultipoles(atom_multipole_file)

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
                dipole_vec = b_to_m*np.array(atom_dict["Dipole_Moment"])
                #convention of vector pointing towards atom of interest, so positive charges exert positive Efield
                dist_vec = A_to_m*np.array([(xs[idx] - xo), (ys[idx] - yo), (zs[idx] - zo)])
                dist_arr = np.outer(dist_vec, dist_vec)
                quadrupole_arr = b_to_m*b_to_m*atom_dict['Quadrupole_Moment']
                
                # Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                Monopole_E = Monopole_E - inv_eps*k*C_e*atom_dict["Atom_Charge"]*(1/(r**3))*dist_vec  #neg for sign convention so Efield points toward neg charge
                Ex_quad = inv_eps*k*C_e*(1/(r**3))*(A_to_m*(xs[idx] - xo))*(1/r**4)*dist_arr[1:, 1:]*quadrupole_arr[1:,1:]
                Ey_quad = inv_eps*k*C_e*(1/(r**3))*(A_to_m*(ys[idx] - yo))*(1/r**4)*dist_arr[0:2:3, 0:2:3]*quadrupole_arr[0:2:3,0:2:3]
                Ez_quad = inv_eps*k*C_e*(1/(r**3))*(A_to_m*(zs[idx] - zo))*(1/r**4)*dist_arr[0:2, 0:2]*quadrupole_arr[0:2,0:2]

                Ex = Ex + inv_eps*k*(1/(r**3))*dist_vec[0]*( -C_e*atom_dict["Atom_Charge"] + C_e*(1/r**2)*np.dot(dipole_vec[1:], dist_vec[1:]))-(1/3)*Ex_quad.sum()
                Ey = Ey + inv_eps*k*(1/(r**3))*dist_vec[1]*( -C_e*atom_dict["Atom_Charge"] + C_e*(1/r**2)*np.dot(dipole_vec[0:2:3], dist_vec[0:2:3])) -(1/3)*Ey_quad.sum()
                Ez = Ez + inv_eps*k*(1/(r**3))*dist_vec[2]*( -C_e*atom_dict["Atom_Charge"] + C_e*(1/r**2)*np.dot(dipole_vec[0:2], dist_vec[0:2])) -(1/3)*Ez_quad.sum()
        E_vec = [Ex, Ey, Ez]


        #For QMMM calculation, include point charges in E field calc
        if self.ptChgs:
            MM_xs = list(df_ptchg['x'])
            MM_ys = list(df_ptchg['y'])
            MM_zs = list(df_ptchg['z'])
            init_MM_charges = list(df_ptchg['charge'])
 
            #make new MM_charges by scaling!
            mm_coords = np.column_stack((np.array(MM_xs), np.array(MM_ys), np.array(MM_zs)))
            mm_charges = np.array(init_MM_charges)
            qm_charges = np.array(QM_charges)
            print(f'Size of mm_coords: {np.shape(mm_coords)} and qm coords: {np.shape(QM_coords)}; mm charges: {np.shape(mm_charges)} and qm_charges: {np.shape(qm_charges)}')
            MM_charges = Electrostatics.correct_mm_charges(QM_coords, qm_charges, mm_coords, mm_charges)

            MM_charges = mm_charges*(1/(math.sqrt(self.dielectric_scale)))
            #change charge range to include these new partial charges!
            charge_range = range(0, len(MM_xs))
            for chg_idx in charge_range:
                r = (((MM_xs[chg_idx] - xo)*A_to_m)**2 + ((MM_ys[chg_idx] - yo)*A_to_m)**2 + ((MM_zs[chg_idx] - zo)*A_to_m)**2)**(0.5)
                dist_vec = A_to_m*np.array([(MM_xs[chg_idx] - xo), (MM_ys[chg_idx] - yo), (MM_zs[chg_idx] - zo)])
                E_to_add = -inv_eps*k*C_e*MM_charges[chg_idx]*(1/(r**3))*dist_vec
                Monopole_E = Monopole_E +  E_to_add
                E_vec[0] += E_to_add[0]
                E_vec[1] += E_to_add[1]
                E_vec[2] += E_to_add[2]
            #Add contributions to Monopole E from point charges to total E
        return [Vm_to_VA*np.array(E_vec), position_vec, df['Atom'][idx_atom], Vm_to_VA*np.array(Monopole_E) ]
    
    def ESPfromMultipole(self, xyfilepath, atom_multipole_file, charge_range, idx_atom):
        '''

        Input: idx_atom: integer of atom index
        charge_range: list of integers of atom indices
        xyfilepath: string of xyz filename
        atom_multipole_file: string of multipole filename
        Output: list of ESP and atomic symbol Units of ESP is Volts!
        '''
        df = self.getGeomInfo(xyfilepath)
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2) also N*m^2/C^2

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
        A_to_m = 10**(-10)
        #KJ_J = 10**-3
        #faraday = 23.06   #kcal/(mol*V)
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        #cal_J = 4.184
        dielectric = self.dielectric
        lst_multipole_dict = Electrostatics.getmultipoles(atom_multipole_file)
        #create list of bound atoms, these are treated with a different dielectric
        bound_atoms = self.getBondedAtoms(xyfilepath, idx_atom)
        total_esp= 0
        for idx in charge_range:
            atom_dict = lst_multipole_dict[idx]
            if idx == idx_atom:
                continue
            elif idx in bound_atoms:
                #default it to exclude bound atoms
                 r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                 total_esp = total_esp + (atom_dict["Atom_Charge"]/r)
            else:
                # Calculate esp and convert to units (A to m)
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                total_esp = total_esp + (1/dielectric)*(atom_dict["Atom_Charge"]/r)

        final_esp = k*total_esp*((C_e))  #Units: N*m^2/(C^2)*(C/m) = N*m/C = J/C = Volt
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
        A_to_m = 10**(-10)
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
            
    def esp_first_coord(self, metal_idx, charge_file, path_to_xyz):
        ''' Calculate ESP accounting for electrostatic contributions for atoms bound to ESP center
        Input: metal_idx: integer of atom index
        charge_file: string of charge filename
        path_to_xyz: string of xyz filename
        Output: list of ESP and atomic symbol
        '''

        print('The index of the metal atom is: ' + str(metal_idx))
        lst_bonded_atoms = self.getBondedAtoms(path_to_xyz, metal_idx)
        [First_coord_ESP, atom_type] = self.calcesp(path_to_xyz, metal_idx, lst_bonded_atoms, charge_file)
        return First_coord_ESP
    
    def esp_second_coord(self, metal_idx, charge_file, path_to_xyz):
        '''A Function to calculate ESP including contributions only from first and second coordination spheres
        Input: metal_idx: integer of atom index
            charge_file: string of charge filename
            path_to_xyz: string of xyz filename
        Output: list of ESP and atomic symbol
            '''
        lst_first_and_second = []
        lst_first_coor_atoms = self.getBondedAtoms(path_to_xyz, metal_idx)
        lst_first_and_second.extend(lst_first_coor_atoms)
        for coor_atom_idx in lst_first_coor_atoms:
            second_coor = self.getBondedAtoms(path_to_xyz, coor_atom_idx)
            lst_first_and_second.extend(second_coor)
        set_second_coor = set(lst_first_and_second)
        final_lst = list(set_second_coor)
        [second_coord_ESP, atom_type] = self.calcesp(path_to_xyz, metal_idx, final_lst, charge_file)
        return second_coord_ESP

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
            dict_multipole = Electrostatics.getmultipoles(chg_filename)  #Index, Element, Atom_Charge, Dipole_Moment....
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
        dielectric = self.dielectric
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2) = N*m^2/(C^2)

        # Unit conversion
        A_to_m = 10**(-10)
        KJ_J = 10**-3
        faraday = 23.06   #kcal/(mol*V)
        C_e = 1.6023*(10**-19)
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
        if self.ptChgs:
            ptchg_filename = self.ptChgfp
            init_file_path = path_to_xyz[0:-len('scr/final_optim.xyz')]
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
        if self.changeDielectBoundBool:
            bound_atoms = self.getBondedAtoms(path_to_xyz, idx_atom)
        for idx in range(0, total_atoms):
            if idx == idx_atom:
                continue
            elif idx in bound_atoms:
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                esps.append(k*C_e*charges[idx]/r) #units of N*m/C =Volt
                distances.append(r)
            else:
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                distances.append(r)
                esps.append(k*(1/dielectric)*C_e*charges[idx]/r) #untis of N*m/C = Volt
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

    def getESPDecay(self, charge_types, ESPdata_filename, multiwfn_module, multiwfn_path, atmrad_path, dielectric=1):
        '''
        Function computes a series of ESP data using the charge scheme specified in charge types.
        Inputs:
        ------- 
        charge_types: list of strings, possible choices include: 'Hirshfeld', 'Voronoi', 'Mulliken',  'Lowdin', 'SCPA', 'Becke', 'ADCH', 'CHELPG', 'MK', 'AIM', 'Hirshfeld_I', 'CM5', 'EEM', 'RESP', 'PEOE'}
        ESPdata_filename: string
            Name of the output file name
        multiwfn_module: string
            Name of the module that contains the multiwfn executable
        multiwfn_path: string
            Path to the multiwfn executable
        atmrad_path: string
            Path to the atmrad executable
        dielectric: float
            Dielectric constant of the solvent
    
        Outputs:
        -------
        ESPdata_filename: string
            Name of the output file name    '''

        self.dielectric = dielectric
       # Access Class Variables
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # Old working directory
        allspeciesdict = []
        counter = 0  # Iterator to account for atomic indices of interest
        for f in list_of_file:
            print('-----------------' + str(f) + '------------------')
            atom_idx = metal_idxs[counter]
            counter = counter + 1
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            subprocess.call(multiwfn_module, shell=True)
            command_A = f"{multiwfn_path} final_optim.molden"
            results_dir = os.getcwd() + '/'
            
            results_dict = {}
            results_dict['Name'] = f

            for key in charge_types:
                print('Partial Charge Scheme:' + str(key))
                try:
                    full_file_path = f"{os.getcwd()}/final_optim_{key}.txt"
                    path_to_xyz = f"{os.getcwd()}/final_optim.xyz"
                    if key == "Hirshfeld_I":
                        atmrad_src = atmrad_path
                        copy_tree(atmrad_src, results_dir + 'atmrad/')
                    try: 
                        [ESP_all, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx)

                        [total_charge,partial_charge_atom] = Electrostatics.charge_atom(full_file_path, atom_idx)
                        [sorted_distances, sorted_esps, cum_esps, sorted_cum_idx, sorted_cum_chg, sorted_atomTypes] = self.esp_bydistance(path_to_xyz, atom_idx, full_file_path)
                        ESP_fcoord = self.esp_first_coord(atom_idx, full_file_path, path_to_xyz)
                        ESP_scoord = self.esp_second_coord(atom_idx, full_file_path, path_to_xyz)
                    
                    except Exception as e:
                        print('The Exception is: ' + str(e))
                        print(traceback.format_exc())
                        print('Error when trying to access electrostatic information: Attemtping to re-compute partial charges of type: ' + str(key))

                        # Re-run multiwfn computation of partial charge 
                        proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                        calc_command = self.dict_of_calcs[key]
                        commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                        if key == 'CHELPG':
                            commands = ['7', calc_command, '1','\n', 'y', '0', 'q']
                        output = proc.communicate("\n".join(commands).encode())
                        new_name = 'final_optim_' +key+'.txt'
                        os.rename('final_optim.chg', new_name)
             
                        [ESP_all, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx, self.inGaCageBool)
                        
                        [total_charge,partial_charge_atom] = Electrostatics.charge_atom(full_file_path, atom_idx)
                        [sorted_distances, sorted_esps, cum_esps, sorted_cum_idx, sorted_cum_chg, sorted_atomTypes] = self.esp_bydistance(path_to_xyz, atom_idx, full_file_path)
                        ESP_fcoord = self.esp_first_coord(atom_idx, full_file_path, path_to_xyz)
                        ESP_scoord = self.esp_second_coord(atom_idx, full_file_path, path_to_xyz)

                    # At this point, all calculations shouldbe complete and succesfull: Add ESP data to dictionary
                    results_dict[str(key) + ' ESP Second Coor Shell (kcal/mol)'] = ESP_scoord
                    results_dict[str(key) + ' ESP First Coor Shell (kcal/mol)'] = ESP_fcoord
                    results_dict['Atoms'] = atom_type
                    results_dict['Total Charge'] = total_charge
                    results_dict['Partial Charge '+str(key)] = partial_charge_atom
                    results_dict['ESP '+ str(key)] = ESP_all
                    results_dict['Sorted Distances'] = sorted_distances
                    results_dict['Sorted ESP '+ str(key)] = sorted_esps
                    results_dict['Cumulative ESP ' + str(key)] = cum_esps
                    results_dict['Dist Sorted Idxs' + str(key)] = sorted_cum_idx
                    results_dict['Dist Sorted Partial Charges' + str(key)] = sorted_cum_chg
                    results_dict['Dist Sorted Atom Types' + str(key)] = sorted_atomTypes

                except Exception as e:
                    logging.exception('An Exception was thrown')
                    continue
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(ESPdata_filename +'.csv')
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
            print(f'Num solvents: {num_solvent_res}')
            for solv_num in range(0, num_solvent_res):
                res_dict[f'Solvent_{solv_num}'] = solvent_idxs[solv_num*num_atoms_solvent:(solv_num+1)*(num_atoms_solvent)]
        return res_dict



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
        print(f'Currently the owd is: {owd}')
        molden_filename = self.molden_filename
        final_structure_file = self.xyzfilename

        for f in list_of_folders:
            comp_cost = 'Na'
            num_atoms = 0
            print(f'Current folder: {f}')
            need_to_run_calculation = True
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            #subprocess.call(multiwfn_module, shell=True)
            file_path_multipole = f"{os.getcwd()}/Multipole{charge_type}.txt"
            file_path_monopole = f"{os.getcwd()}/Charges{charge_type}.txt"
            file_path_xyz = f"{os.getcwd()}/{final_structure_file}"

            #check if previous calculations fully converged for desired multipole/charge scheme 
            if multipole_bool:
                if os.path.exists(file_path_multipole):
                    with open(file_path_multipole, 'r') as file:
                        contents = file.read()
                        if "Calculation took up" in contents:
                            print(f" > Previously completed multipole {charge_type} calculation for {f}")
                            need_to_run_calculation = False
            else:
                if os.path.exists(file_path_monopole):
                    need_to_run_calculation = False
                    #Dynamically get path to package settings.ini file
                    #path_to_ini_file = str(ini_path)
            

            #If you need to run the calculations, urn either the multipole or the monopole calculation!
            if need_to_run_calculation:
                print(f'Running Calculation')
                try: 
                    start = time.time()
                    if multipole_bool:
                        chg_prefix, _ = os.path.splitext(molden_filename)
                        #Dynamically get path to package settings.ini file
                        with resources.path('pyef.resources', 'settings.ini') as ini_path:
                            path_to_init_file = str(ini_path)
                            Command_Multipole = f"{multiwfn_path} {molden_filename} -set {path_to_init_file} > {file_path_multipole}"
                        proc = subprocess.Popen(Command_Multipole, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                        multiwfn_commands = ['15', '-1'] + self.dict_of_multipole[charge_type] + ['0', 'q']
                        if charge_type == 'Hirshfeld_I' and  num_atoms > 320:
                            multiwfn_commands = ['15', '-1'] + ['4', '-2', '1', '2'] + ['0', 'q'] 
                        proc.communicate("\n".join(multiwfn_commands).encode())

                    else:
                        command_A = f"{multiwfn_path} {molden_filename}"
                        chg_prefix,  _ = os.path.splitext(molden_filename)
                        proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                        calc_command = self.dict_of_calcs[charge_type]
                        commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                        if charge_type == 'CHELPG':
                            commands = ['7', calc_command, '1','\n', 'y', '0', 'q']
                        elif charge_type == 'Hirshfeld_I':
                            num_atoms = Electrostatics.mapcount(final_structure_file) - 2
                            print(f'Number of atoms: {num_atoms}')
                            if num_atoms > 720:
                                commands = ['7', '15', '-2', '1', '\n', 'y', '0', 'q']
                            atmrad_src = atmrad_path
                            copy_tree(atmrad_src, os.getcwd() + '/atmrad/')
                            #if too many atoms will need to change calc to run with reasonable memory 
                        output = proc.communicate("\n".join(commands).encode())
                        os.rename(f'{chg_prefix}.chg', file_path_monopole)
                    end = time.time()
                    comp_cost = end - start
                except Exception as e:
                    print(e)
                    #Issue could be from lost memorry
                    os.chdir(owd)
                    continue
                    
            if multipole_bool:
                xyz_fp =  final_structure_file
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

            sorted_dips = df_dip['Dipole']
            avg_dips = np.average(sorted_dips[1:])
            
            avg_first_5A = np.average(df_dip.loc[(df_dip['Distance from solute'] < 5), 'Dipole'][1:])
            avg_5to7A = np.average(df_dip.loc[(df_dip['Distance from solute'] < 7) & (df_dip['Distance from solute'] > 5), 'Dipole'])
            avg_7to11A = np.average(df_dip.loc[(df_dip['Distance from solute'] < 11) & (df_dip['Distance from solute'] > 7), 'Dipole'])
            avg_above11A = np.average(df_dip.loc[df_dip['Distance from solute'] > 11, 'Dipole'])


            #return a list of the solvents
             #return an average of the dipoles in first x closests; next x closest, etc.

            dip_dict = {'DipoleSolute': sorted_dips[0], 'AvgDipSolv': avg_dips, 'Avg5Ang': avg_first_5A, 'Avg5to7Ang': avg_5to7A , 'Avg7to11Ang': avg_7to11A, 'Avgabove11Ang': avg_above11A, 'CompCost': comp_cost, 'Total_atoms': total_atoms}
            print(dip_dict)
            lst_dicts.append(dip_dict)
            os.chdir(owd)
        all_file_df = pd.DataFrame(lst_dicts)
        all_file_df.to_csv('Dipoles.csv')
        return all_file_df


    def getESPData(self, charge_types, ESPdata_filename, multiwfn_module, multiwfn_path, atmrad_path, dielectric=1):
        '''
        Function computes a series of ESP data using the charge scheme specified in charge types. All ESPs in Units of Volts
        Inputs:
        ------- 
        charge_types: list of strings, possible choices include: 'Hirshfeld', 'Voronoi', 'Mulliken',  'Lowdin', 'SCPA', 'Becke', 'ADCH', 'CHELPG', 'MK', 'AIM', 'Hirshfeld_I', 'CM5', 'EEM', 'RESP', 'PEOE'}
        ESPdata_filename: string
            Name of the output file name
        multiwfn_module: string
            Name of the module that contains the multiwfn executable
        multiwfn_path: string
            Path to the multiwfn executable
        atmrad_path: string
            Path to the atmrad executable
        dielectric: float
            Dielectric constant of the solvent
    
        Outputs:
        -------
        ESPdata_filename: string
            Name of the output file name    '''

        self.dielectric = dielectric
       # Access Class Variables
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # Old working directory
        allspeciesdict = []
        counter = 0  # Iterator to account for atomic indices of interest
        for f in list_of_file:
            start_time = time.time()
            print('-----------------' + str(f) + '------------------')
            atom_idx = metal_idxs[counter]
            counter = counter + 1
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            subprocess.call(multiwfn_module, shell=True)
            command_A = f"{multiwfn_path} final_optim.molden"
            results_dir = os.getcwd() + '/'
            
            results_dict = {}
            results_dict['Name'] = f

            for key in charge_types:
                print('Partial Charge Scheme:' + str(key))
                try:
                    full_file_path = f"{os.getcwd()}/final_optim_{key}.txt"
                    path_to_xyz = f"{os.getcwd()}/final_optim.xyz"
                    if key == "Hirshfeld_I":
                        atmrad_src = atmrad_path
                        copy_tree(atmrad_src, results_dir + 'atmrad/')
                    try: 
                        [ESP_all, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx)
                        [total_charge,partial_charge_atom] = Electrostatics.charge_atom(full_file_path, atom_idx)
                    
                    except Exception as e:
                        print('The Exception is: ' + str(e))
                        print(traceback.format_exc())
                        print('Error when trying to access electrostatic information: Attempting to re-compute partial charges of type: ' + str(key))

                        # Re-run multiwfn computation of partial charge 
                        proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                        calc_command = self.dict_of_calcs[key]
                        commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                        if key == 'CHELPG':
                            commands = ['7', calc_command, '1','\n', 'y', '0', 'q']
                        output = proc.communicate("\n".join(commands).encode())
                        new_name = 'final_optim_' +key+'.txt'
                        os.rename('final_optim.chg', new_name)
             
                        [ESP_all, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx)
                        [total_charge,partial_charge_atom] = Electrostatics.charge_atom(full_file_path, atom_idx)


                    results_dict['Atoms'] = atom_type
                    results_dict['Total Charge'] = total_charge
                    results_dict['Partial Charge '+str(key)] = partial_charge_atom
                    results_dict['ESP '+ str(key)] = ESP_all


                except Exception as e:
                    logging.exception('An Exception was thrown')
                    continue
            allspeciesdict.append(results_dict)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"Elapsed time: {elapsed_time:.4f} seconds")
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(ESPdata_filename +'.csv')
        return df
    
    def getESPMultipole(self, ESP_filename, multiwfn_module, multiwfn_path, polarization_scheme='Hirshfeld_I'):
        '''
        Function computes a series of ESP data using the charge scheme specified in charge types. All ESPs are in units of Volts
        Inputs:
        -------
        ESP_filename: string
            Name of the output file name
        multiwfn_module: string
            Name of the module that contains the multiwfn executable
        multiwfn_path: string
            Path to the multiwfn executable
        polarization_scheme: string
            The scheme to use for polarization, default is 'Hirshfeld_I'. Other options can be 'Hirshfeld', 'Becke'

        Outputs:
        -------
        ESP_filename: string
            Name of the output file name

            '''

        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # Old working directory
        allspeciesdict = []
        counter = 0
        for f in list_of_file:
            start_time = time.time()
            atom_idx = metal_idxs[counter] 
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            subprocess.call(multiwfn_module, shell=True)

            # First For this to work, the .molden file should be named: f.molden
            results_dict = {}
            results_dict['Name'] = f
            multiwfn_path = multiwfn_path
            molden_filename = "final_optim.molden"
            final_structure_file = "final_optim.xyz"
            polarization_file = "final_optim_multipole" + polarization_scheme + ".txt"
            
            # Dynamically get path to package settings.ini file
            with resources.path('pyef.resources', 'settings.ini') as ini_path:
                path_to_ini_file = str(ini_path)
                Command_Polarization = f"{multiwfn_path} {molden_filename} -set {path_to_ini_file} > {polarization_file}"

            # Check if the atomic polarizations have been computed
            path_to_pol = os.path.join(os.getcwd(), polarization_file)
            backup_path_to_pol = os.path.join(os.getcwd(), "final_optim_polarization.txt")
            xyz_file_path = os.path.join(os.getcwd(), final_structure_file)
            print(f"Attempting polarization file path: {path_to_pol}")
            
            #Pick lines to include, can exclude atom indices from calculation by calling function excludeAtomsFromEfieldCalc
            total_lines = Electrostatics.mapcount(xyz_file_path)
            init_all_lines = range(0, total_lines - 2)
            all_lines = [x for x in init_all_lines if x not in self.excludeAtomfromEcalc]

            # Check the contents of the polarization file to see if it finished
            need_to_run_calculation = True
            if os.path.exists(path_to_pol):
                with open(path_to_pol, 'r') as file:
                    contents = file.read()
                    if "Calculation took up" in contents:
                        print(f"   > Polarization file: {f}!!")
                        need_to_run_calculation = False

            if need_to_run_calculation:
                print('Starting to run polarization calculation!')
                # Now Run the calculation for atomic dipole and quadrupole moment
                print(f"   > Submitting Multiwfn job using: {Command_Polarization}")
                proc = subprocess.Popen(Command_Polarization, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                polarization_commands = ['15', '-1'] + self.dict_of_multipole[polarization_scheme] + ['0', 'q']
                proc.communicate("\n".join(polarization_commands).encode())
            try:
                [final_esp, atom_name] = self.ESPfromMultipole(xyz_file_path, path_to_pol, all_lines, atom_idx)
                results_dict['Total ESP'] = final_esp
                results_dict['Atom'] = atom_name 
            except Exception as e:
                print(repr(e))
                traceback_str = ''.join(traceback.format_tb(e.__traceback__))
                print(traceback_str)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f'Time for Efield calc: {elapsed_time}')
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(f"{ESP_filename}.csv")



    # input_bond_indices is a list of a list of tuples
    def getEFieldMultipole(self, Efield_data_filename, multiwfn_module, multiwfn_path, atmrad_path, input_bond_indices=[], excludeAtoms=[], polarization_scheme='Hirshfeld_I'):
        '''
        Function computes a series of Efield data using the charge scheme specified in charge types; All in Units of Volts/Angstrom

        Inputs:
        -------
        Efield_data_filename: string
            Name of the output file name
        multiwfn_module: string
            Name of the module that contains the multiwfn executable
        multiwfn_path: string
            Path to the multiwfn executable   
        polarization_scheme: string
            The scheme to use for polarization, default is 'Hirshfeld_I'. Other options can be 'Hirshfeld', 'MBIS', 'Becke'
            
        Outputs:    
        -------        
        Efield_data_filename: string
             Name of the output file name
        '''
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # Old working directory
        allspeciesdict = []
        counter = 0
        
        # If indices of bonds of interest are entered then switch to manual mode
        bool_manual_mode = False
        if len(input_bond_indices) > 0:
            bool_manual_mode = True
        
        for f in list_of_file:
            start_time = time.time()
            print(f'--------------file {f}---------------------------')
            atom_idx = metal_idxs[counter] 
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            subprocess.call(multiwfn_module, shell=True)

            # First For this to work, the .molden file should be named: f.molden
            results_dict = {}
            results_dict['Name'] = f
            multiwfn_path = multiwfn_path
            molden_filename = "final_optim.molden"
            final_structure_file = "final_optim.xyz"
            polarization_file = "Multipole" + polarization_scheme + ".txt"
            
            # Dynamically get path to package settings.ini file
            with resources.path('pyef.resources', 'settings.ini') as ini_path:
                path_to_ini_file = str(ini_path)
                Command_Polarization = f"{multiwfn_path} {molden_filename} -set {path_to_ini_file} > {polarization_file}"

            # Check if the atomic polarizations have been computed
            path_to_pol = os.path.join(os.getcwd(), polarization_file)
            backup_path_to_pol = os.path.join(os.getcwd(), "final_optim_polarization.txt")

            xyz_file_path = os.path.join(os.getcwd(), final_structure_file)
            print(f"Looking for atomic multipole calculations here: {path_to_pol}")
           
            num_pt_chgs = 0
            if self.ptChgs:
                ptchg_filename = self.ptChgfp
                df_ptchg = self.getPtChgs(ptchg_filename)
                print(f'Columns in the ptchg dataframe: {df_ptchg.columns}')
                num_pt_chgs = num_pt_chgs + len(df_ptchg['Atom']) 


            #Pick lines to include, can exclude atom indices from calculation by calling function excludeAtomsFromEfieldCalc
            total_lines = Electrostatics.mapcount(xyz_file_path)
            init_all_lines = range(0, total_lines + num_pt_chgs - 2)

            if len(excludeAtoms) > 0:
                to_exclude = excludeAtoms[counter]
                self.excludeAtomsFromEfieldCalc(to_exclude)
                print(f'Excluding atoms: {to_exclude}')
            all_lines = [x for x in init_all_lines if x not in self.excludeAtomfromEcalc]
            print(f'At this point in the calculation I have {len(all_lines)} and am exlcuding: {self.excludeAtomfromEcalc}')
            # Check the contents of the polarization file to see if it finished
            need_to_run_calculation = True
            if os.path.exists(path_to_pol):
                with open(path_to_pol, 'r') as file:
                    contents = file.read()
                    if "Calculation took up" in contents:
                        print(f"   > Atomic Multipoles already calculated here: {f}!!")
                        need_to_run_calculation = False
            elif polarization_scheme == 'Hirshfeld_I' and os.path.exists(backup_path_to_pol):
                path_to_pol = backup_path_to_pol
                with open(path_to_pol, 'r') as file:
                    contents = file.read()
                    if "Calculation took up" in contents:
                        print(f"   > Atomic Multipoles already calculated here: {f}!!")
                        need_to_run_calculation = False

            if need_to_run_calculation:
                print('Atomic Multipole Calculation initialized')
                # Now Run the calculation for atomic dipole and quadrupole moment
                atmrad_src = atmrad_path
                copy_tree(atmrad_src, os.getcwd() + '/atmrad/')
                print(f"   > Submitting Multiwfn job using: {Command_Polarization}")
                proc = subprocess.Popen(Command_Polarization, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                polarization_commands = ['15', '-1'] + self.dict_of_multipole[polarization_scheme] + ['0', 'q']
                proc.communicate("\n".join(polarization_commands).encode())
    
            #if self.makePDB:
            #    pdbName = f + '.pdb'
            #    self.makePDB(xyz_file_path, path_to_pol, 'Multipole', pdbName)
            try:
                # If bond_indices is longer then one, default to manually entry mode
                if bool_manual_mode:
                    file_bond_indices = input_bond_indices[counter]
                    [proj_Efields, bondedAs, bonded_idx, bond_lens, monopole_proj] = self.E_proj_bondIndices(file_bond_indices, xyz_file_path, path_to_pol, all_lines)
                # Otherwise, automatically sense atoms bonded to metal and output E-fields of those
                else:
                    [proj_Efields, bondedAs, bonded_idx, bond_lens, monopole_proj] = self.E_proj_first_coord(atom_idx,xyz_file_path, path_to_pol, all_lines)
                
                results_dict['Max Eproj'] = max(abs(np.array(proj_Efields)))
                results_dict['Projected_Efields V/Angstrom'] = proj_Efields
                results_dict['First order Proj Efield V/A'] = monopole_proj
                results_dict['Bonded Atoms'] = bondedAs
                results_dict['Bonded Indices'] = bonded_idx
                results_dict['Bond Lengths']= bond_lens

            except Exception as e:
                print(repr(e))
                traceback_str = ''.join(traceback.format_tb(e.__traceback__))
                print(traceback_str)
            counter = counter + 1

            # Probably want to add other bonds to this list!
            allspeciesdict.append(results_dict)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"Elapsed time: {elapsed_time:.4f} seconds")
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(f"{Efield_data_filename}.csv")


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
            command_A = f"{multiwfn_path} final_optim.molden"
            results_dir = os.getcwd() + '/'

            results_dict = {}
            results_dict['Name'] = f
            
            for key in charge_types:
                print('Partial Charge Scheme:' + str(key))
                try:
                    full_file_path = os.getcwd() +'/final_optim_' +key+'.txt'
                    path_to_xyz = os.getcwd() + '/final_optim.xyz'
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
                        calc_command = self.dict_of_calcs[key]
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

    def compute_dipole(res_coords, res_chgs):
        '''
           res_coords in units of angstroms, res_chgs are in a.u. units
           returns dipole in units of Debye
        '''
        convert_Debye = 4.80320427 

        #center coordinates
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
        arr_coords = np.array(df[['x', 'y', 'z']])
        arr_chgs = np.array(df['charge'])
        for res_name in res_dict.keys():
            res_labels.append(res_name)
            res_indices = res_dict[res_name]
            res_coords = arr_coords[res_indices]
            res_centroid = res_coords.mean(axis=0)
            res_chgs = arr_chgs[res_indices]
            dip_vec, dip_mag = Electrostatics.compute_dipole(res_coords, res_chgs)
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
            final_structure_file = "final_optim.xyz"
            polarization_file = "Multipole" + polarization_scheme + ".txt"
            molden_filename = "final_optim.molden"

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
                    polarization_commands = ['15', '-1'] + self.dict_of_multipole[polarization_scheme] + ['0', 'q']
                    proc.communicate("\n".join(polarization_commands).encode())
              


                results_dict = {}
                results_dict['Name'] = f


                lst_multipole_dict = Electrostatics.getmultipoles(path_to_pol)
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
                            calc_command = self.dict_of_calcs[key]
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
