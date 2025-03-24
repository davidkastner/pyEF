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
from collections import deque
from importlib import resources
from distutils.dir_util import copy_tree
import openbabel
from biopandas.pdb import PandasPdb


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
    inGaCage: boolean
        indicates whether files deal with TMCs or TMCsin Ga Cage
    hasECP: boolean
        indicates if an effective core potential was used... in this case, molden file will need to be re-formatted to be compatible multiwfn!
    '''
    def __init__(self, lst_of_folders, lst_of_tmcm_idx, folder_to_file_path, inGaCage=False, hasECP=False, includePtChgs=False):
        self.lst_of_folders = lst_of_folders
        self.lst_of_tmcm_idx = lst_of_tmcm_idx
        self.folder_to_file_path = folder_to_file_path
        self.dict_of_calcs =  {'Hirshfeld': '1', 'Voronoi':'2', 'Mulliken': '5', 'Lowdin': '6', 'SCPA': '7', 'Becke': '10', 'ADCH': '11', 'CHELPG': '12', 'MK':'13', 'AIM': '14', 'Hirshfeld_I': '15', 'CM5':'16', 'EEM': '17', 'RESP': '18', 'PEOE': '19'}
        self.inGaCageBool = inGaCage
        self.dielectric = 1
        self.ptChgs = includePtChgs

        #defauly setting does not generate PDB files
        self.makePDB = False
        self.excludeAtomfromEcalc = []
        #To avoid over-estimating screening from bound atoms, set dielectric to 1 for primary bound atoms in ESP calv
        self.changeDielectBoundBool = False
        # Dictionary is originally from molsimplify, # Data from http://www.webelements.com/ (last accessed May 13th 2015)
        # Palladium covalent radius seemed to be under-estimated in original implementation, so changed to 1.39 per https://webelements.com/palladium/atom_sizes.html
        # Dictionary from molsimplify, https://molsimplify.readthedocs.io/en/latest/_modules/molSimplify/Classes/globalvars.html
        # Some covalent radii updated with 2008 updated values from webelements.com accessed 1/18/24 
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

    def includePtChgs(self):
        self.ptChgs = True

    def excludeAtomsFromEfieldCalc(self, atom_to_exclude):
        self.excludeAtomfromEcalc = atom_to_exclude

    def minDielecBonds(self, bool_bonds):
        self.changeDielectBoundBool = bool_bonds

    def changeDielectric(self, dlc):
        self.dielectric = dlc

    def makePDB(self):
        self.makePDB = True

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
            with open('final_optim.molden', 'r') as file:
                content = file.read()
            pattern_au = re.compile(r'(Au\s+\d+\s+)(\d+)')
            content = pattern_au.sub(r'\g<1>19', content)

            pattern_fe = re.compile(r'(Fe\s+\d+\s+)(\d+)')
            content = pattern_fe.sub(r'\g<1>8', content)

            pattern_i = re.compile(r'(I\s+\d+\s+)(\d+)')
            content = pattern_i.sub(r'\g<1>7', content)

            with open('final_optim.molden', 'w') as file:
                file.write(content)
            print("      > Molden file is fixed\n")
        os.chdir(owd)

    def prepData(self):
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        owd = os.getcwd()
        print('   > Pre-processing data')

        for f in list_of_folders:
            folder_path = os.path.join(owd, f + folder_to_molden)
            print('      > optim_path: ' + folder_path)

            # Processing optim.xyz to create final_optim.xyz
            optim_file_path = os.path.join(folder_path, 'optim.xyz')
            final_optim_xyz = os.path.join(folder_path, 'final_optim.xyz')

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
                    backup_file = os.path.join(folder_path, 'xyz.xyz')
                    if os.path.exists(backup_file):
                        shutil.copy2(backup_file, final_optim_xyz)
                        logging.info('Single point data found. Using xyz.xyz as fallback.')
                    else:
                        logging.exception(f'An unexpected error occurred while processing optim.xyz: {e}')
            else:
                print(f'      > {final_optim_xyz} already exists.')            

            # Copying .molden files to final_optim.molden
            final_optim_molden = os.path.join(folder_path, 'final_optim.molden')
            
            if not os.path.exists(final_optim_molden):
                try:
                    files = glob.iglob(os.path.join(folder_path, "*.molden"))
                    for file in files:
                        if os.path.abspath(file) != os.path.abspath(final_optim_molden):
                            shutil.copy2(file, final_optim_molden)
                except Exception as e:
                    logging.exception('An Exception was thrown while copying molden files.')
            else:
                print(f"      > {final_optim_molden} already exists.")

        os.chdir(owd)

    def createpdbfromcoords(xyz_coords, charges, atoms, output_filename):
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("pdb")
        mol = openbabel.OBMol()
        for coords, charge, symbol in zip(xyz_coords, charges, atoms):
             x, y, z = coords

             # Create a new atom and set its properties
             atom = mol.NewAtom()
             atom.SetType(symbol)  # Use the atomic symbol, e.g., 'C' for carbon
             atom.SetVector(x, y, z)
             atom.SetPartialCharge(charge)

        # Write the molecule to a PDB file
        obConversion.WriteFile(mol, output_filename)

    def makePDB(self, xyzfilename, output_filename, type_charge, pdbName):
        print(f'Final XYZ filename: {xyzfilename}')
        df = self.getGeomInfo(xyzfilename)
        #get out the xyz coords, charges, and atoms!
        #run the outpute_filename
        obConversion = openbabel.OBConversion()

        if type_charge == 'Monopole':
            df = pd.read_csv(output_filename, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
            atoms = df['Atom']
            charges = df['charge']
            xs = df['x']
            ys = df['y']
            zs = df['z']
            for idx in range(0, len(charges)):
                atom = mol.NewAtom()
                atom.SetType(atoms[idx])
                atom.SetVector(xs[idx], ys[idx], zs[idx])
                atom.SetPartialCharge(charges[idx])
            obConversion.WriteFile(mol, pdbName)
        elif type_charge == 'Multipole':
            charge_lst = []
            lst_multipole_dict = Electrostatics.getmultipoles(output_filename)
            obConversion.SetInAndOutFormats("xyz", "pdb")
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, xyzfilename)

            if len(lst_multipole_dict) != mol.NumAtoms():
                print('Unable to generate .pdb with partial charges')
            else:
                for idx in range(0, len(lst_multipole_dict)):
                    atom = mol.GetAtom(idx + 1)
                    atom_dict = lst_multipole_dict[idx] 
                    charge = atom_dict["Atom_Charge"]
                    atom.SetPartialCharge(charge)
                    charge_lst.append(charge)
                obConversion.WriteFile(mol, pdbName)
                ppdb = PandasPdb()
                ppdb.read_pdb(pdbName)
                ppdb.df['HETATM']['b_factor'] = charge_lst
                ppdb.to_pdb(path=pdbName,records=None,gz=False,append_newline=True)

    #Accepts path to the xyz file and returns a dataframe containing the atoms names and the coordinates
    def getGeomInfo(self, filepathtoxyz):
        data = []
        counter_idx = 0
        with open(filepathtoxyz, 'r') as file:
            # Skip the first two lines since they contain meta-deta
            next(file)
            next(file)
            for line in file:
                
                tokens = line.split()
                if len(tokens) == 4:  # Assuming atom name and x, y, z coordinates are present
                    atom_name = tokens[0]
                #for xyz in QMMM simulation, pnt charges at end of file, skip them!
                if atom_name == 'pnt':
                    break
                x, y, z = map(float, tokens[1:])
                rad = self.amassdict[atom_name][2]
                data.append([counter_idx, atom_name, x, y, z, rad])
                counter_idx += 1

        columns = ['Index', 'Atom', 'X', 'Y', 'Z', 'Radius']
        df = pd.DataFrame(data, columns=columns)
        # Define upper limit for bond cutoff depending on the two atoms involved
        return df

    # Accepts an atom and will determine indices of atoms bound, based on implementation in molsimp
    def getBondedAtoms(self, filepathtoxyz, atomidx):
        bonded_atom_indices = []
        df_mol = self.getGeomInfo(filepathtoxyz)
        atm_rad = df_mol['Radius'][atomidx]
        atm_X = df_mol['X'][atomidx]
        atm_Y = df_mol['Y'][atomidx]
        atm_Z = df_mol['Z'][atomidx]
        df_mol['BondCutoff'] = df_mol['Radius'].apply(lambda y: y + atm_rad)
        distsq = df_mol['X'].apply(lambda x: (x - atm_X)**2) + df_mol['Y'].apply(lambda y: (y - atm_Y)**2) + df_mol['Z'].apply(lambda z: (z - atm_Z)**2)
        dist = distsq.apply(lambda d: d**(1/2))
        for dist_idx in range(0, len(dist)):
            if df_mol['BondCutoff'][dist_idx] > dist[dist_idx]:
                bonded_atom_indices.append(dist_idx)
        bonded_atom_indices.remove(atomidx)
        
        return bonded_atom_indices


    def getmultipoles(multipole_name):
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
        chg_df = pd.read_table(filename_pt, skiprows=2, delim_whitespace=True, names=['charge', 'x', 'y', 'z'])
        atm_name = ['pnt']
        atoms = atm_name*len(chg_df['charge'])
        chg_df['Atom'] = atoms
        return chg_df

    # Define the functions to calculate the ESP:
    def mapcount(filename):
        """Function to rapidly count the number of lines in a file"""

        f = open(filename, "r+")
        buf = mmap.mmap(f.fileno(), 0)
        lines = 0
        readline = buf.readline

        while readline():
            lines += 1
        return lines


    def calcesp(self, path_to_xyz, espatom_idx, charge_range, charge_file):
        """
        Calculate the esp

        Notes
        -----
        Espatom_idx should be the index of the atom at which the esp should be calculated
        Charge_range should provide an array with the index of each atom to be considered in calculation of the ESP
        Charge_file is the filepath for the .txt file containing the hirshfeld charges (generated by multiwfn)

        """
        dielectric = self.dielectric
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2)

        # Convert each column to list for quicker indexing
        atoms = list(df['Atom'])
        charges = list(df['charge'])
        xs = list(df['x'])
        ys = list(df['y'])
        zs = list(df['z'])

        #For QMMM calculation, include point charges in ESP calculation
        if self.ptChgs:
            ptchg_filename = 'ptchrg.xyz'
            init_file_path = path_to_xyz[0:-len('scr/final_optim.xyz')]
            full_ptchg_fp = init_file_path + ptchg_filename
            df_ptchg = self.getPtChgs(full_ptchg_fp)
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
        chargeo = charges[idx_atom]
        total_esp = 0

        # Unit conversion
        A_to_m = 10**(-10)
        KJ_J = 10**-3
        faraday = 23.06   #kcal/(mol*V)
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        cal_J = 4.184

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

        final_esp = k*total_esp*((C_e))*cal_J*faraday   # Note: cal/kcal * kJ/J gives 1
        return [final_esp, df['Atom'][idx_atom]]


    def calc_firstTermE(espatom_idx, charge_range, charge_file):
        # E in units of V/(ansgrom) = N*m/(C*Angstrom)
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

        # Convert each column to list for quicker indexing
        atoms = df['Atom']
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
        chargeo = charges[idx_atom]
        position_vec = [xo, yo, zo]
        Ex = 0
        Ey = 0
        Ez = 0

        # Unit conversion
        A_to_m = 10**(-10)
        KJ_J = 10**-3
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        
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

    def calc_fullE(self, idx_atom, charge_range, xyz_file, atom_multipole_file):

        df = self.getGeomInfo(xyz_file)
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2)

        # Convert each column to list for quicker indexing
        atoms = df['Atom']
        xs = df['X']
        ys = df['Y']
        zs = df['Z']

        # Determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        position_vec = [xo, yo, zo]
        Ex = 0
        Ey = 0
        Ez = 0
        
        # Following derivation of -field strenght from TITAN and TUPÃ
        Shaik_E = np.array([0, 0, 0])

        # Unit conversion
        A_to_m = 10**(-10)
        # Bohr to meters (atomic units)
        b_to_m = 5.291772109*(10**-11)
        b_to_A = 0.529177
        KJ_J = 10**-3
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)


        inv_eps = 1/self.dielectric

        #load multipole moments from processed outputs 
        lst_multipole_dict = Electrostatics.getmultipoles(atom_multipole_file)

        #make a list for each term in multipole expansion
        lst_multipole_idxs = list(range(0, len(lst_multipole_dict)))

        #This ensure that atoms we would like to exclude from Efield calculation are excluded from every term in the multipole expansion
        multipole_chg_range = [i for i in lst_multipole_idxs if i in charge_range]
        for idx in multipole_chg_range:
            #If the atom idx is outside of the charge range then skip
            atom_dict = lst_multipole_dict[idx] 
            if idx == idx_atom:
                continue
            else:
                # Units of dipole_vec: 
                dipole_vec = atom_dict["Dipole_Moment"]

                #convention of vector pointing towards atom of interest, so positive charges exert positive Efield
                dist_vec = np.array([(xs[idx] - xo), (ys[idx] - yo), (zs[idx] - zo)])
                dist_arr = np.outer(dist_vec, dist_vec)
                quadrupole_arr = atom_dict['Quadrupole_Moment']
                
                # Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                Shaik_E = Shaik_E + A_to_m*inv_eps*k*C_e*(-atom_dict["Atom_Charge"]*(1/(r**2))*dist_vec/la.norm(dist_vec))
                Ex_quad = inv_eps*k*C_e*(1/(r**3))*((xs[idx] - xo))*(A_to_m**6)*(1/r**4)*dist_arr[1:, 1:]*quadrupole_arr[1:,1:]
                Ey_quad = inv_eps*k*C_e*(1/(r**3))*((ys[idx] - yo))*(A_to_m**6)*(1/r**4)*dist_arr[0:2:3, 0:2:3]*quadrupole_arr[0:2:3,0:2:3]
                Ez_quad = inv_eps*k*C_e*(1/(r**3))*((zs[idx] - zo))*(A_to_m**6)*(1/r**4)*dist_arr[0:2, 0:2]*quadrupole_arr[0:2,0:2]

                Ex = Ex + inv_eps*k*C_e*(1/(r**3))*((xs[idx] - xo))*(-atom_dict["Atom_Charge"]*(A_to_m**2) + (A_to_m**4)*(1/r**2)*np.dot(dipole_vec[1:], dist_vec[1:]))-(1/3)*Ex_quad.sum()
                Ey = Ey + inv_eps*k*C_e*(1/(r**3))*((ys[idx] - yo))*(-atom_dict["Atom_Charge"]*(A_to_m**2) + (A_to_m**4)*(1/r**2)*np.dot(dipole_vec[0:2:3], dist_vec[0:2:3])) -(1/3)*Ey_quad.sum()
                Ez = Ez + inv_eps*k*C_e*(1/(r**3))*((zs[idx] - zo))*(-atom_dict["Atom_Charge"]*(A_to_m**2) + (A_to_m**4)*(1/r**2)*np.dot(dipole_vec[0:2], dist_vec[0:2])) -(1/3)*Ez_quad.sum()
        E_vec = [Ex, Ey, Ez]

        #For QMMM calculation, include point charges in Shaik E field calc
        if self.ptChgs:
            ptchg_filename = 'ptchrg.xyz'
            init_file_path = xyz_file[0:-len('scr/final_optim.xyz')]
            full_ptchg_fp = init_file_path + ptchg_filename
            df_ptchg = self.getPtChgs(full_ptchg_fp)
            MM_xs = list(df_ptchg['x'])
            MM_ys = list(df_ptchg['y'])
            MM_zs = list(df_ptchg['z'])
            MM_charges = list(df_ptchg['charge'])
            #change charge range to include these new partial charges!
            charge_range = range(0, len(MM_xs))
            for chg_idx in charge_range:
                r = (((MM_xs[chg_idx] - xo)*A_to_m)**2 + ((MM_ys[chg_idx] - yo)*A_to_m)**2 + ((MM_zs[chg_idx] - zo)*A_to_m)**2)**(0.5)
                dist_vec = np.array([(MM_xs[chg_idx] - xo), (MM_ys[chg_idx] - yo), (MM_zs[chg_idx] - zo)])
                Shaik_E = Shaik_E +  A_to_m*inv_eps*k*C_e*(-MM_charges[chg_idx]*(1/(r**2))*dist_vec/la.norm(dist_vec))
            #Add contributions to Shaik Efield from point charges 
           
        return [E_vec, position_vec, df['Atom'][idx_atom], Shaik_E]
    
    def ESPfromMultipole(self, xyfilepath, atom_multipole_file, charge_range, idx_atom):
        df = self.getGeomInfo(xyfilepath)
        #df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2)

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
        KJ_J = 10**-3
        faraday = 23.06   #kcal/(mol*V)
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        cal_J = 4.184
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
                #now account for bound atoms
                 r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                 total_esp = total_esp + (atom_dict["Atom_Charge"]/r)
            else:
                # Calculate esp and convert to units (A to m)
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                total_esp = total_esp + (1/dielectric)*(atom_dict["Atom_Charge"]/r)

        final_esp = k*total_esp*((C_e))*cal_J*faraday   # Note: cal/kcal * kJ/J gives 1
        return [final_esp, df['Atom'][idx_atom] ]
 

    # Bond_indices is a list of tuples where each tuple contains the zero-indexed values of location of the atoms of interest
    def E_proj_bondIndices(self, bond_indices, xyz_filepath, atom_multipole_file, all_lines):
        bonded_atoms = []
        E_projected = []
        E_shaik_proj = []
        bonded_positions = []
        # Determine the Efield vector at point of central metal stom
        bond_lens = []
        print(f'Here are the bond indices: {bond_indices}')
        for atomidxA, atomidxB in bond_indices:
            [A_bonded_E, A_bonded_position, A_bonded_atom, A_Shaik_E_bonded]  =  self.calc_fullE(atomidxA, all_lines, xyz_filepath, atom_multipole_file)   
            [B_bonded_E, B_bonded_position, B_bonded_atom, B_Shaik_E_bonded]  =  self.calc_fullE(atomidxB, all_lines, xyz_filepath, atom_multipole_file)  
            bond_vec_unnorm = np.subtract(np.array(A_bonded_position), np.array(B_bonded_position)) 
            bond_len = np.linalg.norm(bond_vec_unnorm)
            bond_vec = bond_vec_unnorm/(bond_len)

            # Initialized a bond_dipole_vec as the (bond_vec_unnorm )*(sum of the partial charges).. can just use dipole! 
            # Compute E-field projected along this bond!
            E_proj = (1/2)*np.dot((np.array(A_bonded_E) + np.array(B_bonded_E)), bond_vec)
            E_proj_Shaik = (1/2)*np.dot((A_Shaik_E_bonded + B_Shaik_E_bonded), bond_vec)
            E_projected.append(E_proj)
            E_shaik_proj.append(E_proj_Shaik)
            bonded_atoms.append((A_bonded_atom, B_bonded_atom))
            bonded_positions.append((A_bonded_position, B_bonded_position))
            bond_lens.append(bond_len)
        return [E_projected, bonded_atoms, bond_indices, bond_lens, E_shaik_proj]
    
    #Calculate Efield projection accounting only for electrostatic effects of directly bound atoms (proxy for inductive effect)
    def E_proj_first_coord(self, metal_idx, xyz_file_path, atom_multipole_file, all_lines):
        bonded_atoms = []
        E_projected = []
        E_shaik_proj = []
        bonded_positions = []
        # Determine the Efield vector at point of central metal stom
        [center_E, center_position, center_atom, Shaik_E_center]  =  self.calc_fullE(metal_idx, all_lines, xyz_file_path, atom_multipole_file)
        lst_bonded_atoms = self.getBondedAtoms(xyz_file_path, metal_idx) 
        
        bond_lens = []
        for bonded_atom_idx in lst_bonded_atoms:
            [bonded_E, bonded_position, bonded_atom, Shaik_E_bonded]  =  self.calc_fullE(bonded_atom_idx, all_lines, xyz_file_path, atom_multipole_file)    
            bond_vec_unnorm = np.subtract(np.array(center_position), np.array(bonded_position)) 
            bond_len = np.linalg.norm(bond_vec_unnorm)
            bond_vec = bond_vec_unnorm/(bond_len)

            # Initialized a bond_dipole_vec as the (bond_vec_unnorm )*(sum of the partial charges).. can just use dipole! 
            # Compute E-field projected along this bond!
            E_proj = (1/2)*np.dot((np.array(bonded_E) + np.array(center_E)), bond_vec)
            E_proj_Shaik = (1/2)*np.dot((Shaik_E_center + Shaik_E_bonded), bond_vec)
            E_projected.append(E_proj)
            E_shaik_proj.append(E_proj_Shaik)
            bonded_atoms.append(bonded_atom)
            bonded_positions.append(bonded_position)
            bond_lens.append(bond_len)
        return [E_projected, bonded_atoms, lst_bonded_atoms, bond_lens, E_shaik_proj]
            
    #Calc ESP accounting only for electrostatic contributions for atoms bound to ESP center
    def esp_first_coord(self, metal_idx, charge_file, path_to_xyz):
        print('The index of the metal atom is: ' + str(metal_idx))
        lst_bonded_atoms = self.getBondedAtoms(path_to_xyz, metal_idx)
        [First_coord_ESP, atom_type] = self.calcesp(path_to_xyz, metal_idx, lst_bonded_atoms, charge_file)
        return First_coord_ESP
    #Calc ESP accounting for first and second coordinating spheres of atom the ESP center atom
    def esp_second_coord(self, metal_idx, charge_file, path_to_xyz):
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

    # Boolean CageTrue
    def charge_atom(filename, atom_idx):
        df = pd.read_csv(filename, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']
        total_charge = np.sum(df['charge'])
        partial_charge_atom = charges[atom_idx]
        return [total_charge, partial_charge_atom]

    def getAtomInfo(filename, atom_idx):
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


    def ESP_all_calcs(self, path_to_xyz, filename, atom_idx, cageTrue):
        # Get the number of lines in the txt file
        dlc = self.dielectric
        total_lines =Electrostatics.mapcount(filename)
        all_lines = range(0, total_lines)
        [ESP_all, atom_type] = self.calcesp(path_to_xyz, atom_idx, all_lines, filename)
        if cageTrue:
            cage_lines = range(0, 280)
            guest_lines = range(280, total_lines)
            # print('Cage indices: ' + str(cage_lines))
            # print('Guest indices: ' + str(guest_lines))
            ESP_just_ligand = self.calcesp(path_to_xyz, atom_idx, guest_lines, filename)[0]
            ESP_just_cage = self.calcesp(path_to_xyz, atom_idx, cage_lines, filename)[0]
            print('ESP for all atoms: ' + str(ESP_all) + ' kJ/(mol*e)')
            print('ESP just ligand: ' + str(ESP_just_ligand) + ' kJ/(mol*e)')
            print('ESP just cage: ' + str(ESP_just_cage) + ' kJ/(mol*e)')
            return [ESP_all, ESP_just_ligand, ESP_just_cage, atom_type]
        else:
            print('ESP at '+str(atom_type) + ' (index = ' + str(atom_idx) + ') : ' + str(ESP_all) + ' kJ/(mol*e)') 
            return [ESP_all, atom_type]


    def esp_bydistance(self, path_to_xyz, espatom_idx,  charge_file):
        dielectric = self.dielectric
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  # Coulombic constant in kg*m**3/(s**4*A**2)

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
            ptchg_filename = 'ptchrg.xyz'
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
                esps.append(k*C_e*cal_J*faraday*charges[idx]/r)
                distances.append(r)
            else:
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                distances.append(r)
                esps.append(k*(1/dielectric)*C_e*cal_J*faraday*charges[idx]/r)
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
       
        new_bool_still = idx_atom in sorted_idx
        #final_sorted_idx = np.delete(init_sorted_idx, init_sorted_idx[idx_atom])
        return [sorted_dist, sorted_esps, cumulative_esps, sorted_idx, sorted_partial_charges, sorted_atomTypes]
  

    # list_of_folders = the list of the folders that contain the desired files
    # new_dir: the [post-folder path to the scr folder that contains the .molden and optim.xyz file themselfs
    # dict of calcs, calculations to be performed by multiwavefunction with the corresponding keys
    # newfilanme: desired name of the .csv fiole that will be createcd in getData cotnaining all of the ESP/other data extracted un the file

    def getESPData(self, charge_types, ESPdata_filename, multiwfn_module, multiwfn_path, atmrad_path, dielectric=1):
        '''
        Function computes a series of ESP data using the charge scheme specified in charge types.

        Attributes
        ----------
        charge_types: list of strings
        ESPdata_filename: string
            Name of the output file name

        '''
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
                        if self.inGaCageBool:
                            # With newly analyzed partial charges, re-compute ESP data
                            [ESP_all, ESP_just_ligand, ESP_just_cage, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx, self.inGaCageBool)
                        else:
                            [ESP_all, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx, self.inGaCageBool)

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
                        
                        if self.inGaCageBool:                        
                            # With newly analyzed partial charges, re-compute ESP data
                            [ESP_all, ESP_just_ligand, ESP_just_cage, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx, self.inGaCageBool)
                        else:
                            [ESP_all, atom_type] = self.ESP_all_calcs(path_to_xyz, full_file_path, atom_idx, self.inGaCageBool)
                        
                        [total_charge,partial_charge_atom] = Electrostatics.charge_atom(full_file_path, atom_idx)
                        [sorted_distances, sorted_esps, cum_esps, sorted_cum_idx, sorted_cum_chg, sorted_atomTypes] = self.esp_bydistance(path_to_xyz, atom_idx, full_file_path)
                        ESP_fcoord = self.esp_first_coord(atom_idx, full_file_path, path_to_xyz)
                        ESP_scoord = self.esp_second_coord(atom_idx, full_file_path, path_to_xyz)

                    # At this point, all calculations shouldbe complete and succesfull: Add ESP data to dictionary
                    results_dict[str(key) + ' ESP Second Coor Shell (kcal/mol)'] = ESP_scoord
                    results_dict[str(key) + ' ESP First Coor Shell (kcal/mol)'] = ESP_fcoord
                    if self.inGaCageBool:
                        results_dict['ESP_ligand '+ str(key)] = ESP_just_ligand
                        results_dict['ESP_just_cage ' +str(key)] = ESP_just_cage
                    else:
                        pass
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

                    # If .molden files deal with encapsulated TMCS, complete an additional set of analyses
                    if self.inGaCageBool:
                        [distoGa, Ga_selfdist]=Electrostatics.calcdist(atom_idx, full_file_path)
                        [cageAtomdist]=Electrostatics.calcNearestCageAtom(atom_idx, full_file_path)
                        results_dict['MetaltoGa_dist'] = distoGa
                        results_dict['GatoGa_dist'] = Ga_selfdist
                        results_dict['NearestCageAtom'] = cageAtomdist
                    else:
                        continue
                except Exception as e:
                    logging.exception('An Exception was thrown')
                    continue
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(ESPdata_filename +'.csv')
        return df
    
    def getESPLargeData(self, ESP_filename, multiwfn_module, multiwfn_path):
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # Old working directory
        allspeciesdict = []
        counter = 0
        
        for f in list_of_file:
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
            polarization_file = "final_optim_polarization.txt"
            
            # Dynamically get path to package settings.ini file
            with resources.path('pyef.resources', 'settings.ini') as ini_path:
                path_to_ini_file = str(ini_path)
                Command_Polarization = f"{multiwfn_path} {molden_filename} -set {path_to_ini_file} > {polarization_file}"

            # Check if the atomic polarizations have been computed
            path_to_pol = os.path.join(os.getcwd(), polarization_file)
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
                polarization_commands = ['15', '-1', '3', '2', '0', 'q'] # For atomic dipole and quadrupole moment of Hirshfeld-I type
                proc.communicate("\n".join(polarization_commands).encode())
            try:
                [final_esp, atom_name] = self.ESPfromMultipole(xyz_file_path, path_to_pol, all_lines, atom_idx)
                results_dict['Total ESP'] = final_esp
                results_dict['Atom'] = atom_name 
            except Exception as e:
                print(repr(e))
                traceback_str = ''.join(traceback.format_tb(e.__traceback__))
                print(traceback_str)
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(f"{ESP_filename}.csv")



    # input_bond_indices is a list of a list of tuples
    def getEFieldData(self, Efield_data_filename, multiwfn_module, multiwfn_path, atmrad_path, input_bond_indices=[], excludeAtoms=[]):

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
            polarization_file = "final_optim_polarization.txt"
            
            # Dynamically get path to package settings.ini file
            with resources.path('pyef.resources', 'settings.ini') as ini_path:
                path_to_ini_file = str(ini_path)
                Command_Polarization = f"{multiwfn_path} {molden_filename} -set {path_to_ini_file} > {polarization_file}"

            # Check if the atomic polarizations have been computed
            path_to_pol = os.path.join(os.getcwd(), polarization_file)
            xyz_file_path = os.path.join(os.getcwd(), final_structure_file)
            print(f"Attempting polarization file path: {path_to_pol}")
            
            #Pick lines to include, can exclude atom indices from calculation by calling function excludeAtomsFromEfieldCalc
            total_lines = Electrostatics.mapcount(xyz_file_path)
            init_all_lines = range(0, total_lines - 2)
            print(f'Initially I have {len(init_all_lines)}')
            all_lines = [x for x in init_all_lines if x not in self.excludeAtomfromEcalc]
            print(f'At this point in the calculation I have {len(all_lines)} and am exlcuding: {self.excludeAtomfromEcalc}')
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
                atmrad_src = atmrad_path
                copy_tree(atmrad_src, os.getcwd() + '/atmrad/')
                print(f"   > Submitting Multiwfn job using: {Command_Polarization}")
                proc = subprocess.Popen(Command_Polarization, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                polarization_commands = ['15', '-1', '4', '1', '2','0', 'q'] # For atomic dipole and quadrupole moment of Hirshfeld-I type
                proc.communicate("\n".join(polarization_commands).encode())
    
            if self.makePDB:
                pdbName = f + '.pdb'
                self.makePDB(xyz_file_path, path_to_pol, 'Multipole', pdbName)
            try:
                # If bond_indices is longer then one, default to manually entry mode
                if bool_manual_mode:
                    print(f' IN BOOL MANUAL MODE')
                    file_bond_indices = input_bond_indices[counter]
                    if len(excludeAtoms) > 0:
                        to_exclude = excludeAtoms[counter]
                        self.excludeAtomsFromEfieldCalc(to_exclude)
                        print(f'Ecluding atoms:m {to_exclude}')
                    [proj_Efields, bondedAs, bonded_idx, bond_lens, shaik_proj] = self.E_proj_bondIndices(file_bond_indices, xyz_file_path, path_to_pol, all_lines)
                # Otherwise, automatically sense atoms bonded to metal and output E-fields of those
                else:
                    [proj_Efields, bondedAs, bonded_idx, bond_lens, shaik_proj] = self.E_proj_first_coord(atom_idx,xyz_file_path, path_to_pol, all_lines)
                
                results_dict['Max Eproj'] = max(abs(np.array(proj_Efields)))
                results_dict['Projected_Efields V/Angstrom'] = proj_Efields
                results_dict['Shaik proj Efields V/A'] = shaik_proj
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

        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(f"{Efield_data_filename}.csv")


    def getpartialchgs(self, charge_types, lst_atom_idxs, partial_chg_filename, multiwfn_path, multiwfn_module, atmrad_path):
        '''
        Function computes partial charges on a select set of atoms using the charge scheme specified in charge types. Note atom indices will be carried over between csvs
        
        Attributes
        ----------
        charge_types: list(str)
            list of strings
        lst_of_atom_idxs: list(int)
            list of integers denoting atom indices (0 indexed!)
        partial_chg_filename: string
            Name of the output file name
        
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

    def getcharge_residues(self, charge_types, res_dict, partial_chg_filename, multiwfn_path, multiwfn_module, atmrad_path):
        '''
        Function computes partial charges on a select set of atoms using the charge scheme specified in charge types. Note atom indices will be carried over between csvs

        Attributes
        ----------
        charge_types: list(str)
            list of strings
        res_dict: dictionary with strings mapped to list(int)
            strings denote name of residues which are mapped to the associated atom indices in the lists(0 indexed!)
        partial_chg_filename: string
            Name of the output file name

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
                            res_indices = res_dict[res_name]
                            [res_charge] = Electrostatics.getAtomsInfo(full_file_path, res_indices)
                            results_dict[f'{res_name} Charge'] = res_charge

                except Exception as e:
                    logging.exception('An Exception was thrown')
                    continue
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(partial_chg_filename +'.csv')
        return df

