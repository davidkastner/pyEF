import numpy as np
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
from importlib import resources
from distutils.dir_util import copy_tree
import openbabel
from biopandas.pdb import PandasPdb
import math
import time

class Geometry:

    def __init__(self, filexyz):
        self.xyzfile = filexyz
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

        #taken from here: https://periodictable.com/Properties/A/VanDerWaalsRadius.v.html

        self.vdwdict = {'X': (1.0, 0, 0.77, 0), 'H':1.2,
             'D': (2.0141, 1, 0.37, 1), 'He': (1.4),
             'Li': (6.94, 3, 1.33, 1), 'Be': (9.0121831, 4, 1.02, 2), 'B': (10.83, 5, 0.85, 3),
             'C': 1.70, 'N': 1.55, 'O':1.52,
             'F':  1.47, 'Ne': (20.1797, 10, 0.67, 8), 'Na': (22.99, 11, 1.55, 1),
             'Mg': (24.30, 12, 1.39, 2), 'Al': (26.98, 13, 1.26, 3), 'Si': (28.08, 14, 1.16, 4),
             'P': 1.75, 'S': 1.8, 'Cl':  1.75,
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

        #Accepts path to the xyz file and returns a dataframe containing the atoms names and the coordinates
    def getGeomInfo(self):
        '''
        Input: 
        filepathtoxyz: string of xyz filename
        Output: dataframe with atomic symbols and coordinates
        '''
        filepathtoxyz  = self.xyzfile
        data = []
        counter_idx = 0
        with open(filepathtoxyz, 'r') as file:

            # Skip the first two lines since they contain meta-deta
            for line in file:
                tokens = line.split()
                #print(tokens)
                if len(tokens) == 4:  # Assuming atom name and x, y, z coordinates are present
                    atom_name = tokens[0].capitalize()
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



            #Accepts path to the xyz file and returns a dataframe containing the atoms names and the coordinates
    def getGeomInfovDW(self):
        '''
        Input:
        filepathtoxyz: string of xyz filename
        Output: dataframe with atomic symbols and coordinates
        '''
        filepathtoxyz  = self.xyzfile
        data = []
        counter_idx = 0
        with open(filepathtoxyz, 'r') as file:
            # Skip the first two lines since they contain meta-deta
            next(file)
            next(file)
            for line in file:
                tokens = line.split()
                if len(tokens) == 4:  # Assuming atom name and x, y, z coordinates are present
                    atom_name = tokens[0].capitalize()
                #for xyz in QMMM simulation, pnt charges at end of file, skip them!
                if atom_name == 'pnt':
                    break
                x, y, z = map(float, tokens[1:])
                rad = self.vdwdict[atom_name]
                data.append([counter_idx, atom_name, x, y, z, rad])
                counter_idx += 1

        columns = ['Index', 'Atom', 'X', 'Y', 'Z', 'Radius']
        df = pd.DataFrame(data, columns=columns)
        # Define upper limit for bond cutoff depending on the two atoms involved
        return df

    def getCentroidDistance(substrate_xyz, ref_xyz):
        centroid_sub = np.mean(np.array(substrate_xyz), axis=1)
        centroid_ref = np.mean(np.array(ref_xyz), axis=1)
        distance = np.linalg.norm(centroid_ref - centroid_sub)
        return distance

    def getSubstrateFramePosition(self, idx_substrate):
        filepathtoxyz  = self.xyzfile
        df_mol = self.getGeomInfo()

        sub_X = df_mol['X'][idx_substrate]
        sub_Y = df_mol['Y'][idx_substrate]
        sub_Z = df_mol['Z'][idx_substrate]

        sub_xyz = np.array([sub_X, sub_Y, sub_Z])
        all_xyz = np.array([df_mol['X'], df_mol['Y'], df_mol['Z']])
        return Geometry.getCentroidDistance(sub_xyz, all_xyz)

    def getBondedAtoms(self,  atomidx):
        '''
        Input:
        filepathtoxyz: string of xyz filename
        atomidx: integer of atom index
        Output: list of integers of bonded atom indices

        '''
        filepathtoxyz  = self.xyzfile
        bonded_atom_indices = []
        df_mol = self.getGeomInfo()
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


    def weak_Hbond_check(self,donors_weak, hydrogens, acceptors_weak, CH_bond_cutoff, HA_bond_cutoff,  weak_angle_cutoff):
        # ---- Weak C?~@~SH...O H-bonds ----
        hbond_weak = []
        for _, donor in donors_weak.iterrows():
            d_pos = np.array(donor['pos'])
            bonded_Hs = hydrogens.copy()
            bonded_Hs['dist'] = bonded_Hs['pos'].apply(lambda h: np.linalg.norm(d_pos - np.array(h)))
            bonded_Hs = bonded_Hs[bonded_Hs['dist'] < CH_bond_cutoff]

            for _, H in bonded_Hs.iterrows():
                h_pos = np.array(H['pos'])
                DH_vec = d_pos - h_pos 
                DH_unit = DH_vec / np.linalg.norm(DH_vec)

                for _, acc in acceptors_weak.iterrows():
                    a_pos = np.array(acc['pos'])
                    DA_vec = a_pos - h_pos
                    Dfull_vec = a_pos - d_pos
                    DA_dist = np.linalg.norm(DA_vec)
                    if DA_dist > HA_bond_cutoff:
                        continue
                    DA_unit = DA_vec / DA_dist
                    angle = np.degrees(np.arccos(np.clip(np.dot(DH_unit, DA_unit), -1.0, 1.0)))
                    if angle >= weak_angle_cutoff:
                        hbond_weak.append({
                            'Type': 'Weak_C-H...O',
                            'Donor': donor['Atom'], 'H': H['Atom'], 'Acceptor': acc['Atom'],
                            'Donor_Index': donor['Index'], 'Acceptor_Index': acc['Index'],
                            'Distance': DA_dist, 'Angle': angle})

        return hbond_weak




    def Hbond_check(self,donors_strong, hydrogens, acceptors_strong, hconv_cutoff, hcoor_cutoff,  strong_angle_cutoff):
        hbond_strong = []
        for _, donor in donors_strong.iterrows():
            d_pos = np.array(donor['pos'])
            bonded_Hs = hydrogens.copy()
            bonded_Hs['dist'] = bonded_Hs['pos'].apply(lambda h: np.linalg.norm(d_pos - np.array(h)))
            bonded_Hs = bonded_Hs[bonded_Hs['dist'] < hconv_cutoff]

            for _, H in bonded_Hs.iterrows():
                h_pos = np.array(H['pos'])
                DH_vec = h_pos - d_pos

                DH_unit = DH_vec / np.linalg.norm(DH_vec)

                for _, acc in acceptors_strong.iterrows():
                    a_pos = np.array(acc['pos'])
                    DA_vec = a_pos - d_pos
                    DA_dist = np.linalg.norm(DA_vec)

                    if DA_dist > hcoor_cutoff:
                        continue
                    DA_unit = DA_vec / DA_dist
                    angle = np.degrees(np.arccos(np.clip(np.dot(DH_unit, DA_unit), -1.0, 1.0)))
                    if angle >= strong_angle_cutoff:
                        hbond_strong.append({
                            'Type': 'Classical',
                            'Donor': donor['Atom'], 'H': H['Atom'], 'Acceptor': acc['Atom'],
                            'Donor_Index': donor['Index'], 'Acceptor_Index': acc['Index'],
                            'Distance': DA_dist, 'Angle': angle })

        return hbond_strong


    def createpdbfromcoords(xyz_coords, charges, atoms, output_filename):
        '''
        Input: xyz_coords: list of tuples of coordinates
        charges: list of charges
        atoms: list of atomic symbols
        output_filename: string of output filename

        Output: .pdb file with partial charges
        '''
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
    def getCloseHcontacts(self, res_atoms):
        print('Looking for close contacts between res atoms: {res_atoms} and environment')
        #definition of a weak Ch..O coordination
        df = self.getGeomInfo()

        # Geometric parameters
        hconv_cutoff = 1.2 #cutoff for h bonded to donor
        hcoor_cutoff = 2.5 #cutoff for h distance from acceptor atom

        CH_bond_cutoff = 1.1
        HA_bond_cutoff = 3.0

        strong_angle_cutoff = 120
        weak_dist_cutoff = 3.2
        weak_angle_cutoff = 100
        bond_cutoff = 1.2
        df['pos'] = df[['X', 'Y', 'Z']].apply(lambda row: np.array([row['X'], row['Y'], row['Z']]), axis=1)
        all_indices = np.arange(0, len(df))
        env_atoms = [val for val in all_indices if val not in res_atoms]
        #CH and O within 3.2
        solute_df = df[df['Index'].isin(res_atoms)]
        env_df = df[df['Index'].isin(env_atoms)]



        # ---- All hydrogens ----
        hydrogens = df[df['Atom'].str.startswith('H')]

        # Results lists

        hbond_strong = []

        ##---- Classical donors and acceptors ----
        donor_atoms_strong = ('O', 'N', 'S', 'F')
        acceptor_atoms_strong = ('O', 'N', 'S', 'F', 'Cl', 'Br', 'I')

        solute_donors_strong = solute_df[solute_df['Atom'].str.startswith(donor_atoms_strong)]
        env_acceptors_strong = env_df[env_df['Atom'].str.startswith(acceptor_atoms_strong)]

        env_donors_strong = env_df[env_df['Atom'].str.startswith(donor_atoms_strong)]
        solute_acceptors_strong = solute_df[solute_df['Atom'].str.startswith(acceptor_atoms_strong)]

        str_Hbonds_solute = self.Hbond_check(solute_donors_strong, hydrogens, env_acceptors_strong, hconv_cutoff, hcoor_cutoff,  strong_angle_cutoff)
        str_Hbonds_solvent = self.Hbond_check(env_donors_strong, hydrogens, solute_acceptors_strong, hconv_cutoff, hcoor_cutoff, strong_angle_cutoff)


        # ---- Weak donors (C–H) and acceptors (include halogens) ----

        donors_weak = solute_df[solute_df['Atom'].str.startswith('C')]
        acceptors_weak = env_df[env_df['Atom'].str.startswith(('O', 'N', 'S', 'F', 'Cl', 'Br', 'I'))]
        hydrogens = df[df['Atom'].str.startswith('H')]

        solv_donors_weak = env_df[env_df['Atom'].str.startswith('C')]
        solute_acceptors_weak = solute_df[solute_df['Atom'].str.startswith(('O', 'N', 'S', 'F', 'Cl', 'Br', 'I'))]

        Hbonds_C_solute = self.weak_Hbond_check(donors_weak, hydrogens, acceptors_weak, CH_bond_cutoff, HA_bond_cutoff, weak_angle_cutoff)
        Hbonds_C_solvent = self.weak_Hbond_check(solv_donors_weak, hydrogens, solute_acceptors_weak, CH_bond_cutoff, HA_bond_cutoff, weak_angle_cutoff)

        # ---- Classical H-bonds (from solute to environment) ----

        hbond_weak = Hbonds_C_solute + Hbonds_C_solvent
        hbond_strong = str_Hbonds_solute + str_Hbonds_solvent

        print(f"🔗 Classical H-bonds from solute → environment: {len(hbond_strong)}")

        for hb in hbond_strong:
            print(f"  {hb['Donor']}({hb['Donor_Index']})–{hb['H']}...{hb['Acceptor']}({hb['Acceptor_Index']}) "
                    f"| dist = {hb['Distance']:.2f} Å, angle = {hb['Angle']:.1f}°")
            print(f"\n🟡 Weak C–H...O interactions: {len(hbond_weak)}")

        print(f"?~_~T~W weak H-bonds from solute ?~F~R environment: {len(hbond_weak)}")
        for hb in hbond_weak:
            print(f"  {hb['Donor']}({hb['Donor_Index']})–{hb['H']}...{hb['Acceptor']}({hb['Acceptor_Index']}) "
                    f"| dist = {hb['Distance']:.2f} Å, angle = {hb['Angle']:.1f}°")


        all_contacts_dict = {'NumHbonds' : len(hbond_strong), 'NumWeakHbonds': len(hbond_weak)}
        return  all_contacts_dict

    def countCloseContacts(self, res_atoms, factor):
        contacts = []
        total_close_atoms = []
        df = self.getGeomInfovDW()
        df['pos'] = df[['X', 'Y', 'Z']].apply(lambda row: np.array([row['X'], row['Y'], row['Z']]), axis=1)
        all_indices = np.arange(0, len(df))
        env_atoms = [val for val in all_indices if val not in res_atoms]

        solute_df = df[df['Index'].isin(res_atoms)]
        env_df = df[df['Index'].isin(env_atoms)]


        for i, atom1 in solute_df.iterrows():
            r1= atom1['Radius']
            pos1 = atom1['pos']

            for j, atom2 in env_df.iterrows():
                pos2 = atom2['pos']
                r2 = atom2['Radius']

                dist = np.linalg.norm(pos1 - pos2)
                contact_threshold = factor * (r1 + r2)


                if dist < contact_threshold :
                    contacts.append({
                    'Atom1_Index': atom1['Index'],
                    'Atom1_Type': atom1['Atom'],
                    'Atom2_Index': atom2['Index'],
                    'Atom2_Type': atom2['Atom'],
                    'Distance': dist,
                    'Threshold': contact_threshold})

        num_contacts = len(contacts)
        df = pd.DataFrame(contacts)

        if num_contacts > 0:
            percent_dist_thrsh = np.array(df['Distance'])/np.array(df['Threshold'])
            avg_percent_dist = np.average(percent_dist_thrsh)

        else:
            avg_percent_dist = 1.1

        new_dict = {f'Contacts_{factor}vdW':num_contacts, 'Avg_dist_vs{factor}perc': avg_percent_dist }
        return new_dict

class Visualize:

    def __init__(self, filexyz='None'):
        self.xyzfile = filexyz

    def makePDBpercentEfield(self):
        #make PDB with b-factor cols colored by percent of final efield 
        xyzfie = self.filexyz
        from .analysis import Electrostatics

    def makePDBcontributionEfield(self):
        #make PDB with b-factor column color by the total E-field contribution from each atom
        xyzfie = self.filexyz


        ppdb.df['HETATM']['b_factor'] = charges
        ppdb.to_pdb(path=pdbName,records=['HETATM'],gz=False,append_newline=True)
        from .analysis import Electrostatics

    def makePDBcontributionESP(self):
        #make PDB with b-factor columned colored by the total contribution form each atom to the ESP
        xyzfie = self.filexyz
        from .analysis import Electrostatics

    def compareESP(self):
        #compare the ESP at each atom, show how it differs betwen two differ
        xyzfie = self.filexyz
        from .analysis import Electrostatics

    def PDBESP():
        #visualize the ESP of the full system
        xyzfie = self.filexyz
        from .analysis import Electrostatics

    def makePDB(self, output_filename, b_col, pdbName):
        ''' Function to generate PDB files with partial charges
        Input: xyzfilename: string of xyz filename
        output_filename: string of output filename
        type_charge: string of type of charge (Monopole or Multipole)
        pdbName: string of output pdb filename
        Output: .pdb file with partial charges
        '''
        xyzfilename = self.xyzfile
        df = Geometry(self.xyzfile).getGeomInfo()
        #get out the xyz coords, charges, and atoms!
        #run the outpute_filename
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("pdb")
        #mol = openbabel.OBMol()
        #df = pd.read_csv(output_filename, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        #atoms = df['Atom']
        #charges = df['charge']
        #xs = df['x']
        #ys = df['y']
        #zs = df['z']

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("xyz", "pdb")
  
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, xyzfilename)
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        mol.FindRingAtomsAndBonds()

        obConversion.WriteFile(mol, pdbName)
        ppdb = PandasPdb()
        ppdb.read_pdb(pdbName)
        ppdb.df['HETATM']['b_factor'] = b_col
        ppdb.to_pdb(path=pdbName,records=['HETATM'],gz=False,append_newline=True)




    #Here we will issue a new pdb with bfactor associatd with the second - first where the
    #final pdb file will contain the corodinates from the first one 
    def diff_pdbs(self, initial_pdb, change_pdb, new_pdb_name, new_to_old_map='None'):
        ppdb = PandasPdb().read_pdb(initial_pdb)
        df_first_pdb = ppdb.df['HETATM']
# Access ATOM/HETATM records as DataFrames

        next_pdb = PandasPdb().read_pdb(change_pdb)
        df_second_pdb = next_pdb.df['HETATM']

        init_indices_tob_map = dict(zip(df_first_pdb['atom_number'], df_first_pdb['b_factor']))
        second_indices_tob_map = dict(zip(df_second_pdb['atom_number'], df_second_pdb['b_factor']))

        new_b_lst = []

        for index, init_bval in init_indices_tob_map.items():
            if new_to_old_map == 'None':
                second_index = index
            else:
                if index in new_to_old_map:
                    second_index = new_to_old_map[index]
                else:
                    new_b_lst.append(0)
                    continue
            if second_index in second_indices_tob_map:
                second_b = second_indices_tob_map[second_index]
                new_b = second_b - init_bval
                new_b_lst.append(new_b)
            else:
                new_b_lst.append(0)

        new_pdb = PandasPdb()
        new_pdb.read_pdb(initial_pdb)
        new_pdb.df['HETATM']['b_factor'] = new_b_lst

        new_pdb.to_pdb(new_pdb_name,records=['HETATM'],gz=False,append_newline=True)

    def makePDBcharges(self, output_filename, b_col, pdbName, type_charge='Monopole'):
        ''' Function to generate PDB files with partial charges
        Input: xyzfilename: string of xyz filename
        output_filename: string of output filename
        type_charge: string of type of charge (Monopole or Multipole)
        pdbName: string of output pdb filename
        Output: .pdb file with partial charges
        '''
        xyzfilename = self.xyzfile
        print(f'Final XYZ filename: {xyzfilename}')
        df = self.getGeomInfo()
        #get out the xyz coords, charges, and atoms!
        #run the outpute_filename
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("pdb")
        mol = openbabel.OBMol()
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
        elif type_charge == 'ChargeFile':
            # Read the file and strip whitespace
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, xyzfilename)
            mol.ConnectTheDots()
            mol.PerceiveBondOrders()
            mol.FindRingAtomsAndBonds

            with open(output_filename) as f:
                charges = [float(line.strip()) for line in f]
            #atom_names = df['Atom']
            #print(f'Here are the atom names: {atom_names}')
            #xs = df['X']
            #ys = df['Y']
            #zs = df['Z']
            #for idx in range(0, len(charges)):
            #    print(f'Index is: {idx}')
            #    #Use the output from the getGeomInfo to carry out the following:
            #    atom = mol.NewAtom()
            #    atom_num = self.periodic_table[atom_names[idx]]
            #    atom.SetAtomicNum(atom_num)
            #    atom.SetVector(xs[idx], ys[idx], zs[idx])
            #    atom.SetPartialCharge(charges[idx])
            #    #atom.GetResidue().SetAtomID(atom, atom_names[idx])


            #print("Unique atom symbols in input:", set(atom_names))
            obConversion.WriteFile(mol, pdbName)
            #print(f"Added {mol.NumAtoms()} atoms to OBMol")
            #print(f"Number of charges: {len(charges)}")
            ppdb = PandasPdb()
            ppdb.read_pdb(pdbName)
            #df_all = pd.concat([
            #    ppdb.df.get('ATOM', pd.DataFrame()),
            #    ppdb.df.get('HETATM', pd.DataFrame())
            #], ignore_index=True)

            # Check consistency
            #if len(charges) != len(df_all):
            #    raise ValueError(f"Length mismatch: {len(charges)} charges vs {len(df_all)} atoms")
            # Assign charges to b_factor column
            #df_all['b_factor'] = charges
            #atom_mask = df_all['record_name'] == 'ATOM'
            #hetatm_mask = df_all['record_name'] == 'HETATM'
            #ppdb.df['ATOM'] = df_all[atom_mask].copy()
            #ppdb.df['HETATM'] = df_all[hetatm_mask].copy()
            #ppdb.to_pdb(path=pdbName,records=['ATOM','HETATM'],gz=False,append_newline=True)
            #print(ppdb.df.columns)
            #print(f"PDB has {len(ppdb.df['HETATM'])} HETATM entries")
            ppdb.df['HETATM']['b_factor'] = charges
            ppdb.to_pdb(path=pdbName,records=['HETATM'],gz=False,append_newline=True)

        elif type_charge == 'Multipole':
            from .analysis import Electrostatics
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



