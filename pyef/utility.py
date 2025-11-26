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
import basis_set_exchange as bse
from .geometry import Geometry
import json

#file contains all utility functions used for processing and making pyef related decisions

class MoldenObject:
    def __init__(self, xyzfile, moldenfile):
        self.xyzfile = xyzfile
        self.moldenFile = moldenfile
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
             'Np': (237.05, 93, 1.90, 7), 'Pu': (244.06, 94, 1.75, 8), 'Am': (243.06, 95, 1.80, 9)}



    

    def countBasis(self, spherical=True):
        degeneracy = {
        's': 1,
        'p': 3,
        'd': 5 if spherical else 6,
        'f': 7 if spherical else 10,
        'g': 9 if spherical else 15}
        nbf = 0
        in_gto = False
        with open(self.moldenFile) as f:
            for line in f:
                if line.strip().startswith('[GTO]'):
                    in_gto = True
                    continue
                if in_gto:
                    if line.strip().startswith('['):  # end of GTO block
                        break
                    parts = line.split()
                    if parts and parts[0].lower() in degeneracy:
                        nbf += degeneracy[parts[0].lower()]
        return nbf


    def build_core_map(self, family):
        """
        Build a {symbol: core_electrons} dictionary for a given ECP family
        """
        mapping = {}
        for element, Z in self.periodic_table.items():
            try:
                data_str = bse.get_basis(family, elements=[Z], fmt="json")
                # Parse JSON string into dict
                data = json.loads(data_str)
                if str(Z) in data["elements"]:
                    print(f'This is Z: {Z}')
                    mapping[element] = data["elements"][str(Z)]['ecp_electrons']
                else:
                    mapping[element] = 0  # all-electron if no ECP
            except Exception as e:
                print("An error occurred:", e)
                traceback.print_exc()
                mapping[element] = 0  # element not defined
        return mapping

    #utility functions to build properly fix molden format in the presence of ECPs
    def build_hybrid_core_map(self, name, cutoff_Z, heavy_family, base_maps):
        """
        Build core-electron map for hybrid family
        Elements with Z <= cutoff_Z ?~F~R core=0
        Elements with Z > cutoff_Z ?~F~R use heavy_family core mapping
        """
        hybrid_map = {}
        heavy_map = base_maps[heavy_family]
        for element, Z in self.periodic_table.items():
            if Z <= cutoff_Z:
                hybrid_map[element] = 0
            else:
                hybrid_map[element] = heavy_map.get(element, 0)
        return hybrid_map


    def fix_ECPmolden(self, owd, ECP):
        """Prepares output terachem data for analysis, mainly isolating final .xyz frame and naming .molden file appropriotely"""
        FAMILIES = ["lanl2dz", "stuttgart_rsc", "def2", "crenbl"]
        # Hybrid families with cutoff rules
        HYBRID_FAMILIES = {
            "lacvps": {
                "cutoff_Z": 18,          # Z <= 18 ?~F~R all-electron
                "heavy_family": "lanl2dz"},
            "lacvp": {
                "cutoff_Z": 10,          # Z <= 10 ?~F~R all-electron
                "heavy_family": "lanl2dz"
                }
        }

        #check if the ECP json already exists... if not you will need to generate it using the funcitonality above
        json_path = owd + '/' + ECP + "ecp_core_maps.json"
        
        # Check if the file exists
        if os.path.exists(json_path) and os.path.isfile(json_path):
               pass
        #if cannot find the json, will need to look for it:
        else:
            #check if the ECP is part of a hybrid family (in which there is a cutoff in atom numbers below which we will NOT apply ECP)
            if ECP in HYBRID_FAMILIES.keys():
                hname = ECP
                rules = HYBRID_FAMILIES[hname]

                base_maps = {}
                for fam in FAMILIES:
                    print(f"Building core map for {fam}...")
                    base_maps[fam] = self.build_core_map(fam)
                hybrid_maps = {}
                hybrid_maps[hname] = self.build_hybrid_core_map(
                    hname,
                    cutoff_Z=rules["cutoff_Z"],
                    heavy_family=rules["heavy_family"],
                    base_maps=base_maps
                    )
                with open(json_path, "w") as f:
                    json.dump(hybrid_maps, f, indent=2)
            else:
                base_maps = {}
                print(f"Building core map for {ECP}...")
                base_maps[ECP] = self.build_core_map(ECP)
                with open(json_path, "w") as f:
                    json.dump(base_maps, f, indent=2)
                            
        #Now the json should be establish 
    
        with open(json_path, "r") as f:
            print(f'        >loading {ECP} parameters')
            ecp_dict = json.load(f)[ECP]
            change_dict = {}
            #get set of atoms in molecule

            df = Geometry(self.xyzfile).getGeomInfo()
            elems_molden = set(list(df['Atom']))
            for elem in elems_molden:
                if ecp_dict[elem] > 0:
                    #this means we need to change it!
                    Z = self.periodic_table[elem]
                    new_Z = Z - ecp_dict[elem]
                    change_dict[elem] = new_Z
            
        with open(self.moldenFile, 'r') as file:
            content = file.read()
            for change_elem, new_Z in change_dict.items():
                print(f'            Molden file at {self.moldenFile} changed: {change_elem} non-core electrons set to {new_Z}')
                pattern = re.compile(rf'({change_elem}\s+\d+\s+)(\d+)')
                content = pattern.sub(rf'\g<1>{new_Z}', content)
            with open(self.moldenFile, 'w') as file:
                file.write(content)
