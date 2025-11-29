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
#import openbabel
#from biopandas.pdb import PandasPdb
import math
import time
import basis_set_exchange as bse
from .geometry import Geometry
import json

from . import constants

#file contains all utility functions used for processing and making pyef related decisions

class MoldenObject:
    def __init__(self, xyzfile, moldenfile):
        self.xyzfile = xyzfile
        self.moldenFile = moldenfile
        # Use centralized constants from constants module
        self.periodic_table = constants.PERIODIC_TABLE
        self.amassdict = constants.ATOMIC_MASS_DICT



    

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
