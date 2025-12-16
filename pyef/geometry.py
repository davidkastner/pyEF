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

from . import utility as constants  # Backward compatibility alias

class Geometry:

    def __init__(self, filexyz):
        self.xyzfile = filexyz
        # Use centralized constants from constants module
        self.periodic_table = constants.PERIODIC_TABLE
        self.amassdict = constants.ATOMIC_MASS_DICT
        self.vdwdict = constants.VDW_RADII_DICT

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




    def hbond_check(self,donors_strong, hydrogens, acceptors_strong, hconv_cutoff, hcoor_cutoff,  strong_angle_cutoff):
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
        import openbabel
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

        str_Hbonds_solute = self.hbond_check(solute_donors_strong, hydrogens, env_acceptors_strong, hconv_cutoff, hcoor_cutoff,  strong_angle_cutoff)
        str_Hbonds_solvent = self.hbond_check(env_donors_strong, hydrogens, solute_acceptors_strong, hconv_cutoff, hcoor_cutoff, strong_angle_cutoff)


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
        xyzfie = self.xyzfile
        from .analysis import Electrostatics
        # TODO: Implement this method
        pass

    def makePDBcontributionEfield(self):
        #make PDB with b-factor column color by the total E-field contribution from each atom
        xyzfie = self.xyzfile
        from .analysis import Electrostatics
        # TODO: Implement this method
        pass

    def makePDBcontributionESP(self):
        #make PDB with b-factor columned colored by the total contribution form each atom to the ESP
        xyzfie = self.xyzfile
        from .analysis import Electrostatics
        # TODO: Implement this method
        pass

    def compareESP(self):
        #compare the ESP at each atom, show how it differs betwen two differ
        xyzfie = self.xyzfile
        from .analysis import Electrostatics
        # TODO: Implement this method
        pass

    def pdb_esp(self):
        #visualize the ESP of the full system
        xyzfie = self.xyzfile
        from .analysis import Electrostatics
        # TODO: Implement this method
        pass

    def makePDB(self, output_filename, b_col, pdbName):
        ''' Function to generate PDB files with partial charges
        Input: xyzfilename: string of xyz filename
        output_filename: string of output filename
        type_charge: string of type of charge (Monopole or Multipole)
        pdbName: string of output pdb filename
        Output: .pdb file with partial charges
        '''
        import openbabel
        from biopandas.pdb import PandasPdb
        import warnings
        import sys
        import os

        # Suppress Open Babel warnings
        warnings.filterwarnings('ignore', category=DeprecationWarning)

        # Redirect stderr to suppress Open Babel warnings
        stderr_backup = sys.stderr
        sys.stderr = open(os.devnull, 'w')

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

        # IMPORTANT: Handle indexing properly
        # - b_col is 0-indexed (Python convention): b_col[0] is for XYZ atom 0
        # - PDB atoms are 1-indexed (PDB standard): PDB atom 1 corresponds to XYZ atom 0
        # - The DataFrame has rows 0, 1, 2, ... which map to PDB atoms 1, 2, 3, ...
        # - So b_col[i] should go to DataFrame row i, which is PDB atom (i+1), which is XYZ atom i ✓

        # Verify lengths match to catch indexing errors
        if len(b_col) != len(ppdb.df['HETATM']):
            print(f"WARNING: Length mismatch in makePDB!")
            print(f"  b_col length: {len(b_col)}")
            print(f"  PDB DataFrame length: {len(ppdb.df['HETATM'])}")
            print(f"  This may cause incorrect B-factor assignment!")
            # Resize b_col if needed
            if len(b_col) < len(ppdb.df['HETATM']):
                import numpy as np
                b_col_padded = np.zeros(len(ppdb.df['HETATM']))
                b_col_padded[:len(b_col)] = b_col
                b_col = b_col_padded
            else:
                b_col = b_col[:len(ppdb.df['HETATM'])]

        # Explicitly map b_col indices to PDB atom numbers to ensure correct assignment
        # Sort DataFrame by atom_number to ensure proper ordering (just in case)
        ppdb.df['HETATM'] = ppdb.df['HETATM'].sort_values('atom_number').reset_index(drop=True)

        # Now assign: b_col[i] goes to row i, which should be PDB atom (i+1)
        ppdb.df['HETATM']['b_factor'] = b_col
        ppdb.to_pdb(path=pdbName,records=['HETATM'],gz=False,append_newline=True)

        # Restore stderr
        sys.stderr.close()
        sys.stderr = stderr_backup




    #Here we will issue a new pdb with bfactor associatd with the second - first where the
    #final pdb file will contain the corodinates from the first one
    def diff_pdbs(self, initial_pdb, change_pdb, new_pdb_name, new_to_old_map='None'):
        from biopandas.pdb import PandasPdb
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



