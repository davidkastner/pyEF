from distutils.dir_util import copy_tree
from collections import deque
import traceback
import logging
import subprocess
import sys
import glob, os, shutil
import numpy as np
import pandas as pd
import mmap
from GeometryCheck import ErrorAnalysis
from molSimplify.Classes.ligand import ligand_breakdown
from molSimplify.Classes.mol3D import *
from molSimplify.Scripts import *
from collections import deque
from molSimplify.Classes.globalvars import all_angle_refs
import scipy.linalg as la

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
    def __init__(self, lst_of_folders, lst_of_tmcm_idx, folder_to_file_path, inGaCage, hasECP=False):
        self.lst_of_folders = lst_of_folders
        self.lst_of_tmcm_idx = lst_of_tmcm_idx
        self.folder_to_file_path = folder_to_file_path
        self.dict_of_calcs =  {'Hirshfeld': '1', 'Voronoi':'2', 'Mulliken': '5', 'Lowdin': '6', 'SCPA': '7', 'Becke': '10', 'ADCH': '11', 'CHELPG': '12', 'MK':'13', 'AIM': '14', 'Hirshfeld_I': '15', 'CM5':'16', 'EEM': '17', 'RESP': '18', 'PEOE': '19'}
        self.inGaCageBool = inGaCage
        self.prepData()

    #function prepares output terachem data for analysis, mainly isolating final .xyz frame and naming .molden file appropriotely
    
    def fix_ECPmolden(self):
        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        owd = os.getcwd()
        print('Re-formatting .molden file to fix ECP artifacts')
        for f in list_of_folders:
            os.chdir(owd)
            print('I am going to: ' + str(f + folder_to_molden))
            os.chdir(f + folder_to_molden)
            with open('final_optim.molden', 'r') as file:
                content = file.read()
            pattern_au = re.compile(r'(Au\s+\d+\s+)(\d+)')
            content = pattern_au.sub(r'\g<1>19', content)

            pattern_i = re.compile(r'(I\s+\d+\s+)(\d+)')
            content = pattern_i.sub(r'\g<1>7', content)

            with open('final_optim.molden', 'w') as file:
                file.write(content)
            print("all molden files are now fixed")
        os.chdir(owd)
    def prepData(self):
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_folders = self.lst_of_folders
        owd = os.getcwd()
        print('Pre-processing data')
        for f in list_of_folders:
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            optim_path = os.getcwd() + '/'
            print('optim_path: ' + str(optim_path))            

            if os.path.exists('final_optim.xyz'):
                pass  #file already parsed
            else:
                try:
                    #first copy the last xyz frame of the optimization into a new file with standard name
                    path_to_init = optim_path[:-4]
                    os.chdir(optim_path)
                    optim_file = 'optim.xyz'
                    full_traj = open(optim_file, 'r')
                    print('At optim path: ' + str(optim_path))
                    num_atoms = int(full_traj.readline())
                    num_lines = num_atoms + 2
                    with open(optim_file) as input_file:
                        head = [next(input_file) for _ in range(num_lines)]
                    with open('initial_' + optim_file, 'w') as initxyz:
                        initxyz.writelines(head)
                    #last xyz in trajectory saved as the final file
                    with open('final_'+optim_file, 'w') as finalxyz:
                        finalxyz.writelines(deque(full_traj, num_lines))
                except Exception as e:
                    shutil.copy2(os.path.join(optim_path, 'xyz.xyz'), 'final_optim.xyz')
                    logging.exception('An Exception was thrown')

            if os.path.exists('final_optim.molden'):
                pass #.molden file already has correct name
            else:
                try:    
                   #now make a copy of the .molden file with the same name is the directory
                    files = glob.iglob(os.path.join(optim_path, "*.molden"))
                    for file in files:
                        shutil.copy2(file, 'final_optim.molden')

                except Exception as e:
                    logging.exception('An Exception was thrown')
          

        os.chdir(owd)
     
    def getmultipoles(multipole_name):
        # Read the text file
        with open(multipole_name, 'r') as file:
            text = file.read()
            # Split the text into sections based on the delimiter
            sections = re.split(r'\n\s*[*]+\s*', text)

        # Define a pattern to extract the index, element name, and atomic dipole and quadrupole moment values
        pattern = r'Atomic multipole moments of\s+(\d+)\s*\(([^)]+)\s*\).*?X=\s*([-+]?\d+\.\d+)\s+Y=\s*([-+]?\d+\.\d+)\s+Z=\s*([-+]?\d+\.\d+).*?XX=\s*([-+]?\d+\.\d+)\s+XY=\s*([-+]?\d+\.\d+)\s+XZ=\s*([-+]?\d+\.\d+).*?YX=\s*([-+]?\d+\.\d+)\s+YY=\s*([-+]?\d+\.\d+)\s+YZ=\s*([-+]?\d+\.\d+).*?ZX=\s*([-+]?\d+\.\d+)\s+ZY=\s*([-+]?\d+\.\d+)\s+ZZ=\s*([-+]?\d+\.\d+)'

        atomicDicts = []
        # Iterate through sections and extract the information
        for section in sections:
            match = re.search(pattern, section, re.DOTALL)
            if match:
                index = match.group(1)
                element = match.group(2)
                x_value = match.group(3)
                y_value = match.group(4)
                z_value = match.group(5)
                dipole_moment = [float(x_value), float(y_value), float(z_value)]

                xx_value = match.group(6)
                xy_value = match.group(7)
                xz_value = match.group(8)
                yx_value = match.group(9)
                yy_value = match.group(10)
                yz_value = match.group(11)
                zx_value = match.group(12)
                zy_value = match.group(13)
                zz_value = match.group(14)

                # Create a 3x3 matrix for the quadrupole moment
                quadrupole_moment = np.array([
                [float(xx_value), float(xy_value), float(xz_value)],
                [float(yx_value), float(yy_value), float(yz_value)],
                [float(zx_value), float(zy_value), float(zz_value)]
                ])

                atomDict = {"Index": index, "Element": element, 'Dipole_Moment': dipole_moment, 'Quadrupole_Moment': quadrupole_moment}
                atomicDicts.append(atomDict)
        return atomicDicts


    #Define the functions to calculate the ESP:
    #function to rapidly count the number of lines in a file
    def mapcount(filename):
        f = open(filename, "r+")
        buf = mmap.mmap(f.fileno(), 0)
        lines = 0
        readline = buf.readline
        while readline():
            lines += 1
        return lines
    def calcdist(espatom_idx, charge_file):
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

        #convert each column to list for quicker indexing
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']

        #pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        #determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        chargeo = charges[idx_atom]
        total_esp = 0
        #distance to gallium atoms
        Ga_idxs = [276, 277, 278, 279]
        Ga_self_dist = np.zeros([4,4])
        distances_vertices = []
        counter = 0
        for ga_idx in Ga_idxs:
            Gax = xs[ga_idx]
            Gay = ys[ga_idx]
            Gaz = zs[ga_idx]
            r = (((Gax - xo))**2 + ((Gay - yo))**2 + ((Gaz - zo))**2)**(0.5)
            distances_vertices.append(r)
            for new_counter in range(0, 3):
                if new_counter == counter:
                    continue
                else:
                    newGa_idx = counter + 276
                    GatoGadist  = (((Gax - xs[newGa_idx]))**2 + ((Gay - ys[newGa_idx]))**2 + ((Gaz - zs[newGa_idx]))**2)**(0.5)
                    Ga_self_dist[counter, new_counter] = GatoGadist
            counter = counter + 1
            
        return [distances_vertices, Ga_self_dist]

    #espatom_idx should be the index of the atom at which the esp should be calculated
    #charge_range should provide an array with the index of each atom to be considered in calculation of the ESP
    #charge_file is the filepath for the .txt file containing the hirshfeld charges (generated by multiwfn)
    def calcesp(espatom_idx, charge_range, charge_file):
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

        #convert each column to list for quicker indexing
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']

        #pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        #determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        chargeo = charges[idx_atom]
        total_esp = 0

        #unit conversion
        A_to_m = 10**(-10)
        KJ_J = 10**-3
        faraday = 23.06   #kcal/(mol*V)
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        cal_J = 4.184

        for idx in charge_range:
            if idx == idx_atom:
                continue
            else:
                #Calculate esp and convert to units (A to m)
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                total_esp = total_esp + (charges[idx]/r)

        final_esp = k*total_esp*((C_e))*cal_J*faraday   #note that cal/kcal * kJ/J gives 1
        #print(str(final_esp) + ' kJ/(mol*e)')
        return [final_esp, df['Atom'][idx_atom]]


    #E in units of V/(ansgrom) = N*m/(C*Angstrom)
    def calc_firstTermE(espatom_idx, charge_range, charge_file):
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

        #convert each column to list for quicker indexing
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']

        #pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        #determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        chargeo = charges[idx_atom]
        position_vec = [xo, yo, zo]
        Ex = 0
        Ey = 0
        Ez = 0

        #unit conversion
        A_to_m = 10**(-10)
        KJ_J = 10**-3
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        
        print('This is the charge range...')
        print(charge_range)
        for idx in charge_range:
            if idx == idx_atom:
                continue
            else:
                #Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                Ex = Ex - k*C_e*(charges[idx]/r)*(1/(xs[idx] - xo))   
                Ey = Ey - k*C_e*(charges[idx]/r)*(1/(ys[idx] - yo))
                Ez = Ez - k*C_e*(charges[idx]/r)*(1/(zs[idx] - zo))

        E_vec = [Ex, Ey, Ez]
        return [E_vec, position_vec, df['Atom'][idx_atom]]

    def calc_fullE(espatom_idx, charge_range, charge_file, atom_multipole_file):
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

        #convert each column to list for quicker indexing
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']

        #pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        #determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        chargeo = charges[idx_atom]
        position_vec = [xo, yo, zo]
        Ex = 0
        Ey = 0
        Ez = 0
        
        #Following derivation of -field strenght from TITAN and TUPÃ
        Shaik_E = np.array([0, 0, 0])

        #unit conversion
        A_to_m = 10**(-10)
        #bohr to meters (atomic units)
        b_to_m = 5.291772109*(10**-11)
        b_to_A = 0.529177
        KJ_J = 10**-3
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        inv_eps = 1/5
        lst_multipole_dict = Electrostatics.getmultipoles(atom_multipole_file)
        for idx in charge_range:
            atom_dict = lst_multipole_dict[idx] 
            if idx == idx_atom:
                continue
            else:
                #units of dipole_vec: 
                dipole_vec = atom_dict["Dipole_Moment"]
                inv_dist_vec = np.array([1/(xs[idx] - xo),1/(ys[idx] - yo), 1/(zs[idx] - zo)])
                dist_vec = np.array([(xs[idx] - xo), (ys[idx] - yo), (zs[idx] - zo)])
                dist_arr = np.outer(dist_vec, dist_vec)
                quadrupole_arr = atom_dict['Quadrupole_Moment']
                #r_2A = ((xs[idx] - xo))**2 + ((ys[idx] - yo))**2 + ((zs[idx] - zo))**2
                #Shaik_E = Shaik_E + inv_eps*k*C_e*charges[idx]*(1/r_2A)*dist_vec 
                #Calculate esp and convert to units (A to m); Calc E-field stenth in kJ/mol*e*m
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                Shaik_E = Shaik_E + A_to_m*inv_eps*k*C_e*charges[idx]*(1/(r**2))*dist_vec/la.norm(dist_vec)
                #Ex = Ex + inv_eps*k*C_e*(1/(r**3))*((xs[idx] - xo))*(-charges[idx]*(A_to_m**2) + (1/r**2)*A_to_m*np.dot(dipole_vec, dist_vec))
                #Ey = Ey + inv_eps*k*C_e*(1/(r**3))*((ys[idx] - yo))*(-charges[idx]*(A_to_m**2) + (1/r**2)*A_to_m*np.dot(dipole_vec, dist_vec))
                #Ez = Ez + inv_eps*k*C_e*(1/(r**3))*((zs[idx] - zo))*(-charges[idx]*(A_to_m**2) + (1/r**2)*A_to_m*np.dot(dipole_vec, dist_vec))
                Ex_quad = inv_eps*k*C_e*(1/(r**3))*((xs[idx] - xo))*(A_to_m**6)*(1/r**4)*dist_arr[1:, 1:]*quadrupole_arr[1:,1:]
                Ey_quad = inv_eps*k*C_e*(1/(r**3))*((ys[idx] - yo))*(A_to_m**6)*(1/r**4)*dist_arr[0:2:3, 0:2:3]*quadrupole_arr[0:2:3,0:2:3]
                Ez_quad = inv_eps*k*C_e*(1/(r**3))*((zs[idx] - zo))*(A_to_m**6)*(1/r**4)*dist_arr[0:2, 0:2]*quadrupole_arr[0:2,0:2]

                Ex = Ex + inv_eps*k*C_e*(1/(r**3))*((xs[idx] - xo))*(-charges[idx]*(A_to_m**2) + (A_to_m**4)*(1/r**2)*np.dot(dipole_vec[1:], dist_vec[1:]))-(1/3)*Ex_quad.sum()
                Ey = Ey + inv_eps*k*C_e*(1/(r**3))*((ys[idx] - yo))*(-charges[idx]*(A_to_m**2) + (A_to_m**4)*(1/r**2)*np.dot(dipole_vec[0:2:3], dist_vec[0:2:3])) -(1/3)*Ey_quad.sum()
                Ez = Ez + inv_eps*k*C_e*(1/(r**3))*((zs[idx] - zo))*(-charges[idx]*(A_to_m**2) + (A_to_m**4)*(1/r**2)*np.dot(dipole_vec[0:2], dist_vec[0:2])) -(1/3)*Ez_quad.sum()
        E_vec = [Ex, Ey, Ez]
        return [E_vec, position_vec, df['Atom'][idx_atom], Shaik_E]

    def E_proj_first_coord(mol_simp_obj, metal_idx, charge_file, atom_multipole_file):
        ang_to_m = 10**(-10)
        bonded_atoms = []
        E_projected = []
        E_shaik_proj = []
        bonded_positions = []
        total_lines = Electrostatics.mapcount(charge_file)
        all_lines = range(0, total_lines)
        #Determine the Efield vector at point of central metal stom
        [center_E, center_position, center_atom, Shaik_E_center]  =  Electrostatics.calc_fullE(metal_idx, all_lines, charge_file, atom_multipole_file)
        print('Computed E field!')
        print('Efield on metal atom is: ' +str(center_E))
        print('Shaik Efield on metal atom is: ' +str(Shaik_E_center))
        lst_bonded_atoms = mol_simp_obj.getBondedAtoms(metal_idx)
        
        for bonded_atom_idx in lst_bonded_atoms:
            [bonded_E, bonded_position, bonded_atom, Shaik_E_bonded]  =  Electrostatics.calc_fullE(bonded_atom_idx, all_lines, charge_file, atom_multipole_file)    
            bond_vec_unnorm = np.subtract(np.array(center_position), np.array(bonded_position)) 
            bond_vec = bond_vec_unnorm/(np.linalg.norm(bond_vec_unnorm))
            #initialized a bond_dipole_vec as the (bond_vec_unnorm )*(sum of the partial charges).. can just use dipole! 
            #Compute E-field projected along this bond!
            E_proj = (1/2)*np.dot((np.array(bonded_E) + np.array(center_E)), bond_vec)
            E_proj_Shaik = (1/2)*np.dot((Shaik_E_center + Shaik_E_bonded), bond_vec)
            E_projected.append(E_proj)
            E_shaik_proj.append(E_proj_Shaik)
            bonded_atoms.append(bonded_atom)
            bonded_positions.append(bonded_position)
        print('Projected E!')
        print(E_projected)
        print('Projected Shaik!')
        print(E_shaik_proj)
        print('Bonded ATOMS:')
        print(bonded_atoms)
        print(bonded_positions)
        return [E_projected, bonded_atoms, bonded_positions]
            
    def esp_first_coord(mol_simp_obj, metal_idx, charge_file):
        print('The index of the metal atom is: ' + str(metal_idx))
        lst_bonded_atoms = mol_simp_obj.getBondedAtoms(metal_idx);
        [First_coord_ESP, atom_type] = Electrostatics.calcesp(metal_idx, lst_bonded_atoms, charge_file)
        return First_coord_ESP

    def esp_second_coord(mol_simp_obj, metal_idx, charge_file):
        lst_first_and_second = []
        lst_first_coor_atoms = mol_simp_obj.getBondedAtoms(metal_idx);
        lst_first_and_second.extend(lst_first_coor_atoms)
        for coor_atom_idx in lst_first_coor_atoms:
            second_coor = mol_simp_obj.getBondedAtoms(coor_atom_idx)
            lst_first_and_second.extend(second_coor)
        set_second_coor = set(lst_first_and_second)
        final_lst = list(set_second_coor)
        final_lst.remove(metal_idx)
        [second_coord_ESP, atom_type] = Electrostatics.calcesp(metal_idx, final_lst, charge_file)
        return second_coord_ESP

    #Boolean CageTrue
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

    def ESP_all_calcs(filename, atom_idx, cageTrue):
        #get the number of lines in the txt file
        total_lines =Electrostatics.mapcount(filename)
        all_lines = range(0, total_lines)
        [ESP_all, atom_type] = Electrostatics.calcesp(atom_idx, all_lines, filename)
        if cageTrue:
            cage_lines = range(0, 280)
            guest_lines = range(280, total_lines)
            #print('Cage indices: ' + str(cage_lines))
            #print('Guest indices: ' + str(guest_lines))
            ESP_just_ligand = Electrostatics.calcesp(atom_idx, guest_lines, filename)[0]
            ESP_just_cage = Electrostatics.calcesp(atom_idx, cage_lines, filename)[0]
            print('ESP for all atoms: ' + str(ESP_all) + ' kJ/(mol*e)')
            print('ESP just ligand: ' + str(ESP_just_ligand) + ' kJ/(mol*e)')
            print('ESP jut cage: ' + str(ESP_just_cage) + ' kJ/(mol*e)')
            return [ESP_all, ESP_just_ligand, ESP_just_cage, atom_type]
        else:
            #print('ESP for all atoms: ' + str(ESP_all) + ' kJ/(mol*e)') 
            return [ESP_all, atom_type]


    def esp_bydistance(espatom_idx,  charge_file):
        df = pd.read_csv(charge_file, sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
        k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

        #unit conversion
        A_to_m = 10**(-10)
        KJ_J = 10**-3
        faraday = 23.06   #kcal/(mol*V)
        C_e = 1.6023*(10**-19)
        one_mol = 6.02*(10**23)
        cal_J = 4.184

        #convert each column to list for quicker indexing
        atoms = df['Atom']
        charges = df['charge']
        xs = df['x']
        ys = df['y']
        zs = df['z']

        #pick the index of the atom at which the esp should be calculated
        idx_atom = espatom_idx

        print("The charge at: "+ str(df['Atom'][idx_atom]) +" atom will be calculated")
        #determine position and charge of the target atom
        xo = xs[idx_atom]
        yo = ys[idx_atom]
        zo = zs[idx_atom]
        chargeo = charges[idx_atom]
        total_esp = 0
        #create an ordering of the atoms based on distance from the central atom
        total_atoms = len(xs)
        distances = []
        esps = []
        for idx in range(0, total_atoms):
            if idx == idx_atom:
                continue
            else:
                r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
                distances.append(r)
                esps.append(k*C_e*cal_J*faraday*charges[idx]/r)
        #Now we sort the distance list, and use sorted indices to sort the
        dist_arr = np.array(distances)
        sorted_idx = np.argsort(dist_arr)
        esp_arr = np.array(esps)
        sorted_esps = esp_arr.take(sorted_idx)
        cumulative_esps = np.cumsum(sorted_esps)
        sorted_dist = dist_arr.take(sorted_idx)

        return [sorted_dist, sorted_esps, cumulative_esps]
  
    #Function that can be called on a class object to compute key error analysis metrics for a transition metal complex
    def errorAnalysis(self, csvName):
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders
        owd = os.getcwd() # old working directory
        allspeciesdict = []
        counter = 0
        for f in list_of_file:
            print(f)
            atom_idx = metal_idxs[counter]
            counter = counter + 1
            molsimp_obj = mol3D()
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            results_dir = os.getcwd() + '/'
            try:
                [results_dict, final_mol] = ErrorAnalysis.optim_rmsd_file(results_dir,atom_idx, False,[], self.inGaCageBool)
                molsimp_obj.copymol3D(final_mol)
            except Exception as e:
                results_dict = {}
            results_dict['Name'] = f
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(csvName + '.csv')
        return df

    #list_of_folders = the list of the folders that contain the desired files
    #new_dir: the [post-folder path to the scr folder that contains the .molden and optim.xyz file themselfs
    #dict of calcs, calculations to be performed by multiwavefunction with the corresponding keys
    #newfilanme: desired name of the .csv fiole that will be createcd in getData cotnaining all of the ESP/other data extracted un the file
    ''' Function computes a series of ESP data using the charge scheme specified in charge types
    Accepts:
    charge_types: list of strings
    ESPdata_filename: string
        Name of the output file name
    
    '''
    def getESPData(self, charge_types, ESPdata_filename):
       #Access Class Variables
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # old working directory
        allspeciesdict = []
        counter = 0  #iterator to account for atomic indices of interest
        for f in list_of_file:
            print(f)
            atom_idx = metal_idxs[counter]
            counter = counter + 1
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            subprocess.call("module load multiwfn/noGUI", shell=True)
            command_A = '/opt/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn '+ 'final_optim.molden'
            results_dir = os.getcwd() + '/'
            molsimp_obj= mol3D()
            molsimp_obj.readfromxyz('final_optim.xyz')
            
            results_dict = {}
            results_dict['Name'] = f

            for key in charge_types:
                print('Current key:' + str(key))
                try:
                    full_file_path = os.getcwd() +'/final_optim_' +key+'.txt'
                    if key == "Hirshfeld_I":
                        atmrad_src = "/opt/Multiwfn_3.7_bin_Linux_noGUI/examples/atmrad"
                        copy_tree(atmrad_src, results_dir + 'atmrad/')
                    try: 
                        if self.inGaCageBool:
                            #With newly analyzed partial charges, re-compute ESP data
                            [ESP_all, ESP_just_ligand, ESP_just_cage, atom_type] = Electrostatics.ESP_all_calcs(full_file_path, atom_idx, self.inGaCageBool)
                        else:
                            [ESP_all, atom_type] = Electrostatics.ESP_all_calcs(full_file_path, atom_idx, self.inGaCageBool)

                        [total_charge,partial_charge_atom] = Electrostatics.charge_atom(full_file_path, atom_idx)
                        [sorted_distances, sorted_esps, cum_esps] = Electrostatics.esp_bydistance(atom_idx,  full_file_path)
                        ESP_fcoord = Electrostatics.esp_first_coord(molsimp_obj, atom_idx, full_file_path)
                        ESP_scoord = Electrostatics.esp_second_coord(molsimp_obj, atom_idx, full_file_path)
                    except Exception as e:
                        print('The Exception is: ' + str(e))
                        print(traceback.format_exc())
                        print('Error when trying to access electrostatic information: Attemtping to re-compute partial charges of type: ' + str(key))
                        #Re-run multiwfn computation of partial charge 
                        proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                        calc_command = self.dict_of_calcs[key]
                        commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                        if key == 'CHELPG':
                            commands = ['7', calc_command, '1','\n', 'y', '0', 'q']
                        output = proc.communicate("\n".join(commands).encode())
                        new_name = 'final_optim_' +key+'.txt'
                        os.rename('final_optim.chg', new_name)
                        if self.inGaCageBool:                        
                            #With newly analyzed partial charges, re-compute ESP data
                            [ESP_all, ESP_just_ligand, ESP_just_cage, atom_type] = Electrostatics.ESP_all_calcs(full_file_path, atom_idx, self.inGaCageBool)
                        else:
                            [ESP_all, atom_type] = Electrostatics.ESP_all_calcs(full_file_path, atom_idx, self.inGaCageBool)
                        [total_charge,partial_charge_atom] = Electrostatics.charge_atom(full_file_path, atom_idx)
                        [sorted_distances, sorted_esps, cum_esps] = Electrostatics.esp_bydistance(atom_idx,  full_file_path)
                        ESP_fcoord = Electrostatics.esp_first_coord(molsimp_obj, atom_idx, full_file_path)
                        ESP_scoord = Electrostatics.esp_second_coord(molsimp_obj, atom_idx, full_file_path)

                    #At this point, all calculations shouldbe complete and succesfull: Add ESP data to dictionary
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

                    #if .molden files deal with encapsulated TMCS, complete an additional set of analyses
                    if self.inGaCageBool:
                        [distoGa, Ga_selfdist]=Electrostatics.calcdist(atom_idx, full_file_path)
                        results_dict['MetaltoGa_dist'] = distoGa
                        results_dict['GatoGa_dist'] = Ga_selfdist
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

    def getEFieldData(self, Efield_data_filename):
        metal_idxs = self.lst_of_tmcm_idx
        folder_to_molden = self.folder_to_file_path
        list_of_file = self.lst_of_folders

        owd = os.getcwd() # old working directory
        allspeciesdict = []
        counter = 0

    
        for f in list_of_file:
            atom_idx = metal_idxs[counter]
            counter = counter + 1
            os.chdir(owd)
            os.chdir(f + folder_to_molden)
            subprocess.call("module load multiwfn/noGUI", shell=True)
            #First For this to work, the .molden file should be named: f.molden

            command_A = '/opt/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn '+ 'final_optim.molden'

 

            results_dir = os.getcwd() + '/'
            molsimp_obj = mol3D()
            molsimp_obj.readfromxyz('final_optim.xyz')
            results_dict = {}
            results_dict['Name'] = f
            Command_Polarization = '/opt/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn '+ 'final_optim.molden >' +  'final_optim_polarization.txt'
            #Check if the atomic polarizations have been computed
            path_to_pol = os.getcwd() + '/' + 'final_optim_polarization.txt'
            if os.path.exists(path_to_pol):
                print('I located polarization file for' + f + "!!")
            else:
                print('Starting to run polarization calculation!')
                #Now Run the calculation for atomic dipole and quadrupole moment
                proc = subprocess.Popen(Command_Polarization, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                polarization_commands = ['15', '-1', '3', '2', '0', 'q'] # for atomic dipole and quadrupole moment of Hirshfeld-I type
                output = proc.communicate("\n".join(polarization_commands).encode())
        
            try:
                full_file_path = os.getcwd() +'/' + 'final_optim_Hirshfeld_I.txt'
                atmrad_src = "/opt/Multiwfn_3.7_bin_Linux_noGUI/examples/atmrad"
                copy_tree(atmrad_src, results_dir + 'atmrad/')        
                #This only works if the current key is Hirshfeld I, otherwise unavailable since polarization path should be in hirshfeld-I
                [proj_Efields, bondedAs, bond_pos] = Electrostatics.E_proj_first_coord(molsimp_obj, atom_idx, full_file_path, path_to_pol)
            
            except Exception as e:
                proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
                calc_command = self.dict_of_calcs['Hirshfeld_I']
                commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
                output = proc.communicate("\n".join(commands).encode())
                new_name = 'final_optim'+'_Hirshfeld_I.txt'
                os.rename('final_optim.chg', new_name)

            results_dict['Max Eproj'] = max(abs(np.array(proj_Efields)))
            #probably want to add other bonds to this list!
            allspeciesdict.append(results_dict)
        os.chdir(owd)
        df = pd.DataFrame(allspeciesdict)
        df.to_csv(Efield_data_filename +'.csv')


