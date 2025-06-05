import os
import logging
import numpy as np
import pandas as pd
from collections import deque
from molSimplify.Scripts import *
from molSimplify.Classes.mol3D import *
from molSimplify.Classes.ligand import ligand_breakdown
from molSimplify.Classes.globalvars import all_angle_refs













class ErrorAnalysis:
    def __init__(self, struct_type):
        self.struct_type = struct_type

    def goodstruct(dict_results):
        """Function will return one if all conditions are met for valid geometry, otherwise will return 0"""
        oct_angle_devi_max = dict_results['oct_angle_devi_max']
        numcoord = dict_results['num_coord_metal']
        max_del_sig_angle = dict_results['max_del_sig_angle']
        dist_del_all = dict_results['dist_del_all']
        dist_del_all_relative = dict_results['dist_del_all_relative']

        if numcoord != 4:
           return 0
        elif oct_angle_devi_max > 22:
            return 0
        elif max_del_sig_angle > 22:
            return 0
        elif dist_del_all > 1:
            return 0
        elif dist_del_all_relative > 0.30:
            return 0
        else:
            return 1
  
    def rmsd_noHs(init_mol, opt_mol):
        new_init = mol3D()
        new_init.copymol3D(init_mol)
        new_opt = mol3D()
        new_opt.copymol3D(opt_mol)
        new_init.deleteHs()
        new_opt.deleteHs()
        opt_rmsd_noHs = rmsd.rigorous_rmsd(new_opt, new_init)
        return opt_rmsd_noHs

    def rmsd_Hs(init_mol, opt_mol):
        opt_rmsd = rmsd.rigorous_rmsd(init_mol, opt_mol)
        return opt_rmsd

    def molsimp_struct(init_mol3D, mol_test,skip_tests):
        # Make a copt of the mol3D object so you don't over-write previous mol3D
        new_test = mol3D()
        new_test.copymol3D(mol_test)
        flag_results = new_test.IsStructure(init_mol=init_mol3D, num_coord=4, angle_ref=all_angle_refs["square planar"], skip=skip_tests)
        dict_results = flag_results[2]
        # Now check if the structure is good using in-script definitions
        boolgoodstruct = ErrorAnalysis.goodstruct(dict_results)
        return [boolgoodstruct, dict_results]
      

    # Determine the smiles string associated with ligands in molecule
    def ligand_info(final_mol):
        # First check if the ligand is inside of cage..
        all_ligand_smiles = []
        total_atoms = len(final_mol.getAtoms())

        if total_atoms > 280: 
            list_guest_atoms = range(280, total_atoms)
            guest_mol = final_mol.create_mol_with_inds(list_guest_atoms)
        else:
            guest_mol = final_mol
        ligList, ligDent, ligConIndices = ligand_breakdown(guest_mol)
        for lig in ligList:
            lig_mol = guest_mol.create_mol_with_inds(lig)
            smi = lig_mol.get_smiles()
            all_ligand_smiles.append(smi)
        return [all_ligand_smiles, ligDent]

    def cage_rmsd(final_mol):
        lst_cage_mols = range(0, 280)
        cage_mol = mol3D()
        cage_mol = final_mol.create_mol_with_inds(lst_cage_mols)

        #Create a mol associated with the actual cage structure
        init_cage = mol3D()
        init_cage.readfromxyz('/home/manets12/Nanocages/docking_unfrozencage/empty_cage.xyz')
        #print('Comaprison of cage coordinates:')
        #print('Coordinate of initial cage')
        #print(init_cage.coords())
        #print('Coordinates of the final cage')
        #print(cage_mol.coords())

        noH_rmsd = ErrorAnalysis.rmsd_noHs(cage_mol, init_cage)
        H_rmsd = ErrorAnalysis.rmsd_Hs(cage_mol, init_cage)
        print('H RMSD cage: ' + str(H_rmsd))
        print('no H RMSD cage: ' + str(noH_rmsd))
        return [H_rmsd, noH_rmsd]
    
    def rmsd_TMC_cage(final_mol, path_to_init):
        total_atoms = len(final_mol.getAtoms())
        list_guest_atoms = range(280, total_atoms)
        guest_mol = final_mol.create_mol_with_inds(list_guest_atoms)
        dft_opt_mol = mol3D()
        dft_opt_mol.readfromxyz(path_to_init + 'optimized_tmc.xyz')
        noH_rmsd = ErrorAnalysis.rmsd_noHs(dft_opt_mol, guest_mol)
        H_rmsd = ErrorAnalysis.rmsd_Hs(dft_opt_mol, guest_mol)
        return [H_rmsd, noH_rmsd]

    #same functionality as optim_rmsd_listfiles but only for a single file
    def optim_rmsd_file(optim_path, atom_idx, isCSD, compare_paths=[], isinCage=False):
        ''' Takes in several paths to folder that contain output, optim trajectories and then performs a series of geometry
    checks as specified in checks list. isCSD is a boolean that specifies if the file should be checked against the initial frame of the
    optim.xyz file

    Parameters
    ----------
    optim_paths: list of strings
        contains a list of strings for paths into directory that contain optim.xyz (based on Terachem output formalism)
    isCSD: Boolean
        determine if the original xyz file is from CSD (in this case, the first xyz and last xyz frame are compared to provide rmsd)
    compare_paths: list of strings, default is empty list
       if initial file is not from csd and desirable to compare the last frame of the optim.xyz file with another xyz file (i.e. isCSD = False), input a list of the paths to the other xyz files here
    isinCage: Boolean
        determine if the xyz file is for a TMC inside a cage (if True, is in cage)

    Returns
    -------
    returndict = dictionary containing relevant info for square planar structs
        contains the list of checks and values implemented in the code

        '''
        init_dir = os.getcwd()
        returndict = {}
        try:
            path_to_init = optim_path[:-4]
            os.chdir(optim_path)
            optim_file = 'optim.xyz'
            full_traj = open(optim_file, 'r')
            num_atoms = int(full_traj.readline())
            num_lines = num_atoms + 2

            with open(optim_file) as input_file:
                head = [next(input_file) for _ in range(num_lines)]
            with open('initial_' + optim_file, 'w') as initxyz:
                initxyz.writelines(head)
            #last xyz in trajectory saved as the final file
            with open('final_'+optim_file, 'w') as finalxyz:
                finalxyz.writelines(deque(full_traj, num_lines))
            
            init_mol = mol3D()
            init_mol.readfromxyz('initial_' + optim_file)
            final_mol = mol3D()
            final_mol.readfromxyz('final_' + optim_file)

            #Use the initial mole file if only one step in optimization (i.e. final is same as initial file)
            with open('final_' + optim_file, 'r') as read_final:
                first_line = read_final.readline()
                print(first_line[0])
                if first_line[0] == '-':
                    final_mol = mol3D()
                    final_mol.copymol3D(init_mol)
            #will over-write the final xyz file with the initial file to avoid erroring
                    
            
            #list of geometry checks to skip when implementing IsStructure
            skipGeomChecks = ['atom_dist_max', 'banned_by_user', 'dist_del_eq', 'devi_linear_avrg', 'devi_linear_max', 'dist_del_eq_relative']
            returndict = {}
            try:
                [boolgoodstruct, dict_errs] = ErrorAnalysis.molsimp_struct(init_mol, final_mol, skipGeomChecks)
            
                for skipcheck in skipGeomChecks:
                    dict_errs.pop(skipcheck, None)
                dict_errs['IsGood?'] = str(boolgoodstruct)
                returndict = dict_errs
            except Exception as e:
                print(e)
            dictMLBL = final_mol.getMLBondLengths()[atom_idx]
            #Create a new molecule object with only the inner guest
            [all_ligand_smiles, ligDent] = ErrorAnalysis.ligand_info(final_mol)
            if isinCage:
                 [cage_rmsd_H, cage_rmsd_noHs] = ErrorAnalysis.cage_rmsd(final_mol)
                 returndict['cage_RMSD_withH'] = cage_rmsd_H 
                 returndict['cage_RMSD_noH'] = cage_rmsd_noHs
                 [tmc_rmsd_H, tmc_rmsd_noHs] = ErrorAnalysis.rmsd_TMC_cage(final_mol, path_to_init)
                 returndict['tmc_RMSD_withH'] = tmc_rmsd_H
                 returndict['tmc_RMSD_noHs'] = tmc_rmsd_noHs
            else:
                 pass
            returndict['LigandSmiles'] = all_ligand_smiles;
            returndict['LigandDenticity'] = ligDent;
            returndict['Avg_M-Lbondlengths'] = np.average(dictMLBL['M-L bond lengths'])
            returndict['Avg_M-Lbonddev'] = np.average(dictMLBL['relative bond lengths'])
            returndict['RMSD_noH'] = ErrorAnalysis.rmsd_noHs(init_mol, final_mol)
            returndict['RMSD_withH'] = ErrorAnalysis.rmsd_Hs(init_mol, final_mol)
        except Exception as e:
            logging.exception('An Exception was thrown')

        os.chdir(init_dir)
        return [returndict, final_mol]

            #This function executes the desired geometry checks for a list of stru

    def optim_rmsd_listfiles(self, optim_paths, isCSD, checks, compare_paths=[]):
        
        ''' Takes in several paths to folder that contain output, optim trajectories and then performs a series of geometry 
    checks as specified in checks list. isCSD is a boolean that specifies if the file should be checked against the initial frame of the 
    optim.xyz file

    Parameters
    ----------
        optim_paths: list of strings
        contains a list of strings for paths into directory that contain scr/optim.xyz (based on Terachem output formalism)
        isCSD: Boolean
            determine if the original xyz file is from CSD (in this case, the first xyz and last xyz frame are compared to provide rmsd)
        checks: list of strings
            list of strings with possible geometry checks. Must be one of the following: 'noHs_rmsd', 'rmsd', 'BL_histo', 'BL'
        compare_paths: list of strings, default is empty list
            if initial file is not from csd and desirable to compare the last frame of the optim.xyz file with another xyz file (i.e. isCSD = False), input a list of the paths to the other xyz files here 
    Returns
    -------
      None
        '''

        # First we will just design this for a single optim.xyz file and check that hte ifal and initial one cna be compared
        # Write initial xyz and final xyz into two separate files:
        # Find optim file and get the output of it as the final xyz

        owd = os.getcwd()
        print('Current directory: '+owd)
    
        # Create lists to store the structure parameters to be collected for each output file

        list_error_dicts = []
        list_structure_names = []
        list_rmsd = []
        list_rmsd_noHs = []
        for optim_path in optim_paths:
            try:
                os.chdir(optim_path+'/scr/') 
                in_dir = os.getcwd()
                print('Changed directory to: '+ in_dir)     
                optim_file = 'optim.xyz'
                full_traj =  open(optim_file, 'r')
                num_atoms = int(full_traj.readline())
                num_lines = num_atoms + 2
        
                # Get the first xyz file
                with open(optim_file) as input_file:
                    head = [next(input_file) for _ in range(num_lines)]

                # Last xyz in trajectory saved as the final file
                with open('final_'+optim_file, 'w') as finalxyz:
                    finalxyz.writelines(deque(full_traj, num_lines))

                # Initial xyz in trajectory saved as the initial file
                with open('initial_'+optim_file, 'w') as initxyz:
                    initxyz.writelines(head)
        
                init_mol = mol3D()
                init_mol.readfromxyz('final_' + optim_file)
                final_mol = mol3D()
                final_mol.readfromxyz('initial_'+optim_file)
        
                # List of geometry cehcks to skip when implementing IsStructure
                skipGeomChecks = ['atom_dist_max', 'banned_by_user', 'dist_del_eq', 'devi_linear_avrg', 'devi_linear_max', 'dist_del_eq_relative']
                list_rmsd_noHs.append(rmsd_noHs(init_mol, final_mol))
                list_rmsd.append(rmsd_Hs(init_mol, final_mol))
                [boolgoodstruct, dict_errs] = molsimp_struct(init_mol, final_mol, skipGeomChecks)
                for skipcheck in skipGeomChecks:
                    dict_errs.pop(skipcheck, None)
                dict_errs['IsGood?'] = str(boolgoodstruct)
                list_error_dicts.append(dict_errs)
                list_structure_names.append(optim_path)
                os.chdir(owd)
            except Exception as e:
                # Print("An Exception Occured when trying to procees: " + str(optim_file)) 
                logging.exception('An Exception was thrown')
        df = pd.DataFrame(list_error_dicts)
        df['Molecule'] = list_structure_names
        df = df.set_index('Molecule')
        df['RMSD'] = list_rmsd
        df['RMSD no Hs'] = list_rmsd_noHs
        df.to_csv('error_structs.csv')



