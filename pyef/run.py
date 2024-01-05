import os
from DataAnalysis import Electrostatics
from GeometryCheck import ErrorAnalysis


####Steps to use#######
### 1.) Change the paths to the files in list_of_folders
### 2.) Change the index of the atom (zero-indexed) to the metla atom of interest
### 3.) change folder_to_file path such that the path listed in list_of_folders _ folder_to_file_path should provide math to molden
### 4.) If running ESP calculation, pick desired partial charge schemes (see documentation in DataAnalysis.py at top of possible options, un-comment this line and re-set the name of the csv as desired (ESP calc for ~400atom system takes about 15 minutes on one node)
### 5.) If runing E-field Calculation, retitle the name of the outpt file as needed (Efield Calc for ~400atom system takes ~1 hour on one node)
##Load conda environment which contains access to molsimplify
##Run in command line or use Analysis_Jobscript (after adjusting file path, name of conda environment, and file names) to run on gibraltar
##If jobs ends abruptly, will need to delete the polarization file before restarting




list_of_folders = ['/home/manets12/DavidMimichrome/mc6_1_2best_c0opt', '/home/manets12/DavidMimichrome/mc6_2_6best_c0opt']
print(str(list_of_folders))

#For atoms in molecules corresponding to the folder, identify the index of atom at which to compute ESP
#Note that this uses zero-indexing!
list_of_atoms = len(list_of_folders)*[0]
#adjust for differing indices of the metal number
list_of_atoms[0] = 486
list_of_atoms[1] = 486

#path from each initial directory to the directory containing the .molden file
folder_to_file_path  = '/scr/'

#Initialize Electrostatics Object
incage_bool = False
dataObject = Electrostatics(list_of_folders, list_of_atoms, folder_to_file_path, incage_bool)
#fix/reformat the name and charges on atoms in .molden file to set up for future calculations
dataObject.prepData()
dataObject.fix_ECPmolden()

#Define name of .csv to contain error data
err_csv_name = 'gfn2xtb_Error'

#Method will calculate various RMSD, molsimplify error parameters/flags
#dataObject.errorAnalysis(err_csv_name)

#Determine Filename prefix for output 
ESPdata_filename = 'MinimiChrome_ESP'

#List of partial C... for ESP also need a list of partial charges of interest
lst_charge_types = ['Hirshfeld_I']

#Create CSV with ESP data
#dataObject.getESPData(lst_charge_types, ESPdata_filename)


#Define prefix name  of E-field datafile
Efield_data_filename = 'MinimiChrome_Efield'

#Method to Compute Efield Projections on bonds connected to the atom specified by index in list_of_atoms
dataObject.getEFieldData(Efield_data_filename)

