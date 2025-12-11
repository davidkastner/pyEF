import os
from analysis import Electrostatics
#from geometry import ErrorAnalysis


####Steps to use#######
### 1.) Provide complete paths to your calculation folders in list_of_folders
### 2.) For ESP calculations: Specify metal atom indices (zero-indexed) in lst_of_tmcm_idx
### 3.) For E-field calculations: Specify bond indices in input_bond_indices
### 4.) For stabilization calculations: Specify substrate_idxs and env_idxs
### 5.) Update paths to multiwfn and atmrad executables
##Load conda environment which contains access to required tools
##Run in command line or submit as a batch job
##If job ends abruptly, you may need to delete polarization files before restarting


# NEW RECOMMENDED USAGE: Provide complete folder paths
list_of_folders = [
    '/home/manets12/DavidMimichrome/ABAJOD/scr',
    '/home/manets12/DavidMimichrome/CAXFOX/scr'
]
print(str(list_of_folders))

# For ESP calculations only: specify metal atom indices (zero-indexed)
metal_indices = [5, 7]  # Metal atom indices for each folder

# Initialize Electrostatics Object (new simplified signature)
dataObject = Electrostatics(
    list_of_folders,
    lst_of_tmcm_idx=metal_indices,  # Only needed for ESP calculations
    dielectric=4.0  # Optional: set dielectric constant
)

# Prepare data: fix/reformat .molden files if needed
dataObject.prepData()
dataObject.fix_allECPmolden()  # Use if ECP basis sets were used

# ========================================
# EXAMPLE 1: Calculate ESP at metal centers
# ========================================
esp_df = dataObject.getESP(
    charge_types=['Hirshfeld_I'],
    ESPdata_filename='MinimiChrome_ESP',
    multiwfn_path='/path/to/multiwfn',
    atmrad_path='/path/to/atmrad',
    use_multipole=True
)
print("ESP calculation complete!")

# ========================================
# EXAMPLE 2: Calculate E-field on specific bonds
# ========================================
bond_indices = [(1, 2), (5, 6)]  # Bonds to analyze
efield_df = dataObject.getEfield(
    charge_types='Hirshfeld_I',
    Efielddata_filename='MinimiChrome_Efield',
    multiwfn_path='/path/to/multiwfn',
    atmrad_path='/path/to/atmrad',
    input_bond_indices=bond_indices,
    multipole_bool=True
)
print("E-field calculation complete!")

# ========================================
# EXAMPLE 3: Calculate Electrostatic Stabilization
# ========================================
estab_df = dataObject.getElectrostatic_stabilization(
    multiwfn_path='/path/to/multiwfn',
    atmrad_path='/path/to/atmrad',
    substrate_idxs=[1, 2, 3, 4, 5],  # Substrate atoms
    env_idxs=[6, 7, 8, 9, 10],       # Environment atoms
    charge_type='Hirshfeld_I',
    multipole_order=2
)
print("Electrostatic stabilization calculation complete!")
