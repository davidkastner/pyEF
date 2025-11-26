import os
import sys
import pyef
import argparse
from pyef.analysis import Electrostatics

def main(job_name, jobs, metal_indices, bond_indices, dielectric, geom_flag, esp_flag, multiwfn_module, multiwfn_path, atmrad_path):
    """
    Main function for running the pyEF workflow.

    Notes
    -----
    1.) Change the paths to the files in jobs
    2.) Change the index of the atom (zero-indexed) to the metla atom of interest
    3.) change folder_to_file path such that the path listed in jobs _ folder_to_file_path should provide math to molden
    4.) If running ESP calculation, pick desired partial charge schemes (see documentation in DataAnalysis.py at top of possible options, un-comment this line and re-set the name of the csv as desired (ESP calc for ~400atom system takes about 15 minutes on one node)
    5.) If runing E-field Calculation, retitle the name of the outpt file as needed (Efield Calc for ~400atom system takes ~1 hour on one node)
    Load conda environment which contains access to molsimplify

    Run in command line or use Analysis_Jobscript (after adjusting file path, name of conda environment, and file names) to run on gibraltar
    If jobs ends abruptly, will need to delete the polarization file before restarting

    """

    # Path from each initial directory to the directory containing the .molden file
    folder_to_file_path  = '/scr/'

    # Initialize Electrostatics Object
    incage_bool = False
    dataObject = Electrostatics(jobs, metal_indices, folder_to_file_path, incage_bool=incage_bool)

    # Fix/reformat the name and charges on atoms in .molden file to set up for future calculations
    dataObject.prepData()
    dataObject.fix_allECPmolden()

    # Define name of .csv to contain error data
    err_csv_name = 'Errordata'

    if geom_flag:
        # Method will calculate various RMSD, molsimplify error parameters/flags
        # NOTE: errorAnalysis() method does not exist in Electrostatics class
        # TODO: Implement errorAnalysis() or remove this functionality
        # dataObject.errorAnalysis(err_csv_name)
        pass

    # Determine Filename prefix for output 
    ESPdata_filename = 'ESPdata'

    # List of partial C... for ESP also need a list of partial charges of interest
    lst_charge_types = ['Hirshfeld_I']

    if esp_flag:
        # Create CSV with ESP data
        dataObject.getESPMultipole(lst_charge_types, ESPdata_filename, multiwfn_module, multiwfn_path, atmrad_path, dielectric)

    # Method to Compute Efield Projections on bonds connected to the atom specified by index in metal_indices
    dataObject.getEFieldMultipole(job_name, multiwfn_module, multiwfn_path, atmrad_path, input_bond_indices=bond_indices)

def read_file_lines(file_path):
    """Reads in auxiliary files containing job information"""
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

if __name__ == "__main__":
    # Example: python run.py --geom --esp --jobs_file path/to/jobs.in --metals_file path/to/metals.in --job_name myJob --bonds_file path/to/bonds.in --multiwfn_module multiwfn --multiwfn_path /path/to/multiwfn --atmrad_path /path/to/atmrad --dielectric 1.0 > pyEF.log
    parser = argparse.ArgumentParser(description="Script Description")
    parser.add_argument("--geom", action="store_true", help="Perform a geometry check")
    parser.add_argument("--esp", action="store_true", help="Perform analysis of electrostatics")
    parser.add_argument("--jobs_file", required=True, help="Path to file containing job paths")
    parser.add_argument("--metals_file", required=True, help="Path to file containing metal indices")
    parser.add_argument("--job_name", required=True, help="Name for the job/output files")
    parser.add_argument("--bonds_file", required=True, help="Path to file containing bond indices")
    parser.add_argument("--multiwfn_module", required=True, help="Module name for multiwfn")
    parser.add_argument("--multiwfn_path", required=True, help="Path to multiwfn executable")
    parser.add_argument("--atmrad_path", required=True, help="Path to atmrad file")
    parser.add_argument("--dielectric", type=float, default=1.0, help="Dielectric constant (default: 1.0)")

    args = parser.parse_args()
    geom_flag = args.geom
    esp_flag = args.esp
    jobs = read_file_lines(args.jobs_file)
    metal_indices = [int(idx) for idx in read_file_lines(args.metals_file)]
    bond_indices = [int(idx) for idx in read_file_lines(args.bonds_file)]

    main(args.job_name, jobs, metal_indices, bond_indices, args.dielectric,
         geom_flag, esp_flag, args.multiwfn_module, args.multiwfn_path, args.atmrad_path)
