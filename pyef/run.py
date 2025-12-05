import os
import sys
import pyef
import argparse
from pyef.analysis import Electrostatics

def main(job_name, jobs, metal_indices, bond_indices, dielectric,
         geom_flag, esp_flag, ef_flag, estab_flag,
         multiwfn_module, multiwfn_path, atmrad_path,
         charge_types=['Hirshfeld_I'], multipole_bool=True, use_multipole=False,
         decompose_atomwise=False, substrate_idxs=None, env_idxs=None, multipole_order=2,
         include_ptchgs=False, ptchg_file='', dielectric_scale=1.0):
    """
    Main function for running the pyEF workflow.

    Parameters
    ----------
    job_name : str
        Name for output files
    jobs : list of str
        List of job paths
    metal_indices : list of int or None
        List of metal atom indices (only required for ESP analysis or auto-finding bonds)
    bond_indices : list of tuples
        List of bond pairs (required for E-field if metal_indices not provided)
    dielectric : float
        Dielectric constant
    geom_flag : bool
        Run geometry checks
    esp_flag : bool
        Run ESP analysis
    ef_flag : bool
        Run electric field analysis
    estab_flag : bool
        Run electrostatic stabilization analysis
    multiwfn_module : str
        Module name for Multiwfn
    multiwfn_path : str
        Path to Multiwfn executable
    atmrad_path : str
        Path to atmrad file
    charge_types : list of str, optional
        Charge partitioning schemes (default: ['Hirshfeld_I'])
    multipole_bool : bool, optional
        Use multipole expansion for E-field (default: True)
    use_multipole : bool, optional
        Use multipole expansion for ESP (default: False)
    decompose_atomwise : bool, optional
        Compute atom-wise decomposition (default: False)
    substrate_idxs : list or None, optional
        Substrate atom indices for estab calculation (default: None)
    env_idxs : list or None, optional
        Environment atom indices for estab calculation (default: None)
    multipole_order : int, optional
        Multipole expansion order (default: 2)
    include_ptchgs : bool, optional
        Include MM point charges from QM/MM calculations (default: False)
    ptchg_file : str, optional
        Path to point charge file (default: '')
    dielectric_scale : float, optional
        Dielectric scaling factor for MM charges (default: 1.0)

    Notes
    -----
    - ESP calculation for ~400 atom system takes about 15 minutes on one node
    - E-field calculation for ~400 atom system takes ~1 hour on one node
    - If jobs end abruptly, delete the polarization file before restarting
    """

    # Path from each initial directory to the directory containing the .molden file
    # NOTE: If jobs already contain complete paths, set folder_to_file_path = ''
    folder_to_file_path  = '/scr/'

    # Initialize Electrostatics Object
    incage_bool = False
    # Metal indices only needed for ESP calculations or auto-finding bonds
    if metal_indices:
        dataObject = Electrostatics(jobs, folder_to_file_path, lst_of_tmcm_idx=metal_indices,
                                    incage_bool=incage_bool, dielectric=dielectric)
    else:
        dataObject = Electrostatics(jobs, folder_to_file_path,
                                    incage_bool=incage_bool, dielectric=dielectric)

    # Configure QM/MM if requested
    if include_ptchgs and ptchg_file:
        dataObject.includePtChgs(ptchg_file)
        dataObject.set_dielectric_scale(dielectric_scale)

    # Fix/reformat the name and charges on atoms in .molden file to set up for future calculations
    dataObject.prepData()
    dataObject.fix_allECPmolden()

    if geom_flag:
        # Method will calculate various RMSD, molsimplify error parameters/flags
        # NOTE: errorAnalysis() method does not exist in Electrostatics class
        # TODO: Implement errorAnalysis() or remove this functionality
        # dataObject.errorAnalysis('Errordata')
        pass

    # Run ESP analysis if requested
    if esp_flag:
        print("Performing electrostatic potential (ESP) analysis...")
        ESPdata_filename = f'{job_name}_ESPdata'
        dataObject.getESP(charge_types, ESPdata_filename, multiwfn_module,
                         multiwfn_path, atmrad_path, use_multipole=use_multipole,
                         dielectric=dielectric)

    # Run E-field analysis if requested
    if ef_flag:
        print("Performing electric field analysis...")
        Efielddata_filename = f'{job_name}_Efielddata'
        dataObject.getEfield(charge_types[0] if isinstance(charge_types, list) else charge_types,
                            Efielddata_filename, multiwfn_module, multiwfn_path, atmrad_path,
                            multipole_bool=multipole_bool, input_bond_indices=bond_indices,
                            decompose_atomwise=decompose_atomwise, dielectric=dielectric)

    # Run electrostatic stabilization analysis if requested
    if estab_flag:
        if substrate_idxs is None:
            print("Warning: substrate_idxs not specified for electrostatic stabilization. Skipping...")
        else:
            print("Performing electrostatic stabilization analysis...")
            estab_filename = f'{job_name}_Estab'
            dataObject.getElectrostatic_stabilization(
                multiwfn_path, multiwfn_module, atmrad_path,
                substrate_idxs=substrate_idxs,
                charge_type=charge_types[0] if isinstance(charge_types, list) else charge_types,
                name_dataStorage=estab_filename,
                env_idxs=env_idxs,
                decompose_atomwise=decompose_atomwise,
                multipole_order=multipole_order
            )

def read_file_lines(file_path):
    """Reads in auxiliary files containing job information"""
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

if __name__ == "__main__":
    # Example: python run.py --ef --esp --estab --jobs_file path/to/jobs.in --metals_file path/to/metals.in --job_name myJob --bonds_file path/to/bonds.in --multiwfn_module multiwfn --multiwfn_path /path/to/multiwfn --atmrad_path /path/to/atmrad --dielectric 1.0 --substrate_idxs 1,2,3 > pyEF.log
    parser = argparse.ArgumentParser(description="PyEF Analysis Pipeline")
    parser.add_argument("--geom", action="store_true", help="Perform a geometry check")
    parser.add_argument("--esp", action="store_true", help="Perform electrostatic potential (ESP) analysis")
    parser.add_argument("--ef", action="store_true", help="Perform electric field analysis")
    parser.add_argument("--estab", action="store_true", help="Perform electrostatic stabilization analysis")
    parser.add_argument("--jobs_file", required=True, help="Path to file containing job paths")
    parser.add_argument("--metals_file", required=True, help="Path to file containing metal indices")
    parser.add_argument("--job_name", required=True, help="Name for the job/output files")
    parser.add_argument("--bonds_file", required=True, help="Path to file containing bond indices")
    parser.add_argument("--multiwfn_module", required=True, help="Module name for multiwfn")
    parser.add_argument("--multiwfn_path", required=True, help="Path to multiwfn executable")
    parser.add_argument("--atmrad_path", required=True, help="Path to atmrad file")
    parser.add_argument("--dielectric", type=float, default=1.0, help="Dielectric constant (default: 1.0)")
    parser.add_argument("--charge_types", nargs='+', default=['Hirshfeld_I'], help="Charge partitioning schemes (default: Hirshfeld_I)")
    parser.add_argument("--multipole", action="store_true", default=True, help="Use multipole expansion for E-field")
    parser.add_argument("--use_multipole", action="store_true", help="Use multipole expansion for ESP")
    parser.add_argument("--decompose_atomwise", action="store_true", help="Compute atom-wise decomposition")
    parser.add_argument("--substrate_idxs", type=str, help="Comma-separated substrate atom indices for estab (e.g., 1,2,3)")
    parser.add_argument("--env_idxs", type=str, help="Comma-separated environment atom indices for estab")
    parser.add_argument("--multipole_order", type=int, default=2, help="Multipole expansion order (default: 2)")

    args = parser.parse_args()
    geom_flag = args.geom
    esp_flag = args.esp
    ef_flag = args.ef
    estab_flag = args.estab
    jobs = read_file_lines(args.jobs_file)
    metal_indices = [int(idx) for idx in read_file_lines(args.metals_file)]
    bond_indices = [int(idx) for idx in read_file_lines(args.bonds_file)]

    # Parse substrate and environment indices if provided
    substrate_idxs = [int(x) for x in args.substrate_idxs.split(',')] if args.substrate_idxs else None
    env_idxs = [int(x) for x in args.env_idxs.split(',')] if args.env_idxs else None

    main(args.job_name, jobs, metal_indices, bond_indices, args.dielectric,
         geom_flag, esp_flag, ef_flag, estab_flag,
         args.multiwfn_module, args.multiwfn_path, args.atmrad_path,
         args.charge_types, args.multipole, args.use_multipole,
         args.decompose_atomwise, substrate_idxs, env_idxs, args.multipole_order)
