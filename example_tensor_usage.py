"""
Simple example showing how to use the new tensor-based multipole expansion.

This is a minimal working example - modify the parameters for your system.

Author: Generated for pyEF multipole expansion feature
Date: 2025-12-02
"""

from pyef.analysis import Electrostatics
import numpy as np


def simple_example():
    """
    Simple example of using the new get_Electrostatic_stabilization_tensor() function.
    """

    # ===========================
    # STEP 1: Setup
    # ===========================

    # Define your calculation folders
    list_of_folders = ['AXEDEN']  # MODIFY THIS for your folders

    # One substrate per folder (modify as needed)
    list_of_atoms = [0 for _ in range(len(list_of_folders))]

    # Path to calculation files
    folder_to_file_path = '/scr/'

    # Create Electrostatics object
    estat = Electrostatics(
        list_of_folders,
        list_of_atoms,
        folder_to_file_path,
        molden_filename='input.molden',
        xyzfilename='xyz.xyz',
        hasECP=True
    )

    # ===========================
    # STEP 2: Define Multiwfn paths
    # ===========================

    multiwfn_module = "module load multiwfn/noGUI_3.7"
    multiwfn_path = '/data1/groups/HJKgroup/src/multiwfn/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn'
    atmrad_path = "/data1/groups/HJKgroup/src/multiwfn/Multiwfn_3.7_bin_Linux_noGUI/examples/atmrad"

    # ===========================
    # STEP 3: Define substrate and environment
    # ===========================

    # Example: substrate is atoms 0-10
    substrate_indices = [np.arange(0, 11)]  # MODIFY THIS for your substrate!

    # Environment is automatically all other atoms (or you can specify env_idxs)

    # ===========================
    # STEP 4: Run the calculation
    # ===========================

    print("Running tensor-based multipole expansion...")
    print(f"Substrate atoms: {substrate_indices[0]}")
    print(f"Multipole order: 2 (includes dipole-dipole)")
    print()

    # Run with order 2 (includes monopole, dipole, and all interactions)
    df = estat.get_Electrostatic_stabilization_tensor(
        multiwfn_path=multiwfn_path,
        multiwfn_module=multiwfn_module,
        atmrad_path=atmrad_path,
        substrate_idxs=substrate_indices,
        charge_type='Hirshfeld_I',
        name_dataStorage='output_tensor',
        multipole_order=2,  # 1=monopole, 2=+dipole, 3=+quadrupole
        decompose_atomwise=True,  # Get per-atom breakdown
        visualize=True  # Create PDB visualization
    )

    # ===========================
    # STEP 5: Process results
    # ===========================

    if isinstance(df, tuple):
        df_total, df_atomwise = df
        print("Results:")
        print(f"  Total energy: {df_total['Total_Energy_kcal_mol'].iloc[0]:.4f} kcal/mol")
        print(f"  Number of substrate atoms: {df_total['Num_Substrate_Atoms'].iloc[0]}")
        print(f"  Number of environment atoms: {df_total['Num_Environment_Atoms'].iloc[0]}")
        print()

        # Show top contributing atoms
        df_sorted = df_atomwise.sort_values('Energy_Contribution_kcal_mol', ascending=False)
        print("Top 5 contributing substrate atoms:")
        for idx, row in df_sorted.head(5).iterrows():
            print(f"  Atom {int(row['Atom_Index']):3d} ({row['Atom_Symbol']:2s}): "
                  f"{row['Energy_Contribution_kcal_mol']:8.4f} kcal/mol")
        print()

        print("Output files created:")
        print("  • output_tensor.csv (total energies)")
        print("  • output_tensor_atomwise.csv (per-atom breakdown)")
        print("  • PDB files with contributions in B-factor column")

    else:
        print(f"Total energy: {df['Total_Energy_kcal_mol'].iloc[0]:.4f} kcal/mol")
        print("Output file: output_tensor.csv")


def advanced_example():
    """
    Advanced example showing different multipole orders.
    """

    # Setup (same as simple example)
    list_of_folders = ['AXEDEN']
    list_of_atoms = [0]
    folder_to_file_path = '/scr/'

    estat = Electrostatics(
        list_of_folders, list_of_atoms, folder_to_file_path,
        molden_filename='input.molden', xyzfilename='xyz.xyz', hasECP=True
    )

    multiwfn_module = "module load multiwfn/noGUI_3.7"
    multiwfn_path = '/data1/groups/HJKgroup/src/multiwfn/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn'
    atmrad_path = "/data1/groups/HJKgroup/src/multiwfn/Multiwfn_3.7_bin_Linux_noGUI/examples/atmrad"

    substrate_indices = [np.arange(0, 11)]

    print("=" * 80)
    print("Comparing different multipole orders:")
    print("=" * 80)

    # ===========================
    # Order 1: Monopole only
    # ===========================
    print("\n[1] Monopole only (q×q):")
    df1 = estat.get_Electrostatic_stabilization_tensor(
        multiwfn_path, multiwfn_module, atmrad_path,
        substrate_indices, charge_type='Hirshfeld_I',
        name_dataStorage='output_order1',
        multipole_order=1
    )
    E1 = df1['Total_Energy_kcal_mol'].iloc[0]
    print(f"    Energy = {E1:.4f} kcal/mol")

    # ===========================
    # Order 2: + Dipole terms
    # ===========================
    print("\n[2] Order 2 (q×q + q×μ + μ×μ):")
    df2 = estat.get_Electrostatic_stabilization_tensor(
        multiwfn_path, multiwfn_module, atmrad_path,
        substrate_indices, charge_type='Hirshfeld_I',
        name_dataStorage='output_order2',
        multipole_order=2
    )
    E2 = df2['Total_Energy_kcal_mol'].iloc[0]
    print(f"    Energy = {E2:.4f} kcal/mol")
    print(f"    Dipole contribution = {E2 - E1:.4f} kcal/mol")

    # ===========================
    # Order 3: + Quadrupole terms
    # ===========================
    print("\n[3] Order 3 (all above + quadrupole terms):")
    df3 = estat.get_Electrostatic_stabilization_tensor(
        multiwfn_path, multiwfn_module, atmrad_path,
        substrate_indices, charge_type='Hirshfeld_I',
        name_dataStorage='output_order3',
        multipole_order=3
    )
    E3 = df3['Total_Energy_kcal_mol'].iloc[0]
    print(f"    Energy = {E3:.4f} kcal/mol")
    print(f"    Quadrupole contribution = {E3 - E2:.4f} kcal/mol")

    print("\n" + "=" * 80)
    print("Summary:")
    print("=" * 80)
    print(f"  Monopole:             {E1:10.4f} kcal/mol")
    print(f"  + Dipole terms:       {E2:10.4f} kcal/mol  (Δ = {E2-E1:+8.4f})")
    print(f"  + Quadrupole terms:   {E3:10.4f} kcal/mol  (Δ = {E3-E2:+8.4f})")
    print()


if __name__ == "__main__":
    print("=" * 80)
    print("TENSOR-BASED MULTIPOLE EXPANSION - USAGE EXAMPLES")
    print("=" * 80)
    print()

    print("Running simple example...")
    print("-" * 80)
    try:
        simple_example()
    except Exception as e:
        print(f"ERROR in simple example: {e}")
        print("Make sure to modify the script for your specific system!")

    print("\n" + "=" * 80)
    print("To run the advanced example, uncomment the line below")
    print("=" * 80)
    # Uncomment to run advanced example:
    # advanced_example()
