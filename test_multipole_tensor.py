"""
Test script comparing multipole expansion methods for electrostatic stabilization.

This script demonstrates the differences between:
1. Field-based approach (old): substrate responds to environment-generated fields
2. Tensor-based approach (new): direct multipole-multipole interactions

Author: Generated for pyEF multipole expansion feature
Date: 2025-12-02
"""

from pyef.analysis import Electrostatics
import os
import numpy as np
import pandas as pd


def compare_multipole_methods():
    """
    Compare different multipole expansion approaches.

    Computes electrostatic stabilization using:
    - Method 1: getElectrostatic_stabilization() with multipole_order=1 (monopole only, field-based)
    - Method 2: getElectrostatic_stabilization() with multipole_order=2 (adds dipole term, field-based)
    - Method 3: get_Electrostatic_stabilization_dipole() (dipole-field interaction only)
    - Method 4: get_Electrostatic_stabilization_tensor() order=1 (monopole, tensor)
    - Method 5: get_Electrostatic_stabilization_tensor() order=2 (+ charge-dipole + dipole-dipole, tensor)
    - Method 6: get_Electrostatic_stabilization_tensor() order=3 (+ quadrupoles, tensor)
    """

    print("=" * 80)
    print("MULTIPOLE EXPANSION COMPARISON TEST")
    print("=" * 80)
    print()

    # ===========================
    # SETUP - MODIFY FOR YOUR SYSTEM
    # ===========================

    # List of folders containing calculations
    list_of_folders = []
    for item in os.listdir():
        if os.path.isdir(item):
            list_of_folders.append(item)

    # If you want to test on specific folders, uncomment and modify:
    # list_of_folders = ['AXEDEN']

    if not list_of_folders:
        print("ERROR: No calculation folders found!")
        print("Please run this script from a directory containing calculation folders.")
        print("Each folder should have a /scr/ subdirectory with .molden and .xyz files.")
        return

    print(f"Found {len(list_of_folders)} calculation folder(s): {list_of_folders}")
    print()

    # Initialize (one substrate atom per folder for now)
    list_of_atoms = [0 for _ in range(len(list_of_folders))]
    folder_to_file_path = '/scr/'

    # Create Electrostatics object
    try:
        estat = Electrostatics(
            list_of_folders,
            list_of_atoms,
            folder_to_file_path,
            molden_filename='input.molden',
            xyzfilename='xyz.xyz',
            hasECP=True
        )
        print("✓ Electrostatics object created successfully")
    except Exception as e:
        print(f"ERROR creating Electrostatics object: {e}")
        return

    # Multiwfn settings - MODIFY THESE FOR YOUR SYSTEM
    multiwfn_module = "module load multiwfn/noGUI_3.7"
    multiwfn_path = '/data1/groups/HJKgroup/src/multiwfn/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn'
    atmrad_path = "/data1/groups/HJKgroup/src/multiwfn/Multiwfn_3.7_bin_Linux_noGUI/examples/atmrad"

    # Define substrate atoms - MODIFY FOR YOUR SYSTEM
    # Example: substrate is atoms 0-10, environment is everything else
    substrate_indices = [np.arange(0, 11)]  # Adjust this range for your system!

    print(f"Substrate atoms: {substrate_indices[0]}")
    print(f"Environment: all other atoms")
    print()

    # Charge partitioning scheme
    charge_type = 'Hirshfeld_I'

    print("=" * 80)
    print("RUNNING CALCULATIONS...")
    print("=" * 80)
    print()

    results = {}

    # ===========================
    # METHOD 1: Field-based, Monopole only
    # ===========================
    print("Method 1: Field-based monopole (q·V) - ORIGINAL METHOD")
    print("-" * 80)
    try:
        df1 = estat.getElectrostatic_stabilization(
            multiwfn_path=multiwfn_path,
            multiwfn_module=multiwfn_module,
            atmrad_path=atmrad_path,
            substrate_idxs=substrate_indices,
            charge_type=charge_type,
            name_dataStorage='test_method1_monopole',
            multipole_order=1
        )
        results['Method 1 (Field q·V)'] = df1['Electro_FirstOrder'].iloc[0]
        print(f"✓ Energy = {df1['Electro_FirstOrder'].iloc[0]:.4f} kcal/mol")
        print(f"  Saved to: test_method1_monopole.csv")
    except Exception as e:
        print(f"✗ ERROR: {e}")
        results['Method 1 (Field q·V)'] = None
    print()

    # ===========================
    # METHOD 2: Field-based, Dipole-field only
    # ===========================
    print("Method 2: Field-based dipole (-μ·E)")
    print("-" * 80)
    try:
        df2 = estat.get_Electrostatic_stabilization_dipole(
            multiwfn_path=multiwfn_path,
            multiwfn_module=multiwfn_module,
            atmrad_path=atmrad_path,
            substrate_idxs=substrate_indices,
            charge_type=charge_type,
            name_dataStorage='test_method2_dipole_field'
        )
        results['Method 2 (Field -μ·E)'] = df2['Dipole_Field_Energy'].iloc[0]
        print(f"✓ Dipole energy = {df2['Dipole_Field_Energy'].iloc[0]:.4f} kcal/mol")
        print(f"  (This is just the dipole term, add to Method 1 for total)")
        print(f"  Saved to: test_method2_dipole_field.csv")
    except Exception as e:
        print(f"✗ ERROR: {e}")
        results['Method 2 (Field -μ·E)'] = None
    print()

    # ===========================
    # METHOD 3: Tensor formalism, Monopole only
    # ===========================
    print("Method 3: Tensor monopole (q×q) - should match Method 1")
    print("-" * 80)
    try:
        df3 = estat.get_Electrostatic_stabilization_tensor(
            multiwfn_path=multiwfn_path,
            multiwfn_module=multiwfn_module,
            atmrad_path=atmrad_path,
            substrate_idxs=substrate_indices,
            charge_type=charge_type,
            name_dataStorage='test_method3_tensor_order1',
            multipole_order=1
        )
        results['Method 3 (Tensor q×q)'] = df3['Total_Energy_kcal_mol'].iloc[0]
        print(f"✓ Energy = {df3['Total_Energy_kcal_mol'].iloc[0]:.4f} kcal/mol")
        print(f"  Saved to: test_method3_tensor_order1.csv")
    except Exception as e:
        print(f"✗ ERROR: {e}")
        results['Method 3 (Tensor q×q)'] = None
    print()

    # ===========================
    # METHOD 4: Tensor formalism, Order 2 (monopole + dipole)
    # ===========================
    print("Method 4: Tensor order 2 (q×q + q×μ + μ×μ) - NEW COMPLETE METHOD!")
    print("-" * 80)
    print("  This includes:")
    print("    • Charge-charge (q×q)")
    print("    • Charge-dipole (q×μ) ← NEW!")
    print("    • Dipole-dipole (μ×μ) ← NEW!")
    try:
        df4 = estat.get_Electrostatic_stabilization_tensor(
            multiwfn_path=multiwfn_path,
            multiwfn_module=multiwfn_module,
            atmrad_path=atmrad_path,
            substrate_idxs=substrate_indices,
            charge_type=charge_type,
            name_dataStorage='test_method4_tensor_order2',
            multipole_order=2,
            decompose_atomwise=True  # Get per-atom breakdown
        )
        if isinstance(df4, tuple):
            df4_total, df4_atomwise = df4
            results['Method 4 (Tensor order 2)'] = df4_total['Total_Energy_kcal_mol'].iloc[0]
            print(f"✓ Total energy = {df4_total['Total_Energy_kcal_mol'].iloc[0]:.4f} kcal/mol")
            print(f"  Saved to: test_method4_tensor_order2.csv")
            print(f"  Atom-wise breakdown saved to: test_method4_tensor_order2_atomwise.csv")

            # Show top contributing atoms
            df4_atomwise_sorted = df4_atomwise.sort_values('Energy_Contribution_kcal_mol', ascending=False)
            print(f"\n  Top 5 contributing substrate atoms:")
            for idx, row in df4_atomwise_sorted.head(5).iterrows():
                print(f"    Atom {int(row['Atom_Index']):3d} ({row['Atom_Symbol']:2s}): "
                      f"{row['Energy_Contribution_kcal_mol']:8.4f} kcal/mol")
        else:
            results['Method 4 (Tensor order 2)'] = df4['Total_Energy_kcal_mol'].iloc[0]
            print(f"✓ Total energy = {df4['Total_Energy_kcal_mol'].iloc[0]:.4f} kcal/mol")
    except Exception as e:
        print(f"✗ ERROR: {e}")
        results['Method 4 (Tensor order 2)'] = None
    print()

    # ===========================
    # METHOD 5: Tensor formalism, Order 3 (+ quadrupoles)
    # ===========================
    print("Method 5: Tensor order 3 (includes quadrupoles)")
    print("-" * 80)
    print("  This includes all order 2 terms PLUS:")
    print("    • Charge-quadrupole (q×Q) ← NEW!")
    print("    • Dipole-quadrupole (μ×Q) ← NEW!")
    print("    • Quadrupole-quadrupole (Q×Q) ← NEW!")
    try:
        df5 = estat.get_Electrostatic_stabilization_tensor(
            multiwfn_path=multiwfn_path,
            multiwfn_module=multiwfn_module,
            atmrad_path=atmrad_path,
            substrate_idxs=substrate_indices,
            charge_type=charge_type,
            name_dataStorage='test_method5_tensor_order3',
            multipole_order=3
        )
        results['Method 5 (Tensor order 3)'] = df5['Total_Energy_kcal_mol'].iloc[0]
        print(f"✓ Total energy = {df5['Total_Energy_kcal_mol'].iloc[0]:.4f} kcal/mol")
        print(f"  Saved to: test_method5_tensor_order3.csv")
    except Exception as e:
        print(f"✗ ERROR: {e}")
        results['Method 5 (Tensor order 3)'] = None
    print()

    # ===========================
    # SUMMARY
    # ===========================
    print("=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print()

    # Create results DataFrame
    results_df = pd.DataFrame([
        {
            'Method': 'Field q·V (monopole)',
            'Energy (kcal/mol)': results.get('Method 1 (Field q·V)'),
            'Terms included': 'q×V (substrate charge × env potential)',
            'Comments': 'Original field-based method'
        },
        {
            'Method': 'Field -μ·E (dipole)',
            'Energy (kcal/mol)': results.get('Method 2 (Field -μ·E)'),
            'Terms included': '-μ·E (substrate dipole × env field)',
            'Comments': 'Dipole correction only'
        },
        {
            'Method': 'Tensor q×q',
            'Energy (kcal/mol)': results.get('Method 3 (Tensor q×q)'),
            'Terms included': 'q×q (direct charge-charge)',
            'Comments': 'Should match Method 1'
        },
        {
            'Method': 'Tensor order 2',
            'Energy (kcal/mol)': results.get('Method 4 (Tensor order 2)'),
            'Terms included': 'q×q + q×μ + μ×μ',
            'Comments': '✨ RECOMMENDED: Complete dipole treatment'
        },
        {
            'Method': 'Tensor order 3',
            'Energy (kcal/mol)': results.get('Method 5 (Tensor order 3)'),
            'Terms included': 'All above + q×Q + μ×Q + Q×Q',
            'Comments': 'Most accurate (includes quadrupoles)'
        }
    ])

    print(results_df.to_string(index=False))
    print()

    # Calculate differences
    if results.get('Method 1 (Field q·V)') and results.get('Method 4 (Tensor order 2)'):
        monopole_energy = results['Method 1 (Field q·V)']
        full_order2 = results['Method 4 (Tensor order 2)']
        higher_order_contribution = full_order2 - monopole_energy

        print("KEY INSIGHT:")
        print("-" * 80)
        print(f"Monopole-only energy:              {monopole_energy:10.4f} kcal/mol")
        print(f"Full order-2 energy (tensor):      {full_order2:10.4f} kcal/mol")
        print(f"Higher-order contribution (q×μ+μ×μ): {higher_order_contribution:10.4f} kcal/mol")
        print(f"Relative importance:               {abs(higher_order_contribution/monopole_energy*100):10.2f}%")
        print()

        if abs(higher_order_contribution) > 0.5:
            print("⚠️  WARNING: Higher-order terms are significant!")
            print("    The monopole-only approximation may not be accurate.")
            print("    Consider using the tensor formalism with multipole_order=2 or 3.")
        else:
            print("✓ Higher-order terms are small.")
            print("  Monopole approximation appears reasonable for this system.")

    print()
    print("=" * 80)
    print("TEST COMPLETE")
    print("=" * 80)
    print()
    print("Output files generated:")
    print("  • test_method1_monopole.csv")
    print("  • test_method2_dipole_field.csv")
    print("  • test_method3_tensor_order1.csv")
    print("  • test_method4_tensor_order2.csv")
    print("  • test_method4_tensor_order2_atomwise.csv (per-atom breakdown)")
    print("  • test_method5_tensor_order3.csv")
    print()
    print("Recommendation:")
    print("  For most accurate results, use get_Electrostatic_stabilization_tensor()")
    print("  with multipole_order=2 (includes dipole-dipole interactions)")
    print()


if __name__ == "__main__":
    compare_multipole_methods()
