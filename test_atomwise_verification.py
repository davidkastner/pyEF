#!/usr/bin/env python3
"""
Standalone verification script to check if atom-wise E-field contributions sum to total.
This script directly implements the calculation without requiring full package imports.
"""
import numpy as np
import pandas as pd

def test_monopole_decomposition():
    """
    Test the mathematical relationship for monopole E-field decomposition.

    For a bond between atoms A and B:
    - E_total = (1/2) * (E_A + E_B) · bond_vec
    - E_atomwise[i] = (1/2) * (E_A_from_i + E_B_from_i) · bond_vec

    We verify: sum(E_atomwise) == E_total
    """

    print("=" * 80)
    print("Testing Atom-Wise E-field Decomposition")
    print("=" * 80)
    print()

    # Use test fixture data
    charge_file = "pyef/tests/fixtures/sample_1/ChargesHirshfeld.txt"
    xyz_file = "pyef/tests/fixtures/sample_1/final_optim.xyz"

    # Physical constants
    COULOMB_CONSTANT = 8.9875517923e9  # N⋅m²/C²
    ELEMENTARY_CHARGE = 1.602176634e-19  # C
    ANGSTROM_TO_M = 1e-10  # m/Å
    VM_TO_VA = 1e10  # Convert V/m to V/Å

    # Read charge data
    try:
        df_charges = pd.read_csv(charge_file, sep=r'\s+',
                                 names=["Atom", 'x', 'y', 'z', "charge"])
        print(f"Loaded charges from: {charge_file}")
        print(f"Number of atoms: {len(df_charges)}")
        print()
    except Exception as e:
        print(f"Error loading charge file: {e}")
        return False

    # Define test bond (atoms 0 and 1)
    atom_A_idx = 0
    atom_B_idx = 1

    print(f"Testing bond: {atom_A_idx} - {atom_B_idx}")
    print(f"Atom A: {df_charges.iloc[atom_A_idx]['Atom']} at ({df_charges.iloc[atom_A_idx]['x']:.4f}, {df_charges.iloc[atom_A_idx]['y']:.4f}, {df_charges.iloc[atom_A_idx]['z']:.4f})")
    print(f"Atom B: {df_charges.iloc[atom_B_idx]['Atom']} at ({df_charges.iloc[atom_B_idx]['x']:.4f}, {df_charges.iloc[atom_B_idx]['y']:.4f}, {df_charges.iloc[atom_B_idx]['z']:.4f})")
    print()

    # Calculate bond vector
    pos_A = np.array([df_charges.iloc[atom_A_idx]['x'],
                     df_charges.iloc[atom_A_idx]['y'],
                     df_charges.iloc[atom_A_idx]['z']])
    pos_B = np.array([df_charges.iloc[atom_B_idx]['x'],
                     df_charges.iloc[atom_B_idx]['y'],
                     df_charges.iloc[atom_B_idx]['z']])

    bond_vec_unnorm = pos_A - pos_B
    bond_len = np.linalg.norm(bond_vec_unnorm)
    bond_vec = bond_vec_unnorm / bond_len

    print(f"Bond length: {bond_len:.6f} Å")
    print(f"Bond vector (normalized): [{bond_vec[0]:.6f}, {bond_vec[1]:.6f}, {bond_vec[2]:.6f}]")
    print()

    # Calculate E-field at atom A from all other atoms
    E_A = np.zeros(3)
    E_A_atomwise = []

    for i in range(len(df_charges)):
        if i == atom_A_idx:
            E_A_atomwise.append(np.zeros(3))
            continue

        r_vec = ANGSTROM_TO_M * np.array([
            df_charges.iloc[i]['x'] - pos_A[0],
            df_charges.iloc[i]['y'] - pos_A[1],
            df_charges.iloc[i]['z'] - pos_A[2]
        ])
        r = np.linalg.norm(r_vec)

        # E-field contribution from atom i
        E_contrib = -COULOMB_CONSTANT * ELEMENTARY_CHARGE * df_charges.iloc[i]['charge'] * r_vec / (r**3)
        E_A += E_contrib
        E_A_atomwise.append(E_contrib)

    # Calculate E-field at atom B from all other atoms
    E_B = np.zeros(3)
    E_B_atomwise = []

    for i in range(len(df_charges)):
        if i == atom_B_idx:
            E_B_atomwise.append(np.zeros(3))
            continue

        r_vec = ANGSTROM_TO_M * np.array([
            df_charges.iloc[i]['x'] - pos_B[0],
            df_charges.iloc[i]['y'] - pos_B[1],
            df_charges.iloc[i]['z'] - pos_B[2]
        ])
        r = np.linalg.norm(r_vec)

        # E-field contribution from atom i
        E_contrib = -COULOMB_CONSTANT * ELEMENTARY_CHARGE * df_charges.iloc[i]['charge'] * r_vec / (r**3)
        E_B += E_contrib
        E_B_atomwise.append(E_contrib)

    # Convert to V/Å
    E_A = VM_TO_VA * E_A
    E_B = VM_TO_VA * E_B
    E_A_atomwise = [VM_TO_VA * E for E in E_A_atomwise]
    E_B_atomwise = [VM_TO_VA * E for E in E_B_atomwise]

    # Calculate total projected E-field
    E_proj_total = 0.5 * np.dot(E_A + E_B, bond_vec)

    # Calculate atom-wise projected E-field contributions
    E_proj_atomwise = []
    for i in range(len(df_charges)):
        E_proj_i = 0.5 * np.dot(E_A_atomwise[i] + E_B_atomwise[i], bond_vec)
        E_proj_atomwise.append(E_proj_i)

    # Sum atom-wise contributions
    E_proj_atomwise_sum = np.sum(E_proj_atomwise)

    # Display results
    print("-" * 80)
    print("RESULTS")
    print("-" * 80)
    print(f"Total E-field at A: [{E_A[0]:.6e}, {E_A[1]:.6e}, {E_A[2]:.6e}] V/Å")
    print(f"Total E-field at B: [{E_B[0]:.6e}, {E_B[1]:.6e}, {E_B[2]:.6e}] V/Å")
    print()
    print(f"Total projected E-field:           {E_proj_total:.10f} V/Å")
    print(f"Sum of atom-wise projections:      {E_proj_atomwise_sum:.10f} V/Å")
    print()

    difference = abs(E_proj_total - E_proj_atomwise_sum)
    relative_error = difference / abs(E_proj_total) if E_proj_total != 0 else 0

    print(f"Absolute difference:  {difference:.2e} V/Å")
    print(f"Relative error:       {relative_error:.2e}")
    print()

    # Check if they match (use relative tolerance for large values)
    abs_tolerance = 1e-10
    rel_tolerance = 1e-12

    if difference < abs_tolerance or relative_error < rel_tolerance:
        print(f"✓ PASS: Atom-wise contributions sum to total")
        print(f"  (Relative error {relative_error:.2e} < {rel_tolerance:.0e})")
        result = True
    else:
        print(f"✗ FAIL: Significant difference detected!")
        print(f"  Relative error {relative_error:.2e} exceeds tolerance {rel_tolerance:.0e}")
        result = False

    print()
    print("-" * 80)
    print("Top 10 atom-wise contributions:")
    print(f"{'Atom':<6} {'Element':<8} {'Charge':<12} {'E_proj (V/Å)':<15} {'Cumulative':<15}")
    print("-" * 80)

    cumsum = 0
    sorted_indices = np.argsort(np.abs(E_proj_atomwise))[::-1][:10]

    for idx in sorted_indices:
        cumsum += E_proj_atomwise[idx]
        print(f"{idx:<6} {df_charges.iloc[idx]['Atom']:<8} {df_charges.iloc[idx]['charge']:>11.6f} {E_proj_atomwise[idx]:>14.8f} {cumsum:>14.8f}")

    print()
    return result


if __name__ == "__main__":
    import sys
    success = test_monopole_decomposition()
    sys.exit(0 if success else 1)
