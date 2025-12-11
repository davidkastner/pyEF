#!/usr/bin/env python3
"""
Verification script to check if atom-wise E-field contributions sum to total E-field projection.
"""
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyef.analysis import Electrostatics

def verify_atomwise_decomposition():
    """Test that atom-wise E-field contributions sum to total projected E-field."""

    print("=" * 80)
    print("Verifying atom-wise E-field decomposition")
    print("=" * 80)

    # Setup paths using test fixture
    fixture_path = "tests/fixtures/sample_1/"
    xyz_file = os.path.join(fixture_path, "final_optim.xyz")
    charge_file = os.path.join(fixture_path, "ChargesHirshfeld.txt")

    # Check if files exist
    if not os.path.exists(xyz_file):
        print(f"ERROR: XYZ file not found: {xyz_file}")
        return False
    if not os.path.exists(charge_file):
        print(f"ERROR: Charge file not found: {charge_file}")
        return False

    print(f"Using XYZ file: {xyz_file}")
    print(f"Using charge file: {charge_file}")
    print()

    # Create a minimal config for testing
    config = {
        'dielectric': 1,
        'includePtChgs': False,
        'changeDielectBoundBool': False
    }

    # Create Electrostatics instance
    es = Electrostatics(
        xyzfilename='final_optim.xyz',
        folder_to_file_path='',
        lst_of_folders=[fixture_path],
        lst_of_tmcm_idx=[0],  # Assume first atom is metal
        config=config
    )

    # Define test bond - let's use atoms 0 and 1
    bond_indices = [(0, 1)]

    # Get all atoms (excluding none for this test)
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    n_atoms = int(lines[0].strip())
    all_lines = list(range(n_atoms))

    print(f"Number of atoms: {n_atoms}")
    print(f"Testing bond: {bond_indices[0]}")
    print()

    # Test monopole mode (simpler, easier to verify)
    print("-" * 80)
    print("Testing MONOPOLE mode")
    print("-" * 80)

    try:
        [E_projected, bonded_atoms, bond_idx, bond_lens, E_proj_atomwise,
         E_proj_atomwise_list] = es.bondEfield(
            bond_indices, xyz_file, charge_file, all_lines,
            bool_multipole=False, df_ptchg=None
        )

        print(f"Total E-field projection: {E_projected[0]:.10f} V/Å")
        print(f"Bond length: {bond_lens[0]:.6f} Å")
        print(f"Bonded atoms: {bonded_atoms[0]}")
        print()

        # Sum the atom-wise contributions
        atomwise_sum = np.sum(E_proj_atomwise_list[0])
        print(f"Sum of atom-wise contributions: {atomwise_sum:.10f} V/Å")
        print()

        # Calculate difference
        difference = abs(E_projected[0] - atomwise_sum)
        relative_error = difference / abs(E_projected[0]) if E_projected[0] != 0 else 0

        print(f"Difference: {difference:.2e} V/Å")
        print(f"Relative error: {relative_error:.2e}")
        print()

        # Check if they match within numerical precision
        tolerance = 1e-10
        if difference < tolerance:
            print(f"✓ PASS: Atom-wise contributions sum to total (within {tolerance} tolerance)")
            monopole_pass = True
        else:
            print(f"✗ FAIL: Atom-wise contributions do NOT sum to total")
            print(f"  Expected difference < {tolerance}, got {difference}")
            monopole_pass = False

    except Exception as e:
        print(f"✗ ERROR in monopole test: {e}")
        import traceback
        traceback.print_exc()
        monopole_pass = False

    print()
    print("=" * 80)

    # Show breakdown of first few atom contributions
    if E_proj_atomwise_list:
        print("\nFirst 10 atom-wise contributions:")
        print(f"{'Atom':<6} {'Contribution (V/Å)':<20} {'Cumulative Sum (V/Å)':<20}")
        print("-" * 50)
        cumsum = 0
        for i, contrib in enumerate(E_proj_atomwise_list[0][:10]):
            cumsum += contrib
            print(f"{i:<6} {contrib:>19.10f} {cumsum:>19.10f}")
        if len(E_proj_atomwise_list[0]) > 10:
            print(f"... ({len(E_proj_atomwise_list[0]) - 10} more atoms)")
        print()

    return monopole_pass

if __name__ == "__main__":
    success = verify_atomwise_decomposition()
    sys.exit(0 if success else 1)
