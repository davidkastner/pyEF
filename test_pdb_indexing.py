#!/usr/bin/env python3
"""
Test script to verify PDB indexing fix.
Creates known test data, writes to PDB, reads back and verifies.
"""

import sys
import os
import numpy as np

# Add pyef to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pyef.geometry import Visualize

print("=" * 80)
print("Testing PDB Indexing Fix")
print("=" * 80)

# Use the test fixture
test_dir = "pyef/tests/fixtures/sample_1"
xyz_file = os.path.join(test_dir, "final_optim.xyz")
charge_file = os.path.join(test_dir, "ChargesHirshfeld.txt")

if not os.path.exists(xyz_file):
    print(f"ERROR: Test file not found: {xyz_file}")
    sys.exit(1)

# Count atoms in XYZ file
with open(xyz_file, 'r') as f:
    n_atoms = int(f.readline().strip())

print(f"\nTest file: {xyz_file}")
print(f"Number of atoms: {n_atoms}")

# Create known test data
# Use a simple pattern: atom i gets value i/100.0
# This makes it easy to verify the mapping
test_bfactors = np.array([i / 100.0 for i in range(n_atoms)])

print(f"\nCreating test B-factors:")
print(f"  b_col[0] = {test_bfactors[0]:.4f} (should go to XYZ atom 0, PDB atom 1)")
print(f"  b_col[1] = {test_bfactors[1]:.4f} (should go to XYZ atom 1, PDB atom 2)")
print(f"  b_col[10] = {test_bfactors[10]:.4f} (should go to XYZ atom 10, PDB atom 11)")
print(f"  b_col[20] = {test_bfactors[20]:.4f} (should go to XYZ atom 20, PDB atom 21)")

# Write to PDB using the fixed makePDB function
pdb_output = "test_indexing_fix.pdb"
print(f"\nWriting to PDB: {pdb_output}")

try:
    viz = Visualize(xyz_file)
    viz.makePDB(charge_file, test_bfactors, pdb_output)
    print("  ✓ PDB created successfully")
except Exception as e:
    print(f"  ✗ ERROR creating PDB: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Read back the PDB and verify
print(f"\nVerifying PDB contents:")

with open(pdb_output, 'r') as f:
    lines = f.readlines()

# Extract B-factors by PDB atom number
pdb_data = {}
for line in lines:
    if line.startswith('HETATM') or line.startswith('ATOM'):
        atom_num = int(line[6:11].strip())
        element = line[12:16].strip()
        b_factor = float(line[60:66].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        pdb_data[atom_num] = {
            'element': element,
            'b_factor': b_factor,
            'coords': (x, y, z)
        }

print(f"  Found {len(pdb_data)} atoms in PDB")

# Verify the mapping
print(f"\nVerifying B-factor mapping:")
print(f"  {'PDB Atom':<10} {'Expected':<12} {'Actual':<12} {'Match':<8}")
print(f"  {'-'*50}")

errors = []
test_atoms = [1, 2, 3, 10, 11, 20, 21, 24, 25]  # Test various atoms

for pdb_atom_num in test_atoms:
    if pdb_atom_num not in pdb_data:
        print(f"  {pdb_atom_num:<10} MISSING!")
        errors.append(f"PDB atom {pdb_atom_num} missing")
        continue

    # PDB atom N should have b_factor from b_col[N-1] (0-indexed XYZ atom N-1)
    xyz_atom_idx = pdb_atom_num - 1
    expected = test_bfactors[xyz_atom_idx]
    actual = pdb_data[pdb_atom_num]['b_factor']

    # PDB format truncates to 2 decimal places
    match = abs(expected - actual) < 0.01
    status = "✓" if match else "✗"

    print(f"  {pdb_atom_num:<10} {expected:<12.4f} {actual:<12.2f} {status:<8}")

    if not match:
        errors.append(f"PDB atom {pdb_atom_num}: expected {expected:.4f}, got {actual:.2f}")

print(f"\n" + "=" * 80)
if errors:
    print(f"FAILED: Found {len(errors)} error(s):")
    for error in errors:
        print(f"  - {error}")
    sys.exit(1)
else:
    print("SUCCESS: All B-factors mapped correctly!")
    print("  The indexing fix is working properly.")
    sys.exit(0)
