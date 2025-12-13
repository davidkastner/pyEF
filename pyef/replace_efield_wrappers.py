#!/usr/bin/env python3
"""Script to replace E-field function implementations with simple wrappers."""

import re

# Read the file
with open('analysis.py', 'r') as f:
    content = f.read()

# Define the wrapper for getEfield_acrossBond
efield_acrossbond_wrapper = '''    def getEfield_acrossBond(self, charge_type, Efielddata_filename, multiwfn_path, multipole_bool, input_bond_indices=[], dielectric=1):
        """Legacy wrapper for getEfield() with bond-specific E-field calculation.

        DEPRECATED: Use getEfield() directly instead.

        This function now calls the unified getEfield() method with appropriate parameters
        to maintain backward compatibility.
        """
        return self.getEfield(
            charge_types=charge_type,
            Efielddata_filename=Efielddata_filename,
            multiwfn_path=multiwfn_path,
            multipole_bool=multipole_bool,
            input_bond_indices=input_bond_indices,
            auto_find_bonds=False,
            decompose_atomwise=True,
            visualize=None,
            dielectric=dielectric
        )'''

# Define the wrapper for getEfield_decomposable
efield_decomposable_wrapper = '''    def getEfield_decomposable(self, charge_type, Efielddata_filename, multiwfn_path, multipole_bool, input_bond_indices=[], dielectric=1):
        """Legacy wrapper for getEfield() with decomposable E-field calculation.

        DEPRECATED: Use getEfield() directly instead.

        This function now calls the unified getEfield() method with appropriate parameters
        to maintain backward compatibility.
        """
        return self.getEfield(
            charge_types=charge_type,
            Efielddata_filename=Efielddata_filename,
            multiwfn_path=multiwfn_path,
            multipole_bool=multipole_bool,
            input_bond_indices=input_bond_indices,
            auto_find_bonds=False,
            decompose_atomwise=True,
            visualize=None,
            dielectric=dielectric
        )'''

# Define the wrapper for getEFieldMultipole
efield_multipole_wrapper = '''    def getEFieldMultipole(self, Efield_data_filename, multiwfn_path, input_bond_indices=[], excludeAtoms=[], polarization_scheme='Hirshfeld_I'):
        """Legacy wrapper for getEfield() with multipole E-field calculation.

        DEPRECATED: Use getEfield() with multipole_bool=True instead.

        This function now calls the unified getEfield() method with appropriate parameters
        to maintain backward compatibility.
        """
        # Set excludeAtoms if provided
        if excludeAtoms:
            self.config['excludeAtomfromEcalc'] = excludeAtoms

        return self.getEfield(
            charge_types=polarization_scheme,
            Efielddata_filename=Efield_data_filename,
            multiwfn_path=multiwfn_path,
            multipole_bool=True,
            input_bond_indices=input_bond_indices,
            auto_find_bonds=False,
            decompose_atomwise=False,
            visualize=None,
            dielectric=1
        )'''

# Pattern to match getEfield_acrossBond from definition to just before next function
pattern1 = re.compile(
    r'(    def getEfield_acrossBond\(self,.*?\):)\s*\'\'\'.*?\'\'\'.*?'
    r'(?=\n#this should work both for polarizable)',
    re.DOTALL
)

# Pattern to match getEfield_decomposable
pattern2 = re.compile(
    r'(    def getEfield_decomposable\(self,.*?\):)\s*\'\'\'.*?\'\'\'.*?'
    r'(?=\n\n\n    def getEFieldMultipole)',
    re.DOTALL
)

# Pattern to match getEFieldMultipole
pattern3 = re.compile(
    r'(    def getEFieldMultipole\(self,.*?\):)\s*\'\'\'.*?\'\'\'.*?'
    r'(?=\n    def getElectrostatic_stabilization)',
    re.DOTALL
)

# Replace each function
content = pattern1.sub(efield_acrossbond_wrapper + '\n', content)
content = pattern2.sub(efield_decomposable_wrapper + '\n\n', content)
content = pattern3.sub(efield_multipole_wrapper + '\n', content)

# Write back
with open('analysis.py', 'w') as f:
    f.write(content)

print("✅ Successfully replaced all 3 E-field functions with wrappers!")
print("File updated: analysis.py")
