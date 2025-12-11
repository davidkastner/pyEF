"""Functions for managing files needed for pyEF"""

def parse_job_batch_file(file_path):
    """
    Parse a CSV file and extract specific columns as lists and tuples.
    The input file allows Python style comments on any line.

    Parameters
    ----------
    file_path : str
        The file path of the CSV file to be parsed.

    Returns
    -------
    analysis_types : list of str
        List of analysis types ('ef', 'esp', 'estab', or combinations like 'ef+esp').
    molden_paths : list of str
        List of paths to .molden files.
    xyz_paths : list of str
        List of paths to .xyz files.
    metal_indices : list of int or None
        List containing metal atom indices.
        None if metal indices are not provided.
    bond_indices : list of tuples or None
        List of bond tuples for each job.
        None if bond indices are not provided.

    Notes
    -----
    New format (required):
    - analysis_type, path_to_molden, path_to_xyz, [bond_tuples or metal_index]
    - Examples:
      - ef, /path/to/optim.molden, /path/to/optim.xyz, (25, 26), (25, 27)
      - esp, /path/to/optim.molden, /path/to/optim.xyz, 30
      - estab, /path/to/optim.molden, /path/to/optim.xyz
      - ef+esp, /path/to/optim.molden, /path/to/optim.xyz, 35
    """
    import re

    analysis_types = []
    molden_paths = []
    xyz_paths = []
    metal_indices = []
    bond_indices = []

    # Valid analysis type keywords
    valid_analysis = {'ef', 'esp', 'estab'}

    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, 1):
            # Skip empty lines and comments that start with '#'
            if line.strip() == '' or line.strip().startswith('#'):
                continue

            # Remove comments from the line
            line = line.split('#')[0].strip()

            # Skip the line if it's empty after removing the comment
            if line == '':
                continue

            # Split line into columns (initially by comma)
            columns = [col.strip() for col in line.split(',')]

            # Validate minimum number of columns
            if len(columns) < 3:
                raise ValueError(
                    f"Line {line_num}: Invalid format. Expected at least 3 columns: "
                    f"analysis_type, path_to_molden, path_to_xyz [, atom_indices]. "
                    f"Got {len(columns)} columns."
                )

            # Extract analysis type (column 1)
            analysis_type = columns[0].lower().strip()
            analysis_parts = set(analysis_type.replace('+', ' ').split())
            if not analysis_parts.issubset(valid_analysis) or not analysis_parts:
                raise ValueError(
                    f"Line {line_num}: Invalid analysis type '{columns[0]}'. "
                    f"Must be one of: {', '.join(valid_analysis)} or combinations like 'ef+esp'."
                )
            analysis_types.append(columns[0])  # Keep original case

            # Extract molden path (column 2)
            molden_path = columns[1].strip()
            if not molden_path:
                raise ValueError(f"Line {line_num}: Missing molden file path in column 2.")
            molden_paths.append(molden_path)

            # Extract xyz path (column 3)
            xyz_path = columns[2].strip()
            if not xyz_path:
                raise ValueError(f"Line {line_num}: Missing xyz file path in column 3.")
            xyz_paths.append(xyz_path)

            # Parse remaining columns for metal indices and/or bond tuples
            remainder = ','.join(columns[3:])
            has_tuples = '(' in remainder and ')' in remainder

            if has_tuples:
                # Parse tuple format for bonds
                job_bonds = []
                metal_idx = None

                tuple_pattern = r'\(\s*(\d+)\s*,\s*(\d+)\s*\)'
                matches = re.findall(tuple_pattern, remainder)

                for match in matches:
                    atom1, atom2 = int(match[0]), int(match[1])
                    job_bonds.append((atom1, atom2))
                    # Set metal_idx to first atom of first bond (for compatibility)
                    if metal_idx is None:
                        metal_idx = atom1

                metal_indices.append(metal_idx)
                bond_indices.append(job_bonds)

            elif len(columns) > 3 and columns[3]:
                # Single metal index
                try:
                    metal_index = int(columns[3])
                    metal_indices.append(metal_index)

                    # Check if there's a bonded atom index
                    if len(columns) > 4 and columns[4]:
                        bonded_atom_index = int(columns[4])
                        bond_indices.append([(metal_index, bonded_atom_index)])
                    else:
                        bond_indices.append([])
                except ValueError:
                    raise ValueError(
                        f"Line {line_num}: Invalid atom index '{columns[3]}'. "
                        f"Expected integer or tuple format like (1, 2)."
                    )
            else:
                # No additional data
                metal_indices.append(None)
                bond_indices.append([])

    # Return None for metal/bond indices if all are empty/None
    if metal_indices and all(m is None for m in metal_indices):
        metal_indices = None
    if bond_indices and all(not b for b in bond_indices):
        bond_indices = None

    return analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices
