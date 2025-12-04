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
    jobs : list of str
        List containing entries from the first column of the CSV.
    metal_indices : list of int or None
        List containing entries from the second column of the CSV as integers.
        None if metal indices are not provided (for E-field/stabilization without ESP).
    bond_indices : list of tuples or None
        List of tuples, each containing the second and third column entries of the CSV.
        None if bond indices are not provided.

    Notes
    -----
    Supported formats:
    - 1 column: job_path (for stabilization analysis)
    - 2 columns: job_path, metal_index (for ESP analysis)
    - 3 columns: job_path, metal_index, bonded_atom_index (for E-field analysis)
    """

    jobs = []
    metal_indices = []
    bond_indices = []

    with open(file_path, 'r') as file:
        for line in file:
            # Skip empty lines and comments that start with '#'
            if line.strip() == '' or line.strip().startswith('#'):
                continue

            # Remove comments from the line
            line = line.split('#')[0].strip()

            # Skip the line if it's empty after removing the comment
            if line == '':
                continue

            columns = [col.strip() for col in line.split(',')]

            # Extract job path (always required)
            jobs.append(columns[0])

            # Extract metal and bond indices if provided
            if len(columns) >= 2 and columns[1]:
                metal_index = int(columns[1])
                metal_indices.append(metal_index)

                if len(columns) >= 3 and columns[2]:
                    bonded_atom_index = int(columns[2])
                    bond_indices.append([(metal_index, bonded_atom_index)])
                else:
                    bond_indices.append([])
            else:
                metal_indices.append(None)
                bond_indices.append([])

    # Return None for metal/bond indices if all are empty/None
    if all(m is None for m in metal_indices):
        metal_indices = None
    if all(not b for b in bond_indices):
        bond_indices = None

    return jobs, metal_indices, bond_indices
