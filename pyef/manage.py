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
    metal_indices : list of int
        List containing entries from the second column of the CSV as integers.
    column_pairs : list of tuples
        List of tuples, each containing the second and third column entries of the CSV, with the second column as an integer.

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
            
            # Extracting and appending data to respective lists
            jobs.append(columns[0])
            metal_index = int(columns[1])
            bonded_atom_index = int(columns[2])
            metal_indices.append(metal_index)
            bond_indices.append([(metal_index, bonded_atom_index)])

    return jobs, metal_indices, bond_indices
