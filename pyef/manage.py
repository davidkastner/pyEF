"""Functions for managing files need for pyEF"""

def parse_job_batch_file(file_path):
    """
    Parse a CSV file and extract specific columns as lists and tuples.

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
    column_pairs = []

    with open(file_path, 'r') as file:
        for line in file:
            # Skip empty lines and comments '#'
            if line.strip() == '' or line.strip().startswith('#'):
                continue

            columns = [col.strip() for col in line.strip().split(',')]
            # Extracting and appending data to respective lists
            jobs.append(columns[0])
            metal_index = int(columns[1])
            bonded_atom_index = int(columns[2])
            metal_indices.append(metal_index)
            column_pairs.append([(metal_index, bonded_atom_index)])

    return jobs, metal_indices, column_pairs