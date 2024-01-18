"""Functions for managing files need for pyEF"""

def parse_job_batch_file(file_path):
    """
    Parse a job batch file to extract job paths.

    This function reads a file specified by `file_path`.
    Each line in the file is expected to contain a path to a job.
    The function processes the file line by line, extracting the job paths.
    Lines that are either empty or start with '#' (comments) are ignored.
    Leading and trailing whitespace on each line is also removed.

    Parameters
    ----------
    file_path: str
        A string specifying the path to the job batch file.

    Returns
    -------
    jobs: list
        A list of job paths extracted from the file. Each job path is a string.

    Example
    -------
    >>> parse_job_batch_file('job_batch.txt')
    ['path/to/job1', 'path/to/job2', ...]
    
    """

    with open(file_path, 'r') as file:
        jobs = []
        for line in file:
            # Strip leading and trailing whitespace
            line = line.strip()

            # Ignore empty lines and lines starting with '#'
            if line and not line.startswith('#'):
                jobs.append(line)

        return jobs
