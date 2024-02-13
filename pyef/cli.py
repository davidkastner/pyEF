"""Command-line interface (CLI) entry point"""

import os
import pyef
import click

def welcome():
    """Print first to welcome the user while it waits to load the modules"""
    
    print("\n     ╔═══════════════════════════╗")
    print("     ║ .-----------------------. ║")
    print("     ║ |       ┌──────┐        | ║")
    print("     ║ |       │  ┌───┘        | ║")
    print("     ║ |       │  └───┐        | ║")
    print("     ║ |       │  ┌───┘        | ║")
    print("     ║ |       │  ├┬┬┬┐        | ║")
    print("     ║ |       └──┴┴┴┴┘        | ║")
    print("     ║ |                       | ║")
    print("     ║ |    WELCOME TO PYEF    | ║")
    print("     ║ '-----------------------' ║")
    print("     ╚═══════════════════════════╝\n")

    print("Default programmed actions for the pyEF package.")
    print("GitHub: https://github.com/davidkastner/pyef")
    print("Documentation: https://pyef.readthedocs.io")
    print("• Command for electric field analysis: pyef ef -i path/to/jobs.in")
    print("• Command for electrostatic analysis: pyef esp -i path/to/jobs.in")
    print("• Example annotated job file: github.com/davidkastner/pyEF/tree/main/demo/jobs.in")
    print("• Update the settings.ini file in pyef.resources, especially the `nthreads` variable\n")

# Welcome even if no flags
welcome()

@click.group()
def cli():
    """CLI entry point"""
    pass

@cli.command()
@click.option("-i", "--input", "job_path", required=True, type=str, help="Path to pyEF job file.")
def ef(job_path):
    """Analyzes electric fields"""
    click.echo("Importing dependencies...")
    from pyef.run import main
    from pyef.manage import parse_job_batch_file

    geom_flag = False   # Perform a geometry check
    esp_flag = False    # Perform analysis of electrostatics

    # Now using the job_path provided by the -i option
    jobs, metal_indices, bond_indices = parse_job_batch_file(job_path)
    job_name = os.path.splitext(job_path)[0]
    pyef.run.main(job_name, jobs, metal_indices, bond_indices, geom_flag, esp_flag)

# THE ESP SECTION IS UNDERCONSTRUCTION AND MAY NOT WORK        
@cli.command()
@click.option("-i", "--input", "job_path", required=True, type=str, help="Path to pyEF job file.")
def esp(job_path):
    """Analyzes electrostatic potential"""
    click.echo("Performing electrostatic analysis.")
    from pyef.run import main
    from pyef.manage import parse_job_batch_file

    geom_flag = False   # Perform a geometry check
    esp_flag = True     # Perform analysis of electrostatics

    jobs, metal_indices, bond_indices = parse_job_batch_file(job_path)
    pyef.run.main(jobs, geom_flag, esp_flag, metal_indices)


if __name__ == '__main__':
    cli()
