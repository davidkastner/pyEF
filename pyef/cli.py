"""Command-line interface (CLI) entry point"""

import pyef
import click

@click.group()
def cli():
    """pyEF Command-line Interface."""
    welcome_message()

@click.command()
@click.option("--run", is_flag=True, help="Analyze electric fields.")
def ef(run):
    """Analyzes electric fields"""
    click.echo("Calcute the projected electric field.")
    if run:
        click.echo("Importing dependencies...")
        from pyef.run import main
        from pyef.analysis import Electrostatics
        from pyef.geometry import ErrorAnalysis

        geom_flag = False   # Perform a geometry check
        esp_flag = False    # Perform analysis of electrostatics
        job_paths = input("   > Paths to jobs separated by commas: ")
        jobs = [job.strip() for job in job_paths.split(",")]
        metal_indices = input("   > Indices of your metals, separated by commas: ")
        metal_indices = [metal.strip() for metal in metal_indices.split(",")]
        pyef.run.main(jobs, geom_flag, esp_flag, metal_indices)
        
@click.command()
@click.option("--run", is_flag=True, help="Perform an action")
def esp(run):
    """Analyzes electrostatic potential"""
    click.echo("Performing electrostatic analysis.")
    if run:
        click.echo("Importing dependencies...")
        from pyef.run import main
        from pyef.analysis import Electrostatics
        from pyef.geometry import ErrorAnalysis

        geom_flag = False   # Perform a geometry check
        esp_flag = True    # Perform analysis of electrostatics
        job_paths = input("   > Paths to jobs separated by commas: ")
        jobs = [job.strip() for job in job_paths.split(",")]
        metal_indices = input("   > Indices of your metals, separated by commas: ")
        metal_indices = [metal.strip() for metal in metal_indices.split(",")]
        pyef.run.main(jobs, geom_flag, esp_flag, metal_indices)

def welcome_message():
    print("\n.-------------------------.")
    print("| WELCOME TO THE pyEF CLI |")
    print(".-------------------------.")
    print("Default programmed actions for the pyEF package.")
    print("GitHub: https://github.com/davidkastner/pyef")
    print("Documentation: https://pyef.readthedocs.io")
    print("Command for electric field analysis: pyef ef --run")
    print("Command for electrostatic analysis: pyef esp --run\n")

cli.add_command(ef)
cli.add_command(esp)

if __name__ == '__main__':
    cli()
