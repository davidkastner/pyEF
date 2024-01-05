"""Command-line interface (CLI) entry point"""

import pyef
import click

@click.group()
def cli():
    """pyEF Command-line Interface."""
    welcome_message()

@click.command()
@click.option("--first_function", is_flag=True, help="Perform an action")
def efield(first_function):
    """Calculates electric fields"""
    click.echo("Calcute the projected electric field.")
    if first_function:
        click.echo("Downloading Glue dataset...")
        from pyef.efield import default
        

@click.command()
@click.option("--first_function", is_flag=True, help="Perform an action")
def esp(first_function,):
    """Calculates electrostatic potential."""
    click.echo("Tools and workflows for making predictions")
    if first_function:
        click.echo("Benchmarking AlphaFold against DockQ...")
        from pyef.esp import benchmark_alphafold


def welcome_message():
    print("\n.-------------------------.")
    print("| WELCOME TO THE PYEF CLI |")
    print(".-------------------------.")
    print("Default programmed actions for the pyEF package.")
    print("GitHub: https://github.com/davidkastner/pyef")
    print("Documentation: https://pyef.readthedocs.io\n")


cli.add_command(efield)
cli.add_command(esp)

if __name__ == '__main__':
    cli()
