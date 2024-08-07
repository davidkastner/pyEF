"""Command-line interface (CLI) entry point"""

import os
import pyef
import yaml
import click

def welcome():
    """Print first to welcome the user while it waits to load the modules"""
    print("\n")
    print("             ╔═══════════════════════╗")
    print("             ║        ┌──────┐       ║")
    print("             ║        │  ┌───┘       ║")
    print("             ║        │  └───┐       ║")
    print("             ║        │  ┌───┘       ║")
    print("             ║        │  ├┬┬┬┐       ║")
    print("             ║        └──┴┴┴┴┘       ║")
    print("             ║                       ║")
    print("             ║    WELCOME TO PYEF    ║")
    print("             ╚═══════════╗╔══════════╝")
    print("                 ╔═══════╝╚═══════╗                 ")
    print("                 ║ THE KULIK LAB  ║                 ")
    print("                 ╚═══════╗╔═══════╝                 ")
    print("  ╔══════════════════════╝╚══════════════════════╗  ")
    print("  ║  Code: github.com/davidkastner/pyef          ║  ")
    print("  ║  Docs: pyef.readthedocs.io                   ║  ")
    print("  ║     - Usage: pyef -c pyef.in                 ║  ")
    print("  ║  Example: github.com/davidkastner/pyef/demo  ║  ")
    print("  ╚══════════════════════════════════════════════╝  \n")

# Welcome even if no flags
welcome()

# Read in the configuration yaml file
def read_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

@click.group(invoke_without_command=True)
@click.option("--config", "-c", type=click.Path(exists=True), help="Path to the configuration YAML file")
@click.pass_context
def cli(ctx, config):
    """CLI entry point"""
    if ctx.invoked_subcommand is None and config:
        ctx.invoke(run, config=config)
    elif ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())

@click.command()
@click.option("--config", "-c", required=True, type=click.Path(exists=True), help="Path to the configuration YAML file")
def run(config):
    from pyef.run import main
    from pyef.manage import parse_job_batch_file

    config_data = read_config(config)
    input = config_data.get('input', '')
    dielectric = config_data.get('dielectric', 1)
    run_ef = config_data.get('ef', False)
    run_esp = config_data.get('esp', False)
    run_geometry_check = config_data.get('esp', False)
    multiwfn_module = config_data.get('multiwfn_module', False)
    multiwfn_path = config_data.get('multiwfn_path', False)
    atmrad_path = config_data.get('atmrad_path', False)

    if run_ef:
        """Analyzes electric fields"""
        click.echo("Performing electric field analysis")

        jobs, metal_indices, bond_indices = parse_job_batch_file(input)
        job_name = os.path.splitext(input)[0]
        pyef.run.main(job_name, jobs, metal_indices, bond_indices, dielectric, run_geometry_check, run_esp, multiwfn_module, multiwfn_path, atmrad_path)

    # THE ESP SECTION IS UNDERCONSTRUCTION AND MAY NOT WORK  
    if run_esp:      
        """Analyzes electrostatic potential"""
        click.echo("Performing electrostatic analysis")

        jobs, metal_indices, bond_indices = parse_job_batch_file(input)
        pyef.run.main(job_name, jobs, metal_indices, bond_indices, dielectric, run_geometry_check, run_esp, multiwfn_module, multiwfn_path, atmrad_path)

cli.add_command(run)

if __name__ == '__main__':
    cli()
