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
    run_estab = config_data.get('estab', False)
    run_geometry_check = config_data.get('geometry_check', False)
    multiwfn_module = config_data.get('multiwfn_module', False)
    multiwfn_path = config_data.get('multiwfn_path', False)
    atmrad_path = config_data.get('atmrad_path', False)

    # Additional configuration options
    charge_types = config_data.get('charge_types', ['Hirshfeld_I'])
    multipole_bool = config_data.get('multipole', True)
    use_multipole = config_data.get('use_multipole', False)
    decompose_atomwise = config_data.get('decompose_atomwise', False)
    substrate_idxs = config_data.get('substrate_idxs', None)
    env_idxs = config_data.get('env_idxs', None)
    multipole_order = config_data.get('multipole_order', 2)

    # Parse job batch file
    jobs, metal_indices, bond_indices = parse_job_batch_file(input)
    job_name = os.path.splitext(input)[0]

    # Run requested analyses
    if run_ef or run_esp or run_estab:
        pyef.run.main(
            job_name, jobs, metal_indices, bond_indices, dielectric,
            run_geometry_check, run_esp, run_ef, run_estab,
            multiwfn_module, multiwfn_path, atmrad_path,
            charge_types, multipole_bool, use_multipole,
            decompose_atomwise, substrate_idxs, env_idxs, multipole_order
        )
    else:
        click.echo("No analysis type specified. Please set ef, esp, or estab to true in config file.")

cli.add_command(run)

if __name__ == '__main__':
    cli()
