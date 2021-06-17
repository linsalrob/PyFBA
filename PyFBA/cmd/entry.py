"""
An entry point for the PyFBA command
"""

import os
import sys
import argparse
import PyFBA
from .gapfill_from_roles import gapfill_from_roles

def full_help():
    """
    Just return the help text
    :return: The help
    """

    return f"""
Welcome to PyFBA version {PyFBA.__version__}

Please use one of these commands with their appropriate flags. Use pyfba <command> -h for more help

gapfill_roles\tGapfill Flux Balance Analysis from a list of functional roles
fluxes\tGiven a set of reactions that form a model, report the fluxes through those reactions


help\tThis help menu

    """

def run():
    """
    Run the appropriate pyfba command
    """

    parser = argparse.ArgumentParser(description='Run Flux Balance Analysis')
    parser.add_argument('command', help='the pyfba command you would like to run', required=True)
    args = parser.parse_args()

    if args['command'] == 'help':
        print(full_help())
        sys.exit(0)

    if args['command'] == 'gapfill_roles':
        gapfill_from_roles()


