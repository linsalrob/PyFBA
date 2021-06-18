"""
An entry point for the PyFBA command
"""

import os
import sys
import argparse
import PyFBA
from .citation import cite_me_please
from .fluxes import measure_fluxes
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
version\tPrint the version and exit

citations\tGet the citations for PyFBA and the work that it is built upon
    """


def run():
    """
    Run the appropriate pyfba command
    """


    if sys.argv[1] == 'help' or sys.argv[1] == '-h' or sys.argv[1] == '--help': 
        print(full_help())
        sys.exit(0)
    elif 'version' in sys.argv[1] or '-v' in sys.argv[1]:
        print(PyFBA.__version__)
        sys.exit(0)
    elif 'citation' in sys.argv[1] or 'cite' in sys.argv[1]:
        cite_me_please()
        sys.exit(0)
    elif sys.argv[1] == 'gapfill_roles':
        gapfill_from_roles()
        sys.exit(0)
    elif sys.argv[1] == 'fluxes':
        measure_fluxes()
        sys.exit(0)
    else:
        sys.stderr.write(f"Sorry. Don't understand {args.command}.")
        sys.stderr.write(full_help())
        sys.exit(0)

