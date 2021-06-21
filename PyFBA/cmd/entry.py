"""
An entry point for the PyFBA command
"""

import sys
import PyFBA
from . import *


def full_help():
    """
    Just return the help text
    :return: The help
    """

    return f"""
Welcome to PyFBA version {PyFBA.__version__}

Please use one of these commands with their appropriate flags. Use pyfba <command> -h for more help

fba\tGiven a file with a set of reactions, run an FBA on that set of reactions
fluxes\tGiven a set of reactions that form a model, report the fluxes through those reactions

to_reactions\tConvert a set of functional roles or feature names to a list of reactions
gapfill_roles\tGapfill Flux Balance Analysis from a list of functional roles
multiple_media\tGapfill a set of reactions with multiple media where the organism can/can not grow
create_gaps\tThe opposite of gapfill: Given a media and a set of reactions, reduce them to the smallest set that can run

reactions_to_roles\tGiven a file with a set of reactions, print the roles that implement those reactions

media\tList the names of all the predefined media

help\tThis help menu
version\tPrint the version and exit

citations\tGet the citations for PyFBA and the work that it is built upon
    """


def run():
    """
    Run the appropriate pyfba command
    """

    if len(sys.argv) == 1 or sys.argv[1] == 'help' or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print(full_help())
        sys.exit(0)
    elif 'version' in sys.argv[1] or '-v' in sys.argv[1]:
        print(PyFBA.__version__)
        sys.exit(0)
    elif 'citation' in sys.argv[1] or 'cite' in sys.argv[1]:
        cite_me_please()
    elif sys.argv[1] == 'multiple_media':
        gapfill_multiple_media()
    elif sys.argv[1] == 'gapfill_roles':
        gapfill_from_roles()
    elif sys.argv[1] == 'fluxes':
        measure_fluxes()
    elif sys.argv[1] == 'media':
        list_media()
    elif sys.argv[1] == 'reactions_to_roles':
        convert_reactions_to_roles()
    elif sys.argv[1] == 'fba':
        run_the_fba()
    elif sys.argv[1] == 'to_reactions':
        to_reactions()
    elif sys.argv[1] == 'create_gaps':
        create_reaction_gaps()
    else:
        sys.stderr.write(f"Sorry. Don't understand {sys.argv[1]}.")
        sys.stderr.write(full_help())
        sys.exit(0)

