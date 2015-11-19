import sys

from lp import solve
from metabolism import Compound
from . import compound_bounds
from . import create_stoichiometric_matrix
from . import reaction_bounds


def run_fba(compounds, reactions, reactions_to_run, media, biomass_equation, verbose=False):
    """
    Run an fba for a set of data. We required the reactions object,
    a list of reactions to run, the media, and the biomass equation.

    With all of these we run the fba and return:

    :param compounds: The dict of all compounds
    :type compounds: dict
    :param reactions: The dict of all reactions
    :type reactions: dict
    :param reactions_to_run: the reactions to run
    :type reactions_to_run: set
    :param media: An array of compound.Compound objects representing the media
    :type media: set
    :param biomass_equation: The biomass equation
    :type biomass_equation: network.reaction.Reaction
    :param verbose: Print more output
    :type verbose: bool
    :return: which type of linear resolution, the output value of the model, whether the model grew
    :rtype: (str, float, bool)

    """

    cp, rc, reactions = create_stoichiometric_matrix(reactions_to_run, reactions, compounds, media, biomass_equation, verbose=False)

    rbvals = reaction_bounds(reactions, rc, media, verbose=verbose)
    compound_bounds(cp)

    if verbose:
        sys.stderr.write("Length of the media: {}\n".format(len(media)))
        sys.stderr.write("Number of reactions to run: {}\n".format(len(reactions_to_run)))
        sys.stderr.write("Number of compounds in SM: {}\n".format(len(cp)))
        sys.stderr.write("Number of reactions in SM: {}\n".format(len(rc)))
        sys.stderr.write("Revised number of total reactions: {}\n".format(len(reactions)))
        sys.stderr.write("Number of total compounds: {}\n".format(len(compounds)))
        sys.stderr.write("SMat dimensions: {} x {}\n".format(len(cp), len(rc)))

    status, value = solve()
    
    growth = False
    if value > 1:
        growth = True
    
    return status, value, growth

