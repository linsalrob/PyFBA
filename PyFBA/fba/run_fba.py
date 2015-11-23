import sys

from PyFBA import lp
import PyFBA

def run_fba(compounds, reactions, reactions_to_run, media, biomass_equation, uptake_secretion={}, verbose=False):
    """
    Run an fba for a set of data. We required the reactions object,
    a list of reactions to run, the media, and the biomass_equation equation.

    With all of these we run the fba and return:

    :param uptake_secretion: A hash of uptake and secretion reactions that should be added to the model. Calculated if not provided.
    :type uptake_secretion: dict of Reaction
    :param compounds: The dict of all compounds
    :type compounds: dict
    :param reactions: The dict of all reactions
    :type reactions: dict
    :param reactions_to_run: the reactions to run
    :type reactions_to_run: set
    :param media: An array of compound.Compound objects representing the media
    :type media: set
    :param biomass_equation: The biomass_equation equation
    :type biomass_equation: network.reaction.Reaction
    :param verbose: Print more output
    :type verbose: bool
    :return: which type of linear resolution, the output value of the model, whether the model grew
    :rtype: (str, float, bool)

    """

    cp, rc, reactions = PyFBA.fba.create_stoichiometric_matrix(reactions_to_run, reactions, compounds, media, biomass_equation,
                                                     uptake_secretion, verbose=False)

    rbvals = PyFBA.fba.reaction_bounds(reactions, rc, media, verbose=verbose)
    PyFBA.fba.compound_bounds(cp)

    if verbose:
        sys.stderr.write("Length of the media: {}\n".format(len(media)))
        sys.stderr.write("Number of reactions to run: {}\n".format(len(reactions_to_run)))
        sys.stderr.write("Number of compounds in SM: {}\n".format(len(cp)))
        sys.stderr.write("Number of reactions in SM: {}\n".format(len(rc)))
        sys.stderr.write("Revised number of total reactions: {}\n".format(len(reactions)))
        sys.stderr.write("Number of total compounds: {}\n".format(len(compounds)))
        sys.stderr.write("SMat dimensions: {} x {}\n".format(len(cp), len(rc)))

    status, value = PyFBA.lp.solve()
    
    growth = False
    if value > 1:
        growth = True
    
    return status, value, growth

