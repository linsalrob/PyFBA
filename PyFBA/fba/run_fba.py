import sys

from PyFBA import lp
import PyFBA

def run_fba(modeldata, reactions_to_run, media, biomass_equation, uptake_secretion=None, verbose=False):
    """
    Run an fba for a set of data. We required the reactions object,
    a list of reactions to run, the media, and the biomass_equation equation.

    With all of these we run the fba and return:

    :param modeldata: the model seed object that includes compounds and reactions
    :type modeldata: PyFBA.model_seed.ModelSeed
    :param reactions_to_run: the reactions to run
    :type reactions_to_run: set
    :param media: An array of compound.Compound objects representing the media
    :type media: set
    :param biomass_equation: The biomass_equation equation
    :type biomass_equation: network.reaction.Reaction
    :param uptake_secretion: A hash of uptake and secretion reactions that should be added to the model.
    Calculated if not provided.
    :type uptake_secretion: dict of Reaction
    :param verbose: Print more output
    :type verbose: bool
    :return: which type of linear resolution, the output value of the model, whether the model grew
    :rtype: (str, float, bool)

    """

    cp, rc = PyFBA.fba.create_stoichiometric_matrix(reactions_to_run, modelseed, media,
                                                               biomass_equation, uptake_secretion, verbose=verbose)

    rbvals = PyFBA.fba.reaction_bounds(modelseed.reactions, rc, media)
    PyFBA.fba.compound_bounds(cp)

    if verbose:
        sys.stderr.write("Length of the media: {}\n".format(len(media)))
        sys.stderr.write("Number of reactions to run: {}\n".format(len(reactions_to_run)))
        sys.stderr.write("Number of compounds in SM: {}\n".format(len(cp)))
        sys.stderr.write("Number of reactions in SM: {}\n".format(len(rc)))
        sys.stderr.write("Revised number of total reactions: {}\n".format(len(modelseed.reactions)))
        sys.stderr.write("Number of total compounds: {}\n".format(len(modelseed.compounds)))
        sys.stderr.write("SMat dimensions: {} x {}\n".format(len(cp), len(rc)))

    status, value = PyFBA.lp.solve()
    
    growth = False
    if value > 1:
        growth = True
    
    return status, value, growth

