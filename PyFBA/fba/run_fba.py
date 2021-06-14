from PyFBA import lp, log_and_message
import PyFBA


def run_fba(modeldata, reactions_to_run, media, biomass_equation, uptake_secretion=None, verbose=False):
    """
    Run an fba for a set of data. We required the reactions object,
    a list of reactions to run, the media, and the biomass_equation equation.

    With all of these we run the fba and return:

    :param modeldata: the model seed object that includes compounds and reactions
    :type modeldata: PyFBA.model_seed.ModelData
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

    cp, rc, upsr = PyFBA.fba.create_stoichiometric_matrix(reactions_to_run, modeldata, media,
                                                          biomass_equation, uptake_secretion, verbose=verbose)

    rbvals = PyFBA.fba.reaction_bounds(modeldata.reactions, rc, media, verbose=verbose)
    PyFBA.fba.compound_bounds(cp)

    if verbose:
        log_and_message(f"Length of the media: {len(media)}", stderr=verbose)
        log_and_message(f"Number of reactions to run: {len(reactions_to_run)}", stderr=verbose)
        log_and_message(f"Number of compounds in SM: {len(cp)}", stderr=verbose)
        log_and_message(f"Number of reactions in SM: {len(rc)}", stderr=verbose)
        log_and_message(f"Number of uptake/secretion reactions {len(upsr)}", stderr=verbose)
        log_and_message(f"SMat dimensions: {len(cp)} x {len(rc)}", stderr=verbose)

    status, value = PyFBA.lp.solve()

    growth = False
    if value > 1:
        growth = True

    return status, value, growth
