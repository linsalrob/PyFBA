import copy

import PyFBA


def uptake_and_secretion_reactions(model_compounds, compounds):
    """
    Figure out which compounds can be taken up from the media and/or secreted into the media. We provide an endless
    reaction for these which allows them to be taken up and/or secreted without affecting the rest of the stoichiometric
    matrix

    We also add a reaction for biomass_equation

    :param model_compounds: A set of the identifiers of all the compounds we have identified so far
    :type model_compounds: set
    :param compounds: the dict of all the compounds
    :type compounds: dict
    :return: A hash of new uptake and secretion reactions we have identified
    :rtype: hash
    """

    uptake_sec_reactions = {}
    for c in model_compounds:
        if compounds[c].location == 'e' or compounds[c].name == 'Biomass':
            # this is an uptake or secretion reaction
            # we need to add a new compound like this with a false location
            us_leftside = compounds[c]
            us_rightside = copy.copy(us_leftside)
            us_rightside.location = 'b'
            # this is similar name that they use in the model seed
            # us_reaction = Reaction("EX_" + us_leftside.model_seed_id + "_" + us_leftside.location + "0")
            # but we normally use a different name
            us_reaction = PyFBA.metabolism.Reaction("UPTAKE_SECRETION_REACTION " + us_leftside.model_seed_id)
            us_reaction.equation = '(1) + ' + str(us_leftside) + " <=> (1) + " + str(us_rightside)
            us_reaction.add_left_compounds({us_leftside})
            us_reaction.set_left_compound_abundance(us_leftside, 1)
            us_reaction.add_right_compounds({us_rightside})
            us_reaction.set_right_compound_abundance(us_rightside, 1)
            us_reaction.set_direction('=')
            us_reaction.is_uptake_secretion = True
            uptake_sec_reactions[str(us_reaction)] = us_reaction
            us_leftside.add_reactions({us_reaction})
            us_rightside.add_reactions({us_reaction})

    return uptake_sec_reactions


def remove_uptake_and_secretion_reactions(reactions):
    """
    Remove all the uptake and secretion reactions added to a model, eg. when you are running multiple simulations.
    :param reactions: The reactions dict
    :type reactions: dict
    :return: The enzymes, compounds, and reactions data structure
    :rtype: dict
    """

    toremove = set()
    for r in reactions:
        if r.startswith("UPTAKE_SECRETION_REACTION"):
            toremove.add(r)

    for r in toremove:
        reactions.pop(r)
    return reactions
