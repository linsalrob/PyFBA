import copy
import PyFBA
from PyFBA import log_and_message


def uptake_and_secretion_reactions(model_compounds, media):
    """
    Figure out which compounds can be taken up from the media and/or secreted into the media. We provide an endless
    reaction for these which allows them to be taken up and/or secreted without affecting the rest of the stoichiometric
    matrix

    We also add a reaction for biomass_equation.

    We will eventually set the bounds for these reactions, such that if the compound is in the media, the
    bounds are (-1000,1000) [i.e. the compound can flow into the media] whereas if the compounds are not
    in the media, the bounds are (0,1000) [ie. the compound can flow out of the media, but not into it]

    :param model_compounds: A set of the compounds that are in this model
    :type model_compounds: set[PyFBA.metabolism.CompoundWithLocation]
    :param media: the media we want to grow on
    :type media: set[PyFBA.metabolism.CompoundWithLocation]
    :return: A hash of new uptake and secretion reactions we need to add to the model
    :rtype: hash
    """

    uptake_sec_reactions = {}
    count = 0
    for c in model_compounds:
        if c.location == 'e' or c.name == 'Biomass':
            # this is an uptake or secretion reaction
            # we need to add a new compound like this with a false location
            us_leftside = c
            us_rightside = copy.deepcopy(us_leftside)
            # The uptake and secretion compounds typically have reaction bounds that allow them to be consumed
            # (i.e. diffuse away from the cell) but not produced. However, our media components can also increase
            # in concentration (i.e. diffuse to the cell) and thus the bounds are set higher. Whenever you change the
            # growth media, you also need to adjust the reaction bounds to ensure that the media can be consumed!
            # the b is for boundary and is secretion away from the cell
            us_rightside.location = 'b'
            us_rightside.uptake_secretion = True
            # this is similar name that they use in the model seed
            # us_reaction = Reaction('us_001', "EX_" + us_leftside.model_seed_id + "_" + us_leftside.location + "0")
            # but we normally use a different name
            us_reaction_id = f"upsr_{count}"
            us_reaction = PyFBA.metabolism.Reaction(us_reaction_id, f"UPTAKE_SECRETION_REACTION {count}")
            count += 1
            us_reaction.equation = '(1) + ' + str(us_leftside) + " <=> (1) + " + str(us_rightside)
            us_reaction.add_left_compounds({us_leftside})
            us_reaction.set_left_compound_abundance(us_leftside, 1)
            us_reaction.add_right_compounds({us_rightside})
            us_reaction.set_right_compound_abundance(us_rightside, 1)
            us_reaction.set_direction('=')
            us_reaction.is_uptake_secretion = True
            us_leftside.add_reactions({us_reaction})
            us_rightside.add_reactions({us_reaction})
            # Here we set reaction bounds. If the compound is in the media, we let it flow freely
            # otherwise we only let it diffuse away
            if c in media:
                us_reaction.lower_bound = -1000
                us_reaction.upper_bound = 1000
            else:
                us_reaction.lower_bound = 0
                us_reaction.upper_bound = 1000
            uptake_sec_reactions[us_reaction_id] = us_reaction

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
        if r.startswith('upsr_'):
            toremove.add(r)

    for r in toremove:
        reactions.pop(r)
    return reactions
