import sys

from PyFBA import lp, log_and_message


def reaction_bounds(reactions, reactions_with_upsr, media, lower=-1000.0, mid=0.0, upper=1000.0, verbose=False):
    """
    Set the bounds for each reaction. We set the reactions to run between
    either lower/mid, mid/upper, or lower/upper depending on whether the
    reaction runs <=, =>, or <=> respectively.

    :param reactions: The dict of all reactions we know about
    :type reactions: dict of metabolism.Reaction
    :param reactions_with_upsr: The sorted list of reactions to run
    :type reactions_with_upsr: set
    :param media: The media compounds
    :type media: set
    :param lower: The default lower bound
    :type lower: float
    :param mid: The default mid value (typically 0)
    :type mid: float
    :param upper: The default upper bound
    :type upper: float
    :return: A dict of the reaction ID and the tuple of bounds
    :rtype: dict
    """

    rbvals = {}
    media_uptake_secretion_count = 0
    other_uptake_secretion_count = 0
    for r in reactions_with_upsr:
        if r == 'BIOMASS_EQN':
            rbvals[r] = (mid, upper)
            continue

        # if we already know the bounds, eg from an SBML file or from our uptake/secretion reactions
        if reactions[r].lower_bound != None and reactions[r].upper_bound != None:
            rbvals[r] = (reactions[r].lower_bound, reactions[r].upper_bound)
            continue

        if r in reactions:
            direction = reactions[r].direction
        else:
            sys.stderr.write("Did not find {} in reactions\n".format(r))
            direction = "="

        """
        RAE 16/6/21
        We no longer use this block to check for media components. Instead, we us the uptake_and_secretion_reactions
        in external_reactions.py to do so.
        
        We assume that if you provide uptake_and_secretion_reactions you have already culled them for the media, though
        perhaps we should add a test for that.
        
        """

        if False and (reactions[r].is_uptake_secretion or reactions[r].is_transport or reactions[r].is_input_reaction()):
            in_media = False
            override = False # if we have external compounds that are not in the media, we don't want to run this as a media reaction
            for c in reactions[r].left_compounds:
                if c.location == 'e':
                    if c in media:
                        in_media = True
                    else:
                        override = True
            # in this case, we have some external compounds that we should not import.
            # for example, H+ is used to translocate things
            if override:
                in_media = False

            if in_media:
                # This is what I think it should be:
                rbvals[r] = (lower, upper)
                #rbvals[r] = (0.0, upper)
                media_uptake_secretion_count += 1
            else:
                rbvals[r] = (0.0, upper)
                #rbvals[r] = (lower, upper)
                other_uptake_secretion_count += 1
            continue

        if direction == "=":
            # This is what I think it should be:
            rbvals[r] = (lower, upper)
            # rbvals[r] = (mid, upper)
        elif direction == ">":
            # This is what I think it should be:
            rbvals[r] = (mid, upper)
            # rbvals[r] =  (lower, upper)
        elif direction == "<":
            # This is what I think it should be:
            # rbvals[r] = (lower, mid)
            rbvals[r] = (lower, upper)
        else:
            sys.stderr.write("DO NOT UNDERSTAND DIRECTION " + direction + " for " + r + "\n")
            rbvals[r] = (mid, upper)

    if verbose:
        sys.stderr.write("In parsing the bounds we found {} media uptake ".format(media_uptake_secretion_count) +
                         "and secretion reactions and {} other u/s reactions\n".format(other_uptake_secretion_count))

    rbounds = [rbvals[r] for r in reactions_with_upsr]
    for r in reactions_with_upsr:
        if r in reactions:
            reactions[r].lower_bound, reactions[r].upper_bound = rbvals[r]

    lp.col_bounds(rbounds)
    return rbvals


def compound_bounds(cp, lower=0, upper=0):
    """
    Impose constraints on the compounds. These constraints limit what
    the variation of each compound can be and is essentially 0 for
    most compounds except those that are in the media or otherwise
    external.

    This is the zero flux vector.

    Parameters:
        cp: the list of compound ids
        lower: the default lower value
        upper: the default upper value
    """

    cbounds = [(lower, upper) for c in cp]
    cbvals = {c: (lower, upper) for c in cp}

    lp.row_bounds(cbounds)
    return cbvals
