import sys

from PyFBA import lp, log_and_message


def reaction_bounds(reactions, reactions_to_run, media, uptakesecretionreactions={}, lower=-1000.0, mid=0.0,
                    upper=1000.0, verbose=False):
    """
    Set the bounds for each reaction. We set the reactions to run between
    either lower/mid, mid/upper, or lower/upper depending on whether the
    reaction runs <=, =>, or <=> respectively.

    :param reactions: The dict of all reactions we know about
    :type reactions: dict of metabolism.Reaction
    :param reactions_to_run: The sorted list of reactions to run
    :type reactions_to_run: set
    :param media: The media compounds
    :type media: set
    :param uptakesecretionreactions: The optional uptake and secretion reactions
    :type uptakesecretionreactions: Dict[str, PyFBA.metabolism.Reaction]
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
    for r in reactions_to_run:
        if r == 'BIOMASS_EQN':
            rbvals[r] = (0, 1000)
            continue
        # if we already know the bounds, eg from an SBML file
        if r in uptakesecretionreactions and uptakesecretionreactions[r].lower_bound is not None and \
                uptakesecretionreactions[r].upper_bound is not None:
            rbvals[r] = (uptakesecretionreactions[r].lower_bound, uptakesecretionreactions[r].upper_bound)
            continue

        if r in reactions and reactions[r].lower_bound is not None and reactions[r].upper_bound is not None:
            rbvals[r] = (reactions[r].lower_bound, reactions[r].upper_bound)
            continue

        importer = False # does this reaction import something?

        left_compounds = set()
        if r in uptakesecretionreactions:
            direction = uptakesecretionreactions[r].direction
            importer = True
            left_compounds = uptakesecretionreactions[r].left_compounds
        elif r in reactions:
            direction = reactions[r].direction
            # does this reaction import something [e] -> [c]
            importer = reactions[r].is_transport or reactions[r].is_uptake_secretion or reactions[r].is_input_reaction()
            left_compounds = reactions[r].left_compounds
        else:
            log_and_message(f"Did not find {r} in reactions", stderr=verbose)
            direction = "="

        # this is where we define whether our media has the components
        if importer:
            in_media = False
            # if we have external compounds that are not in the media, we don't want to run this as a media reaction
            override = False
            for c in left_compounds:
                if c.location == 'e':
                    if c in media:
                        in_media = True
                    else:
                        override = True

            if override:
                # in this case, we have some external compounds that we should not import.
                # for example, many reactions have H+[e] on the left side and H+[c] on the
                # right side, so we don't want to allow those reactions unless we want to import everything.
                # e.g. (1) H+[e] + (1) D-Alanine[e] <=> (1) H+[c] + (1) D-Alanine[c]
                in_media = False

            """
            For the bounds:
            (-1000,1000) means that the reaction can proceed in either direction
            (0, 1000) means that you can only go left to right
            (-1000, 0) means you can flux right to left
            """

            if in_media:
                rbvals[r] = (lower, upper)
                media_uptake_secretion_count += 1
            else:
                # We have now defined those bounds in external_reactions.py, and so we should not get here!
                if len(reactions[r].left_compounds) == 1 and len(reactions[r].right_compounds) == 1:
                    log_and_message(f"Probably should do something with {r} {reactions[r].equation}", stderr=verbose)
                rbvals[r] = (0.0, 0.0)
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
            rbvals[r] = (lower, mid)
            # rbvals[r] = (lower, upper)
        else:
            sys.stderr.write("DO NOT UNDERSTAND DIRECTION " + direction + " for " + r + "\n")
            rbvals[r] = (mid, upper)

    log_and_message(f"In parsing the bounds we found {media_uptake_secretion_count} media uptake " +
                    f"and secretion reactions and {other_uptake_secretion_count} other u/s reactions", stderr=verbose)

    if len(rbvals) != len(reactions_to_run):
        log_and_message(f"ERROR: We have {len(rbvals)} boundary values but {len(reactions_to_run)} reactions",
                        stderr=True, loglevel='ERROR')

    rbounds = [rbvals[r] for r in reactions_to_run]

    # set the bounds in the reactions
    for r in rbvals:
        if r in reactions:
            (reactions[r].lower_bound, reactions[r].upper_bound) = rbvals[r]

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
