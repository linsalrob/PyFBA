import sys

from lp import col_bounds
from lp import row_bounds


def reaction_bounds(reactions, reactions_to_run, media, lower=-1000.0, mid=0.0, upper=1000.0):
    """
    Set the bounds for each reaction. We set the reactions to run between
    either lower/mid, mid/upper, or lower/upper depending on whether the
    reaction runs <=, =>, or <=> respectively.

    :param reactions: The dict of all reactions we know about
    :type reactions: dict of metabolism.Reaction
    :param reactions_to_run: The sorted list of reactions to run
    :type reactions_to_run: list
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

    for r in reactions_to_run:
        if r in reactions:
            direction = reactions[r].direction
        elif r == 'bme':
            direction = '>'
        else:
            sys.stderr.write("Did not find " + r + " in reactions\n")
            direction = "="

        # this is where we define whether our media has the components
        if r != 'bme' and reactions[r].is_uptake_secretion:
            in_media = False
            for c in reactions[r].left_compounds:
                if c in media:
                    in_media = True
            if in_media:
                rbvals[r] = (lower, upper)
            else:
                rbvals[r] = (0.0, upper)
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

    rbounds = [rbvals[r] for r in reactions_to_run]
    col_bounds(rbounds)
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

    row_bounds(cbounds)
    return cbvals
