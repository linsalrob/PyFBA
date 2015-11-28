import sys

__author__ = "Rob Edwards"


def suggest_by_compound(compounds, reactions, reactions2run, max_reactions, verbose=False):
    """
    Identify a set of reactions that you should add to your model for growth because they contain orphan compounds

    This is a slightly different approach to suggesting by compound.
    We look for "orphan" compounds that only have a few connections to
    other reactions, and then add those to our network to see what
    happens.

    Note: that our notion of two compounds being equal normally
    includes their location (and str(c) includes the location).
    However, we probably have several instances of:

        cpd [e] -> cpd [c]
        and
        cpd [c] -> products

    These should be considered to be the same, and we probably
    don't want to consider external compounds any way

    :param reactions: The reactions dictionary
    :type reactions: dict
    :param compounds: The compounds dictionary
    :type compounds: dict
    :param reactions2run: The set of reactions that we will already run
    :type reactions2run: set
    :param max_reactions: The maximum number of reactions that a compound can be associated with. Avoids, eg. H2O
    :type max_reactions: int
    :param verbose: Print more output
    :type verbose: bool
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set

    """

    cpd = {}
    for r in reactions2run:
        for c in reactions[r].all_compounds():
            cpd[str(c)] = cpd.get(str(c), 0) + 1

    ikeep = set()
    ekeep = set()

    external = 0
    internal = 0
    for c in compounds:
        if str(compounds[c]) in cpd and cpd[str(compounds[c])] <= max_reactions:
            if compounds[c].location == 'e':
                external += 1
                ekeep.update(compounds[c].all_reactions())
            else:
                internal += 1
                ikeep.update(compounds[c].all_reactions())

    if verbose:
        sys.stdout.write("{} | {} | {} | {} | {}\n".format(max_reactions, internal, len(ikeep), external, len(ekeep)))

    ikeep = {r for r in ikeep if r in reactions and r not in reactions2run}
    ekeep = {r for r in ekeep if r in reactions and r not in reactions2run}

    return ikeep
