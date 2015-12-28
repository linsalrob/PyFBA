

def limit_reactions_by_compound(reactions, reactions2run, suggestions, max_rcts=50):
    """
    Limit the reactions in suggestions based on the compounds present in
    the reactions in reactions2run and the number of reactions that each
    compound is associated with.

    We need to have < max_rcts reactions per compound for it to be
    considered. This is to avoid things like H2O that have a lot of
    connections

    :param reactions: The reactions dict
    :type reactions: dict
    :param reactions2run: our base set of reactions that we will run
    :type reactions2run: set
    :param suggestions: the reactions we are considering adding
    :type suggestions: set
    :param max_rcts: the maximum number of reactions per compound
    :type max_rcts: int
    :return: a set of reactions which is those members of suggestions that meet our criteria
    :rtype: set

    """

    cpd = {}
    for r in reactions2run:
        for c in reactions[r].all_compounds():
            cpd[str(c)] = cpd.get(str(c), 0) + 1

    keep = set()
    for r in suggestions:
        for c in reactions[r].all_compounds():
            if str(c) in cpd and (cpd[str(c)] < max_rcts):
                keep.add(r)

    keep.difference_update(reactions2run)

    return keep
