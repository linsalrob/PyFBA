import sys

from . import suggest_reactions_with_proteins


def compound_probability(reactions, reactions2run, cutoff=0, rxn_with_proteins=True, verbose=False):
    """
    Identify a set of reactions that you should add to your model for growth based on the probability for the reaction
    to run left to right and right to left.

    The probability is basically the fraction of compounds that are present.

    If you set cutoff to zero we calculate the minimum coverage of the compounds in the reactions already in the model
    and use that to determine which other reactions should be added based on the observation that they have a similar
    proportion of compounds already in the network.

    :param reactions: our reactions dict
    :type reactions: dict
    :param reactions2run: The current set of reactions that we will run
    :type reactions2run: set
    :param cutoff: is a minimum probability that must be exceeded for the reaction to be proposed
    :type cutoff: float
    :param rxn_with_proteins: limits to just those reactions that have proteins
    :type rxn_with_proteins: bool
    :param verbose: print more output
    :type verbose: bool
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set
    """

    # get a set of all compounds in all the reactions
    cpds = set()
    for r in reactions2run:
        cpds.update(reactions[r].all_compounds())

    if verbose:
        sys.stderr.write("There are " + str(len(cpds)) + " compounds before assigning probabilities\n")

    r2rscores = {'left': [], 'right': []}

    rtest = set(reactions.keys())
    if rxn_with_proteins:
        # just consider those reactions with proteins
        rtest = suggest_reactions_with_proteins(reactions)

    for r in rtest:
        # score our probability as a similarity so 1 is good, 0 is bad
        # similarity = (1.0 * len(t)-len(t.difference(s)))/len(t)
        if len(reactions[r].left_compounds) > 0:
            left = reactions[r].left_compounds
            lsim = (1.0 * len(left) - len(left.difference(cpds))) / len(left)
        else:
            lsim = 0

        if len(reactions[r].right_compounds) > 0:
            right = reactions[r].right_compounds
            rsim = (1.0 * len(right) - len(right.difference(cpds))) / len(right)
        else:
            rsim = 0

        reactions[r].pLR = lsim
        reactions[r].pRL = rsim

        if r in reactions2run:
            r2rscores['left'].append(lsim)
            r2rscores['right'].append(rsim)

    # what is the average score for our reactions2run?
    existing_scores = [
        1.0 * sum(r2rscores['left']) / len(r2rscores['left']),
        1.0 * sum(r2rscores['right']) / len(r2rscores['right'])
    ]

    if cutoff == 0:
        cutoff = min(existing_scores)

    suggested_reactions = set()
    for r in reactions:
        if r in reactions2run:
            continue
        if reactions[r].pLR > cutoff:
            suggested_reactions.add(r)
        elif reactions[r].pRL > cutoff:
            suggested_reactions.add(r)

    if verbose:
        sys.stderr.write("Average score for reactions we are running:\n")
        sys.stderr.write("Left: " + str(existing_scores[0]) + "\n")
        sys.stderr.write("Right: " + str(existing_scores[1]) + "\n")
        sys.stderr.write("Using a cutoff of " + str(cutoff) + "\n")
        sys.stderr.write("Proposing " + str(len(suggested_reactions)) + " reactions\n")

    return suggested_reactions
