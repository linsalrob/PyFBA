import sys


def reactions_with_no_proteins(reactions, verbose=False):
    """
    Figure out which reactions in our set have no proteins associated with them.

    :param reactions: The reactions dictionary
    :type reactions: dict
    :param verbose: prints out how many reactions have no proteins out of the total
    :type verbose: bool
    :return: a set of reaction ids that have no proteins associated with them.
    :rtype: set
    """

    nopegs = set()
    for r in reactions:
        if reactions[r].number_of_enzymes() == 0:
            nopegs.add(r)

    if verbose:
        sys.stderr.write("REACTIONS WITH NO PROTEINS: {} reactions have no pegs associated ".format(len(nopegs)) +
                         "with them (out of {} reactions)\n".format(len(reactions)))

    return nopegs


def reactions_with_proteins(reactions, verbose=False):
    """
    Figure out which reactions in our set have proteins associated with them.

    :param reactions: The reactions dictionary
    :type reactions: dict
    :param verbose: prints out how many reactions have no proteins out of the total
    :type verbose: bool
    :return: a set of reaction ids that have proteins associated with them.
    :rtype: set
    """

    pegs = set()
    for r in reactions:
        if reactions[r].number_of_enzymes() != 0:
            pegs.add(r)

    if verbose:
        sys.stderr.write("REACTIONS WITH PROTEINS: {} reactions have pegs associated ".format(len(pegs)) +
                         "with them (out of {} reactions)\n".format(len(reactions)))

    return pegs
