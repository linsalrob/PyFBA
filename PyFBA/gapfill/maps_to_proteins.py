import sys

import PyFBA

__author__ = "Rob Edwards"


def suggest_reactions_without_proteins(reactions, verbose=False):
    """
    Add the reactions that don't have any proteins associated with them.

    This is a bad idea to add all of these since there are about 30,000
    and it will break your computer if you try and solve it with FBA


    :param reactions: our reactions dictionary from parsing the model seed
    :type reactions: dict
    :param verbose: add additional output
    :type verbose: bool
    :return: a set of reactions that could be added to test for growth
    :rtype: set
    """

    # Here are the reactions that do not have a protein associated with them
    nopegs = PyFBA.filters.reactions_with_no_proteins(reactions, verbose=verbose)

    if verbose:
        sys.stderr.write("WITHOUT proteins suggesting {} additional reactions\n".format(len(nopegs)))

    return nopegs


def suggest_reactions_with_proteins(reactions, verbose=False):
    """
    Suggest those reactions that we don't have but that have proteins
    associated with them.

    :param reactions: our reactions dictionary from parsing the model seed
    :type enz_react: dict
    :param verbose: add additional output
    :type verbose: bool
    :return: a set of reactions that could be added to test for growth
    :rtype: set
    """

    peg_rx = PyFBA.filters.reactions_with_proteins(reactions, verbose=verbose)

    if verbose:
        sys.stderr.write("WITH proteins suggesting {} additional reactions\n".format(len(peg_rx)))

    return peg_rx
