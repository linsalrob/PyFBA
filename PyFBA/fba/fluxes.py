import os
import sys
from PyFBA import lp
import PyFBA


def reaction_fluxes(verbose=False):
    """
    Return the reaction fluxes from the solved FBA model.

    :param verbose: Print more output
    :type verbose: bool
    :return: A dict of reaction ID and flux through that reaction
    :rtype: dict of str and float
    """

    return lp.col_primal_hash()
