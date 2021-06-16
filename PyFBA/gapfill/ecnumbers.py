import os
import sys
import re

import PyFBA
from PyFBA import log_and_message


def suggest_reactions_using_ec(roles, modeldata, reactions2run, verbose=False):
    """
    Identify a set of reactions that you should add to your model for growth based on the EC numbers
    that may be found in the role names.

    :param roles: A set of all roles to search for EC numbers.
    :type roles: set
    :param modeldata: our modeldata datastructure from parsing the seed models
    :type reactions: PyFBA.model_seed.ModelData
    :param reactions2run: set of reactions that  we are going to run
    :type reactions2run: set
    :param verbose: add additional output
    :type verbose: bool
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set
    """

    ec_to_reactions = {}
    for r in modeldata.enzymes:
        for e in modeldata.enzymes[r].ec_number:
            if e not in ec_to_reactions:
                ec_to_reactions[e] = set()
            ec_to_reactions[e].add(r)
    for r in modeldata.reactions:
        for e in modeldata.reactions[r].ec_numbers:
            if e not in ec_to_reactions:
                ec_to_reactions[e] = set()
            ec_to_reactions[e].add(r)

    # Find all EC numbers in the list of roles
    suggested_reactions = set()
    for role in roles:
        # Extract each EC number
        for ec in re.findall(r"[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+", role):
            # Check if we know what that EC number is
            if ec in ec_to_reactions:
                # Check all reactions mapping to that EC number to make sure
                # we have seen that reaction before
                for rxnid in ec_to_reactions[ec]:
                    if rxnid in modeldata.reactions:
                        suggested_reactions.add(rxnid)

    # Remove reactions we already have
    suggested_reactions = suggested_reactions.difference(reactions2run)
    log_and_message(f"Gapfilling by EC number found {len(suggested_reactions)} new reactions", stderr=verbose)

    return suggested_reactions
