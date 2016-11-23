import os
import sys
import re

MODELSEED_DIR = os.environ['ModelSEEDDatabase']

def suggest_reactions_using_ec(roles, reactions, reactions2run, rf="SOLRDump/Reactions.tsv", verbose=False):
    """
    Identify a set of reactions that you should add to your model for growth based on the EC numbers
    that may be found in the role names.

    :param roles: A set of all roles to search for EC numbers.
    :type roles: set
    :param reactions: our reactions dictionary from parsing the model seed
    :type reactions: dict
    :param reactions2run: set of reactions that  we are going to run
    :type reactions2run: set
    :param rf: a reactions file from the SEED
    :type rf: str
    :param verbose: add additional output
    :type verbose: bool
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set
    """

    if not os.path.exists(os.path.join(MODELSEED_DIR, rf)):
        sys.stderr.write("FATAL: The reactions file {} does not exist from the directory {}.".format(rf, MODELSEED_DIR) +
                         " Please provide a path to that file\n")
        return set()

    # Read the ModelSEED reactions file
    ec_to_reactions = {}
    with open(os.path.join(MODELSEED_DIR, rf), "r") as rin:
        for l in rin:
            if l.startswith("#") or l.startswith("id"):
                # Ignore comment lines
                continue
            ll = l.strip().split("\t")
            rxnid = ll[0]
            ecs = ll[13]
            # EC number might be null
            if ecs == "null":
                continue
            # Multiple EC numbers can be assigned to a reaction
            for e in ecs.split(";"):
                if e not in ec_to_reactions:
                    ec_to_reactions[e] = set()
                ec_to_reactions[e].add(rxnid)

    # Find all EC numbers in the list of roles
    suggested_reactions = set()
    for role in roles:
        # Extract each EC number
        for ec in re.findall("[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+", role):
            # Check if we know what that EC number is
            if ec in ec_to_reactions:
                # Check all reactions mapping to that EC number to make sure
                # we have seen that reaction before
                for rxnid in ec_to_reactions[ec]:
                    if rxnid in reactions:
                        suggested_reactions.add(rxnid)

    if verbose:
        sys.stderr.write("Found " + str(len(suggested_reactions)) + " reactions\n")

    # Remove reactions we already have
    suggested_reactions = suggested_reactions.difference(reactions2run)

    if verbose:
        sys.stderr.write("Suggesting " + str(len(suggested_reactions)) + " reactions\n")

    return suggested_reactions
