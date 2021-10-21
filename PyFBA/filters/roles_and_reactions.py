import sys
import re
import PyFBA
from PyFBA import log_and_message


def reactions_to_roles(reaction_set, organism_type=None, verbose=False):
    """
    Convert between reactions and roles using the model seed data

    For a set of reaction IDs return a hash of the reaction id and the
    roles in that reaction.

    :param reaction_set: A set of reaction IDs that we want to convert to roles
    :type reaction_set: set[str]
    :param verbose: print error reporting
    :type verbose: bool
    :return: a hash of reaction ids and set of the associated roles
    :rtype: dict[str, set(str)]
    """

    if not organism_type:
        log_and_message("WARNING: Organism type is not defined while gleaning the reactions. You should probably "
                        "specify e.g. Gram_Negative, Gram_Positive, etc. We used microbial",
                        stderr=verbose, loglevel="WARNING")
        organism_type = 'microbial'

    if isinstance(reaction_set, list):
        reaction_set = set(reaction_set)
    elif isinstance(reaction_set, str):
        reaction_set = {reaction_set}

    # key is complex and value is all reactions
    cmpxs = PyFBA.parse.model_seed.complexes(organism_type=organism_type, verbose=verbose)
    rct2cmpx = {}
    for c in cmpxs:
        for rxn in cmpxs[c]:
            if rxn not in rct2cmpx:
                rct2cmpx[rxn] = set()
            rct2cmpx[rxn].add(c)

    # key is role and value is all complexes
    roles = PyFBA.parse.model_seed.roles(organism_type=organism_type, verbose=verbose)
    cmpx2role = {}
    for r in roles:
        for c in roles[r]:
            if c not in cmpx2role:
                cmpx2role[c] = set()
            cmpx2role[c].add(r)

    roles = {}
    for r in reaction_set:
        if r not in rct2cmpx:
            if verbose:
                log_and_message(f"Converting reaction {r} to role: reaction not found in model_seed complexes for organism_type={organism_type}", stderr=True)
            continue
        roles[r] = set()
        for c in rct2cmpx[r]:
            if c not in cmpx2role:
                if verbose:
                    log_and_message(f"Complex {c} not found in the complexes", stderr=True)
                continue
            for rl in cmpx2role[c]:
                roles[r].add(rl)

    return roles


def roles_to_reactions(roles, organism_type=None, verbose=False):
    """
    Convert between roles and reactions using the model seed data

    For a set of roles return a hash where the key is the role and
    the value is the set of reactions that role is involved in.

    :param roles: A set of roles that we want to convert to reaction IDs
    :type roles: set
    :param organism_type: what type of organism is this?
    :type organism_type: str
    :param verbose: print error reporting
    :type verbose: bool
    :return: a hash of roles and set of the associated reaction ids
    :rtype: dict of set of str
    """

    if not organism_type:
        log_and_message("WARNING: Organism type is not defined while gleaning the reactions. You should probably "
                        "specify e.g. Gram_Negative, Gram_Positive, etc. We used microbial",
                        stderr=verbose, loglevel="WARNING")
        organism_type = 'microbial'
    if isinstance(roles, list):
        roles = set(roles)
    elif isinstance(roles, str):
        roles = {roles}

    # key is complex and value is all reactions
    cmpxs = PyFBA.parse.model_seed.complexes(organism_type=organism_type, verbose=verbose)
    # key is role and value is all complexes
    seedroles = PyFBA.parse.model_seed.roles(organism_type=organism_type, verbose=verbose)

    rcts = {}
    for r in roles:
        # check to see if it is a multifunctional role
        if '; ' in r or ' / ' in r or ' @ ' in r:
            log_and_message(f"{r} is a multifunctional role. You should separate the roles", stderr=verbose)
        if r not in seedroles:
            # I don't think we should report all missed roles as likely to be many
            # if verbose:
            #    log_and_message(f"Role {r} is not a role we understand. Skipped", stderr=verbose)
            continue

        rcts[r] = set()
        for c in seedroles[r]:
            if c not in cmpxs:
                if verbose:
                    # this occurs because there are reactions like cpx.1898 where we don't yet have a
                    # reaction for the complex
                    log_and_message(f"ERROR: {c} was not found in the complexes file, but is from a reaction",
                                    stderr=verbose)
                continue
            rcts[r].update(cmpxs[c])

    return rcts


def roles_to_ec_reactions(roles, organism_type=None, verbose=False):
    """
    For a list of roles, find those with EC numbers in them, parse out the EC number and see if any of our
    reactions have those EC numbers associated with them
    :param roles: a list of roles
    :type roles: list[str]
    :param organism_type: the type of organism, eg. Gram_Negative
    :type organism_type: str
    :param verbose: more output
    :type verbose: bool
    :return: a dict with roles as key and a set of reactions as value
    """

    rcts = {}
    seedreactions = PyFBA.parse.model_seed.reactions(organism_type=organism_type, verbose=verbose)
    ezs = {}
    for sr in seedreactions:
        if seedreactions[sr].ec_numbers:
            for ec in seedreactions[sr].ec_numbers:
                if ec not in ezs:
                    ezs[ec] = set()
                ezs[ec].add(seedreactions[sr].id)

    for r in roles:
        ecno2rctns = set()
        for ecno in re.findall(r'[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+', r):
            if ecno in ezs:
                ecno2rctns.update(ezs[ecno])
        if ecno2rctns:
            rcts[r] = ecno2rctns

    return rcts

if __name__ == '__main__':
    try:
        rt = sys.argv[1]
        rx = sys.argv[2]
    except IndexError:
        sys.exit(sys.argv[0] + " <type: either role or reaction> <id>")

    hs = set()
    if rt == 'role':
        hs = roles_to_reactions(rx)
    elif rt == 'reaction':
        hs = reactions_to_roles(rx)

    for h in hs:
        print(h + "\t" + ("\n" + h + "\t").join(hs[h]))

