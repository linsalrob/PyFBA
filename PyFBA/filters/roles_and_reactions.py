import sys

import PyFBA


def reactions_to_roles(reaction_set, verbose=False):
    """
    Convert between reactions and roles using the model seed data

    For a set of reaction IDs return a hash of the reaction id and the
    roles in that reaction.

    :param reaction_set: A set of reaction IDs that we want to convert to roles
    :type reaction_set: set
    :param verbose: print error reporting
    :type verbose: bool
    :return: a hash of reaction ids and set of the associated roles
    :rtype: dict of set of str
    """

    if isinstance(reaction_set, list):
        reaction_set = set(reaction_set)
    elif isinstance(reaction_set, str):
        reaction_set = {reaction_set}

    # key is complex and value is all reactions
    cmpxs = PyFBA.parse.model_seed.complexes()
    # key is role and value is all complexes
    roles = PyFBA.parse.model_seed.roles()
    
    rct2cmpx = {}
    for c in cmpxs:
        for r in cmpxs[c]:
            if r not in rct2cmpx:
                rct2cmpx[r] = set()
            rct2cmpx[r].add(c)

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
                sys.stderr.write("ERROR " + r + " not found\n")
            continue
        roles[r] = set()
        for c in rct2cmpx[r]:
            if c not in cmpx2role:
                if verbose:
                    sys.stderr.write("Complex " + c + " not found in the complexes\n")
                continue
            for rl in cmpx2role[c]:
                roles[r].add(rl)

    return roles


def roles_to_reactions(roles, verbose=False):
    """
    Convert between roles and reactions using the model seed data

    For a set of roles return a hash where the key is the role and 
    the value is the set of reactions that role is involved in.

    :param roles: A set of roles that we want to convert to reaction IDs
    :type roles: set
    :param verbose: print error reporting
    :type verbose: bool
    :return: a hash of roles and set of the associated reaction ids
    :rtype: dict of set of str
    """

    if isinstance(roles, list):
        roles = set(roles)
    elif isinstance(roles, str):
        roles = {roles}

    # key is complex and value is all reactions
    cmpxs = PyFBA.parse.model_seed.complexes()
    # key is role and value is all complexes
    seedroles = PyFBA.parse.model_seed.roles()

    rcts = {}
    for r in roles:
        if r not in seedroles:
            if verbose:
                sys.stderr.write(r + " is not a role we understand. Skipped\n")
            continue

        rcts[r] = set()
        for c in seedroles[r]:
            if c not in cmpxs:
                if verbose:
                    # this occurs because there are reactions like cpx.1896 where we don't yet have a
                    # reaction for the complex
                    sys.stderr.write("ERROR: " + c + " was not found in the complexes file, but is from a reaction\n")
                continue
            for rc in cmpxs[c]:
                rcts[r].add(rc)

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

