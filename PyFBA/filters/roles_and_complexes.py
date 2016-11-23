import sys
import PyFBA


def roles_to_complexes(roles, verbose=False):
    """
    Convert between roles and complexes using the model seed data

    For a set of roles return a hash where the key is the role and
    the value is the set of  that role is involved in.

    For a set of roles return a hash containing two keys: complete and
    incomplete.
    "complete" points to a set of complexes where each complex has
    all of its roles present.
    "incomplete" points to a set of complexes where each complex has at least
    one role present but not all of its roles are present.

    :param roles: A set of roles that we want to convert to reaction IDs
    :type roles: set
    :param verbose: print error reporting
    :type verbose: bool
    :return: a hash of sets of complexes containing two keys:
             "complete" and "incomplete"
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
    # Build reverse of seedroles: key is complex and value is all roles
    seedcpxs = {}
    for r, cpxset in seedroles.items():
        for c in cpxset:
            # Check if role's complex is in our modelseed complex set
            if c not in cmpxs:
                if verbose:
                    # this occurs because there are reactions like cpx.1898 where we don't yet have a
                    # reaction for the complex
                    sys.stderr.write("ERROR: " + c + " was not found in the complexes file, but is from a reaction\n")
                continue
            # Record the complex to role mapping
            if c not in seedcpxs:
                seedcpxs[c] = set()
            seedcpxs[c].add(r)

    # Record which roles we have for our complexes
    mycpxs = {}
    for r in roles:
        # check to see if it is a multifunctional role
        if '; ' in r or ' / ' in r or ' @ ' in r:
            sys.stderr.write("It seems that {} is a multifunctional role. You should separate the roles\n".format(r))
        if r not in seedroles:
            if verbose:
                sys.stderr.write(r + " is not a role we understand. Skipped\n")
            continue

        for c in seedroles[r]:
            if c not in mycpxs:
                mycpxs[c] = set()
            mycpxs[c].add(r)

    # Determine which of our complexes are complete and incomplete
    ret_cpx = {"complete": set(), "incomplete": set()}
    for c, roleset in mycpxs.items():
        if c not in seedcpxs:
            continue
        which = "complete"
        for r in seedcpxs[c]:
            if r not in roleset:
                which = "incomplete"
                break
        ret_cpx[which].add(c)

    return ret_cpx
