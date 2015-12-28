import sys

import PyFBA


def suggest_from_roles(roles_file, reactions, threshold=0, verbose=False):
    """
    Identify a set of reactions that we should add based on a roles file.

    We assume that the roles file has the format: [role, probability] separated by tabs. We make no assumption about
    where you got that file, but you might, for example, look at the closely related organisms.

    :param threshold: the threshold for inclusion of the role based on the probability in the file (default = All roles)
    :type threshold: float
    :param roles_file: a file with a list of roles and their probabilities
    :type roles_file: str.
    :param reactions: The reactions dictionary from parsing the model seed
    :type reactions: dict.
    :param verbose: add additional output
    :type verbose: bool.
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set
    """

    role_suggestions = {}
    with open(roles_file, 'r') as rin:
        for l in rin:
            p = l.rstrip().split("\t")
            try:
                if float(p[1]) >= threshold:
                    for role in PyFBA.parse.rast.roles_of_function(p[0]):
                        role_suggestions[role] = p[1]
            except IndexError as e:
                sys.stderr.write("{} does not have enough columns\n".format(l.rstrip()))


    reaction_suggestions = set()
    ro2rx = PyFBA.filters.roles_to_reactions(set(role_suggestions.keys()))

    for rxnset in ro2rx.values():
        reaction_suggestions.update(rxnset)

    # limit this to only those reactions that we know about!
    reaction_suggestions = {x for x in reaction_suggestions if x in reactions}

    if verbose:
        sys.stderr.write("From {} we found {} role suggestions\n".format(roles_file, len(role_suggestions)))
        sys.stderr.write("Those roles converted to {} reactions\n".format(len(reaction_suggestions)))

    return reaction_suggestions
