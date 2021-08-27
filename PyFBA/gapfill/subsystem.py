try:
    from importlib.resources import open_text
except ImportError:
    # this is for python<3.7
    from importlib_resources import open_text

import PyFBA
from PyFBA import log_and_message


def suggest_reactions_from_subsystems(reactions, reactions2run, organism_type=None,
                                      ssfile="SS_functions.txt", threshold=0, verbose=False):
    """
    Identify a set of reactions that you should add to your model for growth based on the subsystems that are present
    in your model and their coverage.

    Read roles and subsystems from the subsystems file (which has role, subsystem, classification 1, classification 2)
    and make suggestions for missing reactions based on the subsystems that only have partial reaction coverage.

    :param organism_type: The type of the organism (gram -ve etc)
    :type organism_type: str
    :param threshold: The minimum fraction of the genes that are already in the subsystem for it to be added (default=0)
    :type threshold: float
    :param reactions: our reactions dictionary from parsing the model seed
    :type reactions: dict
    :param reactions2run: set of reactions that  we are going to run
    :type reactions2run: set
    :param ssfile: a subsystem file (really the output of dump_functions.pl on the seed machines)
    :type ssfile: str
    :param verbose: add additional output
    :type verbose: bool
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set
    """


    # read the ss file
    subsys_to_roles = {}
    roles_to_subsys = {}
    with open_text("PyFBA.Biochemistry.SEED.Subsystems", ssfile) as sin:
        for l in sin:
            if l.startswith('#'):
                continue
            p = l.rstrip().split("\t")
            if len(p) < 2:
                log_and_message(f"Too few columns in subsystem file at line: {l.strip()}", stderr=verbose)
                continue
            if p[1] not in subsys_to_roles:
                subsys_to_roles[p[1]] = set()
            for role in PyFBA.parse.roles_of_function(p[0]):
                if role not in roles_to_subsys:
                    roles_to_subsys[role] = set()
                subsys_to_roles[p[1]].add(role)
                roles_to_subsys[role].add(p[1])

    # now convert our reaction ids in reactions2run into roles
    # we have a hash with keys = reactions and values = set of roles
    reacts = PyFBA.filters.reactions_to_roles(reactions2run, organism_type, verbose=verbose)

    # foreach subsystem we need to know the fraction of roles that are present
    # this is complicated by multifunctional enzymes, as if one function is present they all should be
    # but for the moment (??) we are going to assume that each peg has the multi-functional annotation
    ss_present = {}
    ss_roles = {}
    for r in reacts:
        for rl in reacts[r]:
            if rl in roles_to_subsys:
                for s in roles_to_subsys[rl]:
                    if s not in ss_present:
                        ss_present[s] = set()
                        ss_roles[s] = set()
                    ss_present[s].add(rl)
                    ss_roles[s].add(r)

    ss_fraction = {}
    for s in ss_present:
        ss_fraction[s] = 1.0 * len(ss_present[s]) / len(subsys_to_roles[s])

    if verbose:
        for s in ss_roles:
            print("{}\t{}\t{}".format(s, ss_fraction[s], ss_roles[s], "; ".join(ss_present)))

    # now we can suggest the roles that should be added to complete subsystems.
    suggested_ss = set()
    for s in ss_fraction:
        if ss_fraction[s] >= threshold:
            suggested_ss.add(s)

    log_and_message(f"Subsystems suggesting {len(suggested_ss)} subsystems", stderr=verbose)

    # suggested_ss = {s for s, f in ss_fraction.items() if f>0}
    suggested_roles = set()
    for s in suggested_ss:
        for r in subsys_to_roles[s]:
            if r not in reactions2run:
                suggested_roles.add(r)

    log_and_message(f"Subsystems suggesting {len(suggested_roles)} roles", stderr=verbose)

    # finally, convert the roles to reactions
    new_reactions = PyFBA.filters.roles_to_reactions(suggested_roles, organism_type=organism_type, verbose=verbose)

    log_and_message(f"Subsystems suggesting {len(new_reactions)} total reactions", stderr=verbose)

    suggested_reactions = set()
    for rl in new_reactions:
        suggested_reactions.update(new_reactions[rl])

    log_and_message(f"Subsystems suggested {len(suggested_reactions)} unique reactions", stderr=verbose)

    suggested_reactions = {r for r in suggested_reactions if r in reactions and r not in reactions2run}

    return suggested_reactions
