import sys

from PyFBA import log_and_message


def suggest_from_media(compounds, reactions, reactions2run, media, verbose=False):
    """
    Identify a set of reactions that you should add to your model for growth based on the media compounds

    :param reactions: Our reactions dict object
    :type reactions: dict
    :param verbose: Print more output
    :type verbose: bool
    :param compounds: Our compounds dictionary
    :type compounds: dict
    :param reactions2run: The reactions we are running
    :type reactions2run: set.
    :param media: A set of the compounds in the media
    :type media: set.
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set
    """

    # which compounds are in our media
    suggest = set()
    cpdnames = {}
    cpdaliases = {}
    if type(compounds) == set:
        cpdnames = {c.name: c for c in compounds}
        for c in compounds:
            if c.aliases:
                if 'Name' in c.aliases:
                    cpdaliases.update({x.lower(): c for x in c.aliases['Name']})

    elif type(compounds) == dict:
        cpdnames = {compounds[c].name: compounds[c] for c in compounds}
        for c in compounds:
            if compounds[c].aliases:
                if 'Name' in compounds[c].aliases:
                    cpdaliases.update({x.lower(): compounds[c] for x in compounds[c].aliases['Name']})

    for m in media:
        # can we find it by name
        if m.name in cpdnames:
            rxns = cpdnames[m.name].all_reactions()
            if verbose:
                log_and_message(f"For {m.name} added {len(rxns)} reactions", stderr=True)
            suggest.update(rxns)
        elif m.name in cpdaliases:
            if verbose:
                log_and_message(f"Found {m.name} as an alias. Added {cpdaliases[m.name].name} and reactions " +
                                f"{cpdaliases[m.name].all_reactions()}", stderr=True)
            rxns = cpdaliases[m.name].all_reactions()
            suggest.update(rxns)
        else:
            if verbose:
                log_and_message(f"Compound {m.name} does not exist in the compound database", stderr=True)

    suggest = {r for r in suggest if r in reactions and r not in reactions2run}

    return suggest
