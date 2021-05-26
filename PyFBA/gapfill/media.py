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
    cpdnames = {c.name: c for c in compounds}
    for m in media:
        # can we find it by name
        if m.name in cpdnames:
            rxns = cpdnames[m.name].all_reactions()
            suggest.update(rxns)
        else:
            if verbose:
                log_and_message(f"Compound {m.name} does not exist in the compound database", stderr=True)

    suggest = {r for r in suggest if r in reactions and r not in reactions2run}

    return suggest
