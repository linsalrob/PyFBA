from PyFBA import log_and_message


def suggest_from_media(modeldata, reactions2run, media, verbose=False):
    """
    Identify a set of reactions that you should add to your model for growth based on the media compounds

    :param modeldata: the model seed object that includes compounds and reactions
    :type modeldata: PyFBA.model_seed.ModelData
    :param verbose: Print more output
    :type verbose: bool
    :param reactions2run: The reactions we are running
    :type reactions2run: set.
    :param media: A set of the compounds in the media
    :type media: set.
    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set[str]
    """

    # which compounds are in our media
    suggest = set()
    cpdaliases = {}
    cpdnames = {c.name: c for c in modeldata.compounds}

    for c in modeldata.compounds:
        if c.aliases:
            if 'Name' in c.aliases:
                cpdaliases.update({x.lower(): c for x in c.aliases['Name']})

    for m in media:
        # can we find it by name
        if m.name in cpdnames:
            rxns = set()
            for r in cpdnames[m.name].all_reactions():
                if r not in modeldata.reactions:
                    if not r.id.startswith('upsr'):
                        log_and_message(f"ERROR: {r} was not found in our reactions", stderr=verbose)
                    continue
                for c in modeldata.reactions[r].all_compounds():
                    if c.name == m.name and c.location == 'e':
                        rxns.add(r)
            log_and_message(f"Adding from media: For {m.name} added {len(rxns)} reactions", stderr=verbose)
            suggest.update(rxns)
        elif m.name in cpdaliases:
            log_and_message(f"Adding from media: Found {m.name} as an alias. Added {cpdaliases[m.name].name} and "
                            f"reactions {cpdaliases[m.name].all_reactions()}", stderr=verbose)
            rxns = set()
            for r in cpdaliases[m.name].all_reactions():
                if r not in modeldata.reactions:
                    log_and_message(f"ERROR: {r} was not found in our reactions", stderr=verbose)
                    continue
                for c in modeldata.reactions[r].all_compounds():
                    if c.name == m.name and c.location == 'e':
                        rxns.add(r)
            suggest.update(rxns)
        else:
            log_and_message(f"Compound {m.name} does not exist in the compound database", stderr=verbose)

    suggest = {r for r in suggest if r in modeldata.reactions and r not in reactions2run}

    return suggest
