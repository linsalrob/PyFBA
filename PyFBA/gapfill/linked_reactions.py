"""
Suggest reactions because they are "linked" to reactions that we are all ready going to run
"""
from PyFBA import log_and_message


def suggest_linked_reactions(modeldata, reactions_to_run, verbose=False):
    """
    Propose linked reactions
    :param modeldata: our model data structure
    :type: PyFBA.model_seed.ModelData
    :param reactions_to_run: our reactions to run
    :type reactions_to_run: set[str]
    :param verbose: more output
    :type verbose: bool
    :return: a set of reactions to add
    :rtype: set[str]
    """

    tosuggest = set()
    for r in reactions_to_run:
        if modeldata.reactions[r].linked_reaction:
            for lr in modeldata.reactions[r].linked_reaction.split(";"):
                tosuggest.add(lr)

    tosuggest.difference_update(reactions_to_run)
    log_and_message(f"Gapfill by linked reactions found {len(tosuggest)} new reactions", stderr=verbose)
    return tosuggest