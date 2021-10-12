"""
Run gapfilling and return a set of reactions when we grow
"""
import copy
import sys

import PyFBA
from PyFBA import log_and_message


def update_r2r(old, new, why, verbose=False):
    """
    Update the reactions to run and log the changes
    :param old: the initial reactions to run
    :param new: the new reactions to add
    :param why: the step we are at
    :param verbose: more output
    :return: a set of reactions to run
    :rtype: set[str]
    """
    before = len(old)
    old.update(new)
    msg = f"Before updating reactions from {why}: from {before} reactions to {len(old)} reactions"
    log_and_message(msg, stderr=verbose)
    return old


def run_eqn(why, md, r2r, med, bme, verbose=False):
    """
    Run the fba
    :param why: why are we doing this
    :param md: modeldata
    :param r2r: reactions to run
    :param med: media object
    :param bme: biomass equation
    :param verbose: more output
    :type verbose: bool
    :return: (value, growth)
    """

    status, value, growth = PyFBA.fba.run_fba(md, r2r, med, bme)
    log_and_message(f"FBA run {why} has a biomass flux value of {value} --> Growth: {growth}", stderr=verbose)
    return value, growth


def gapfill(reactions, model_data, growth_media, biomass_eqtn, close, genome_type, r2exclude=None, verbose=False):
    """
    Gapfill a set of reactions and return a tuple of [new reactions that grow, [reason, list of reactions]].

    We do not bisect these reactions, and you should do so

    :param reactions: set of initial reactions
    :param model_data: the model data
    :param growth_media: the growth media object
    :param biomass_eqtn: the biomass equation
    :param close: a list of closely related genomes
    :param genome_type: the genome type
    :param verbose: more output
    :param r2exclude: a set of reactions to exclude (optional)
    :return: a dict of the reactions and why they are there!
    :rtype: dict[str, str]
    """

    if r2exclude is None:
        r2exclude = {}
    added_reactions = []
    original_reactions_to_run = copy.deepcopy(reactions)

    val, growth = run_eqn(f"Initial test to make sure we don't grow!", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return {r:"initial" for r in reactions}

    #############################################################################################
    #                                       Gapfilling                                          #
    #                                                                                           #
    #  We do this in the order:                                                                 #
    #     essential reactions: because you need to have these, but it is stronger evidence if   #
    #              your friends have it too!                                                    #
    #     linked reactions: connections to other reactions in our set                           #
    #     media: because you should be importing everything in the media                        #
    #     closely related organisms: because you should have roles your friends have            #
    #     subsystems: to complete things you already have                                       #
    #     orphans: to make sure everything is produced/consumed                                 #
    #                                                                                           #
    #     At the moment we don't run these approaches, but easy to add back again!              #
    #     probability: because there are other reactions we can add                             #
    #     reactions with proteins: to make sure you can at least grow on the media              #
    #                                                                                           #
    #############################################################################################

    #############################################################################################
    #                                       ESSENTIAL PROTEINS                                  #
    #############################################################################################

    step = "Essential Reactions"
    essential_reactions = PyFBA.gapfill.suggest_essential_reactions()
    essential_reactions.difference_update(r2exclude)
    for r in essential_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("essential", essential_reactions))
    reactions = update_r2r(reactions, essential_reactions, "ESSENTIAL REACTIONS")

    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                                       LINKED REACTIONS                                    #
    #############################################################################################

    step = "Linked Reactions"
    linked_reactions = PyFBA.gapfill.suggest_linked_reactions(model_data, reactions)
    linked_reactions.difference_update(r2exclude)
    for r in linked_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("linked_reactions", linked_reactions))
    reactions = update_r2r(reactions, linked_reactions, "LINKED REACTIONS")

    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)

    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                                       Media import reactions                              #
    #############################################################################################

    step = "media"
    media_reactions = PyFBA.gapfill.suggest_from_media(model_data, reactions, growth_media, verbose)
    media_reactions.difference_update(r2exclude)
    for r in media_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("media", media_reactions))
    reactions = update_r2r(reactions, media_reactions, "MEDIA REACTIONS")

    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                                        Other genomes and organisms                        #
    #############################################################################################

    if close:
        for close_genome in close:
            step = f"Close genome: {close_genome}"
            # add reactions from roles in close genomes
            close_reactions = PyFBA.gapfill.suggest_from_roles(close_genome, model_data.reactions, threshold=0,
                                                               verbose=verbose)
            # find the new reactions
            close_reactions.difference_update(reactions)
            close_reactions.difference_update(r2exclude)
            for r in close_reactions:
                model_data.reactions[r].reset_bounds()
            added_reactions.append((step, close_reactions))
            reactions = update_r2r(reactions, close_reactions, "CLOSE REACTIONS")
            val, growth = run_eqn(f"Test growth after {step}", model_data,
                                  r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
            if growth:
                return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data,
                                                        growth_media,
                                                        biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                                        Subsystems                                         #
    #############################################################################################

    step = "Subsytems"
    subsystem_reactions = PyFBA.gapfill.suggest_reactions_from_subsystems(model_data.reactions, reactions,
                                                                          organism_type=genome_type,
                                                                          threshold=0.5, verbose=verbose)
    subsystem_reactions.difference_update(r2exclude)
    for r in subsystem_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("subsystems", subsystem_reactions))
    reactions = update_r2r(reactions, subsystem_reactions, "SUBSYSTEMS REACTIONS")
    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                                        Orphan compounds                                   #
    #############################################################################################

    step = "Orphan compounds"
    orphan_reactions = PyFBA.gapfill.suggest_by_compound(model_data, reactions, 1)
    orphan_reactions.difference_update(r2exclude)
    for r in orphan_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("orphans", orphan_reactions))
    reactions = update_r2r(reactions, orphan_reactions, "ORPHAN REACTIONS")
    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                                        Probability of inclusion                           #
    #############################################################################################

    step = "Probability"
    # use reactions wtih pLR or pRL > cutoff
    prob_reactions = PyFBA.gapfill.compound_probability(model_data.reactions, reactions, 0, True, True)
    prob_reactions.difference_update(reactions)
    prob_reactions.difference_update(r2exclude)
    for r in prob_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("probability", prob_reactions))
    reactions = update_r2r(reactions, prob_reactions, "PROBABILITY OF REACTIONS")

    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                       Reactions that map to proteins                       #
    #############################################################################################

    step = "Gap filling from ALL OTHER REACTIONS WITH PROTEINS"
    # propose other reactions that we have proteins for
    with_p_reactions = PyFBA.gapfill.suggest_reactions_with_proteins(model_data.reactions, True)
    # find the new reactions
    with_p_reactions.difference_update(reactions)
    with_p_reactions.difference_update(r2exclude)
    for r in with_p_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("With proteins", with_p_reactions))
    reactions = update_r2r(reactions, with_p_reactions, "REACTIONS WITH PROTEINS")

    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    #############################################################################################
    #                       Reactions that do not map to proteins                               #
    #                                                                                           #
    #    This is all other reactions, and should really be avoided, but we include this         #
    #    in case the model does not grow by the time we get this far. If this doesn't work      #
    #    your model probably can not grow on the media that you have given it!                  #
    #                                                                                           #
    #############################################################################################

    log_and_message("Gap filling from ALL OTHER REACTIONS WITHOUT PROTEINS", stderr=verbose)
    # propose other reactions that we have proteins for
    without_p_reactions = PyFBA.gapfill.suggest_reactions_without_proteins(model_data.reactions, True)
    # we have to limit this to things we have compounds in our reaction list, or we will not be able to solve the
    # FBA (We may not be able to anyway)
    without_p_reactions = PyFBA.gapfill.limit_reactions_by_compound(model_data.reactions, reactions,
                                                                    without_p_reactions)
    # find the new reactions
    without_p_reactions.difference_update(reactions)
    without_p_reactions.difference_update(r2exclude)
    for r in without_p_reactions:
        model_data.reactions[r].reset_bounds()
    added_reactions.append(("Without proteins", without_p_reactions))
    reactions = update_r2r(reactions, without_p_reactions, "REACTIONS WITHOUT PROTEINS")

    val, growth = run_eqn(f"Test growth after {step}", model_data,
                          r2r=reactions, bme=biomass_eqtn, med=growth_media, verbose=verbose)
    if growth:
        return PyFBA.gapfill.minimize_reactions(original_reactions_to_run, added_reactions, model_data, growth_media,
                                                biomass_eqtn, verbose=verbose)

    log_and_message(f"ERROR: WE COULD NOT GAPFILL TO GET GROWTH, EVEN WITH ALL REACTIONS", stderr=True, loglevel="WARN")
    sys.exit(0)
