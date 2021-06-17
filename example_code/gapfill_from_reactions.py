"""
This code is designed to exemplify some of the gap-filling approaches. If you start with an ungapfilled set of
reactions, we iteratively try to build on the model until it is complete, and then we use the bisection code
to trim out reactions that are not necessary, to end up with the smallest set of reactions.

"""

import argparse
import copy
import os
import sys

import PyFBA
from PyFBA import log_and_message


def resolve_additional_reactions(ori_reactions, adnl_reactions, modeldata, mediaset, biomass_eqn):
    """
    Iteratively resolve additional reactions that are required.

    :param ori_reactions: original reactions
    :type ori_reactions: Set[str]
    :param modeldata: the model seed object that includes compounds and reactions
    :type modeldata: PyFBA.model_seed.ModelData
    :param adnl_reactions: a list of tuples of how the reactions were suggested, and the set of additional reactions
    :type adnl_reactions: list of tuple
    :param mediaset: our media object
    :type mediaset: Set[PyFBA.metabolism.CompoundWithLocation]
    :param biomass_eqn: our biomass object
    :type biomass_eqn: PyFBA.metabolism.Reaction
    :return: set of additional reactions from all of the added_reactions
    :rtype: set
    """

    reqd_additional = set()

    while adnl_reactions:
        ori = copy.copy(ori_reactions)
        ori.update(reqd_additional)
        (how, new) = adnl_reactions.pop()
        sys.stderr.write("Testing suggestions from " + how + "\n")
        # get all the other reactions we need to add
        for tple in adnl_reactions:
            ori.update(tple[1])
        new_essential = PyFBA.gapfill.minimize_additional_reactions(ori, new, modeldata, mediaset, biomass_eqn,
                                                                    verbose=True)
        for new_r in new_essential:
            modeldata.reactions[new_r].is_gapfilled = True
            modeldata.reactions[new_r].gapfill_method = how
        reqd_additional.update(new_essential)

    return reqd_additional


def minimal_set(rs, orr, arr, md, med, bme):
    # if the model grows then we want to find the minimal set of reactions
    # that we need to add for growth and call it good.

    addtns = resolve_additional_reactions(orr, arr, md, med, bme)
    for r in orr.union(addtns):
        if r not in rs:
            rs[r] = "UNKNOWN??"
        print("{}\t{}".format(r, rs[r]))
    sys.exit(0)


def read_media(mediafile, modeldata, verbose=False):
    """
    Read the media file and return a set of compounds
    :param mediafile: the media file to read
    :param modeldata: the model seed object that includes compounds and reactions
    :type modeldata: PyFBA.model_seed.ModelData
    :param verbose: more output
    :return: a set of media compounds
    :rtype: Set[PyFBA.metabolism.Compound]
    """

    if mediafile in PyFBA.parse.media_files():
        log_and_message(f"parsing media directly from {mediafile}", stderr=verbose)
        media = PyFBA.parse.pyfba_media(mediafile, modeldata)
    elif os.path.exists(mediafile):
        log_and_message(f"parsing media file {mediafile}", stderr=verbose)
        media = PyFBA.parse.read_media_file(mediafile)
    elif 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)):
        log_and_message(f"parsing media file {os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)}", stderr=verbose)
        media = PyFBA.parse.read_media_file(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile))
    else:
        log_and_message(f"Can't figure out how to parse media from {mediafile}", stderr=True, loglevel="CRITICAL")
        sys.exit(-1)
    return media


if __name__ == '__main__':
    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Import a list of reactions and then iterate through our gapfilling'
                                                 ' steps to see when we get growth')
    parser.add_argument('-r', help='reactions file', required=True)
    parser.add_argument('-m', help='media name', required=True)
    parser.add_argument('-c', help='close genomes reactions file')
    parser.add_argument('-g', help='other genera reactions file')
    parser.add_argument('-t', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.t not in orgtypes:
        sys.exit("Sorry, {} is not a valid organism type".format(args.t))

    # read the enzyme data
    # compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes(args.t)
    modeldata = PyFBA.parse.model_seed.parse_model_seed_data(args.t, verbose=args.v)

    reactionsource = {}

    reactions2run = set()
    with open(args.r, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue
            if "biomass" in l:
                if args.v:
                    sys.stderr.write("Biomass reaction was skipped from the list as it is auto-imported\n")
                continue
            r = l.strip()
            if r in modeldata.reactions:
                reactions2run.add(r)

    for r in reactions2run:
        reactionsource[r] = args.r

    media = read_media(args.m, modeldata, args.v)

    if args.t == 'gramnegative':
        biomass_eqtn = PyFBA.metabolism.biomass.biomass_equation('gramnegative')
    else:
        biomass_eqtn = PyFBA.metabolism.biomass.biomass_equation('standard')

    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn, verbose=args.v)
    sys.stderr.write("For the initial run we get growth of {} which is {}\n".format(value, growth))
    if growth:
        sys.exit("No need to gapfill!")

    added_reactions = []
    original_reactions = copy.copy(reactions2run)

    #############################################################################################
    #                                       Gapfilling                                          #
    #                                                                                           #
    #  We do this in the order:                                                                 #
    #     essential reactions: because you need to have these, but it is stronger evidence if   #
    #              your friends have it too!                                                    #
    #     media: because you should be importing everything in the media                        #
    #     linked_reactions: because they make sense!                                            #
    #     closely related organisms: because you should have roles your friends have            #
    #     subsystems: to complete things you already have                                       #
    #     orphans: to make sure everything is produced/consumed                                 #
    #     probability: because there are other reactions we can add                             #
    #     reactions with proteins: to make sure you can at least grow on the media              #
    #                                                                                           #
    #############################################################################################

    #############################################################################################
    #                                       ESSENTIAL PROTEINS                                  #
    #############################################################################################

    sys.stderr.write("Gap filling from ESSENTIAL PROTEINS\n")
    essential_reactions = PyFBA.gapfill.suggest_essential_reactions()
    # find only the new reactions
    essential_reactions.difference_update(reactions2run)
    added_reactions.append(("essential", essential_reactions))
    reactions2run.update(essential_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    log_and_message(f"After adding {len(essential_reactions)} ESSENTIAL reactions we get {value} (growth is {growth})",
                    stderr=args.v)
    for r in essential_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'essential_reactions'
    if growth:
        minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
        sys.exit(0)

    #############################################################################################
    #                                       Media import reactions                              #
    #############################################################################################

    log_and_message("Gap filling from MEDIA\n")
    media_reactions = PyFBA.gapfill.suggest_from_media(modeldata, reactions2run, media, verbose=args.v)
    added_reactions.append(("media", media_reactions))
    reactions2run.update(media_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    log_and_message(f"After adding {len(media_reactions)} MEDIA reactions we get {value} (growth is {growth})",
                    stderr=args.v)
    for r in media_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'media_reactions'

    if growth:
        minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
        sys.exit(0)

    #############################################################################################
    #                                       LINKED REACTIONS                                    #
    #############################################################################################

    sys.stderr.write("Gap filling from Linked Reactions")
    linked_reactions = PyFBA.gapfill.suggest_linked_reactions(modeldata, reactions2run)
    # find only the new reactions
    linked_reactions.difference_update(reactions2run)
    added_reactions.append(("linked", linked_reactions))
    reactions2run.update(linked_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    log_and_message(f"After adding {len(linked_reactions)} EC_NUMBER reactions we get {value} (growth is {growth})",
                    stderr=args.v)

    for r in linked_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'linked_reactions'

    if growth:
        minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
        sys.exit(0)

    #############################################################################################
    #                                        Other genomes and organisms                        #
    #############################################################################################

    sys.stderr.write("Gap filling from CLOSE GENOMES\n")
    close_reactions = set()
    if args.c:
        # add reactions from roles in close genomes
        close_reactions = PyFBA.gapfill.suggest_from_roles(args.c, modeldata.reactions, threshold=0, verbose=True)
        # find the new reactions
        close_reactions.difference_update(reactions2run)
        added_reactions.append(("close genomes ", close_reactions))
        reactions2run.update(close_reactions)
        status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
        sys.stderr.write("After adding {} reactions in {} we get {} (growth is {})\n\n".format(len(close_reactions),
                                                                                               args.c, value, growth))
        for r in close_reactions:
            if r not in reactionsource:
                reactionsource[r] = 'close_reactions'

        if growth:
            minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
            sys.exit(0)

    genus_reactions = set()
    if args.g:
        # add reactions from roles in similar genera
        genus_reactions = PyFBA.gapfill.suggest_from_roles(args.g, modeldata.reactions, threshold=0, verbose=True)
        # find the new reactions
        genus_reactions.difference_update(reactions2run)
        added_reactions.append(("other genera", genus_reactions))
        reactions2run.update(genus_reactions)
        status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
        sys.stderr.write("After adding {} reactions in {} we get {} (growth is {})\n\n".format(len(genus_reactions),
                                                                                               args.g, value, growth))
        if args.v:
            sys.stderr.write("REACTIONS FROM GENUS GENOMES: {}\n".format(genus_reactions))
        for r in genus_reactions:
            if r not in reactionsource:
                reactionsource[r] = 'genus_reactions'

        if growth:
            minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
            sys.exit(0)

    #############################################################################################
    #                                        Subsystems                                         #
    #############################################################################################

    sys.stderr.write("Gap filling from SUBSYSTEMS\n")
    subsystem_reactions = PyFBA.gapfill.suggest_reactions_from_subsystems(modeldata.reactions, reactions2run,
                                                                          threshold=0.5, verbose=args.v)
    added_reactions.append(("subsystems", subsystem_reactions))
    reactions2run.update(subsystem_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    sys.stderr.write("After adding {} SUBSYSTEM reactions we get {} (growth is {})\n\n".format(len(subsystem_reactions),
                                                                                               value, growth))
    for r in subsystem_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'subsystem_reactions'
    if growth:
        minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
        sys.exit(0)

    #############################################################################################
    #                                        Orphan compounds                                   #
    #############################################################################################

    sys.stderr.write("Gap filling from ORPHAN COMPOUNDS\n")
    orphan_compounds = PyFBA.gapfill.suggest_by_compound(modeldata, reactions2run, 1)
    added_reactions.append(("orphans", orphan_compounds))
    reactions2run.update(orphan_compounds)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    sys.stderr.write("After adding {} ORPHAN compounds we get {} (growth is {})\n\n".format(len(orphan_compounds),
                                                                                            value, growth))
    for r in orphan_compounds:
        if r not in reactionsource:
            reactionsource[r] = 'orphan_compounds'

    if growth:
        minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
        sys.exit(0)

    #############################################################################################
    #                                        Probability of inclusion                           #
    #############################################################################################

    sys.stderr.write("Gap filling from PROBABILITY\n")
    # use reactions wtih pLR or pRL > cutoff
    prob_reactions = PyFBA.gapfill.compound_probability(modeldata.reactions, reactions2run, 0, True, True)
    prob_reactions.difference_update(reactions2run)
    added_reactions.append(("probability", prob_reactions))
    reactions2run.update(prob_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    sys.stderr.write("After adding {} PROBABILITY reactions we get {} (growth is {})\n\n".format(len(prob_reactions),
                                                                                                 value, growth))

    for r in prob_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'probable_reactions'

        if growth:
            minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
            sys.exit(0)

    #############################################################################################
    #                       Reactions that [do or do not] map to proteins                       #
    #############################################################################################

    sys.stderr.write("Gap filling from ALL OTHER REACTIONS WITH PROTEINS\n")
    # propose other reactions that we have proteins for
    with_p_reactions = PyFBA.gapfill.suggest_reactions_with_proteins(modeldata.reactions, True)
    # find the new reactions
    with_p_reactions.difference_update(reactions2run)
    added_reactions.append(("With proteins", with_p_reactions))
    reactions2run.update(with_p_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    sys.stderr.write("After adding {} ALL WITH PROTEINS reactions ".format(len(with_p_reactions)) +
                     " we get {} (growth is {})\n\n".format(value, growth))

    for r in with_p_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'reactions_with_proteins'

    if growth:
        minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
        sys.exit(0)

    #############################################################################################
    #                       Reactions that do not map to proteins                               #
    #                                                                                           #
    #    This is all other reactions, and should really be avoided, but we include this         #
    #    in case the model does not grow by the time we get this far. If this doesn't work      #
    #    your model probably can not grow on the media that you have given it!                  #
    #                                                                                           #
    #############################################################################################

    sys.stderr.write("Gap filling from ALL OTHER REACTIONS WITHOUT PROTEINS\n")
    # propose other reactions that we have proteins for
    without_p_reactions = PyFBA.gapfill.suggest_reactions_without_proteins(modeldata.reactions, True)
    # we have to limit this to things we have compounds in our reaction list, or we will not be able to solve the
    # FBA (We may not be able to anyway)
    without_p_reactions = PyFBA.gapfill.limit_reactions_by_compound(modeldata.reactions, reactions2run,
                                                                    without_p_reactions)
    # find the new reactions
    without_p_reactions.difference_update(reactions2run)
    added_reactions.append(("Without proteins", without_p_reactions))
    reactions2run.update(without_p_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
    sys.stderr.write("After adding {} reactions without proteins ... last ditch ".format(len(without_p_reactions)) +
                     " we get {} (growth is {})\n\n".format(value, growth))
    if growth:
        sys.stderr.write("Congratulations, your model grew. You need to figure out which reactions we needed to add\n")
    else:
        sys.stderr.write("Your model does not grow regardless of how many reactions we add to it.\n")
        sys.stderr.write("Your media composition is probably lacking something essential (like C, N, S, P).\n")
        sys.stderr.write("Check that first\n")

    for r in without_p_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'reactions_without_proteins'

    if growth:
        minimal_set(reactionsource, original_reactions, added_reactions, modeldata, media, biomass_eqtn)
        sys.exit(0)

