import argparse
import copy
import sys

import PyFBA

"""
This code is designed to exemplify some of the gap-filling approaches. If you start with an ungapfilled set of
reactions, we iteratively try to build on the model until it is complete, and then we use the bisection code
to trim out reactions that are not necessary, to end up with the smallest set of reactions.

We will handle multiple positive files and negative files and try and resolve in the most meaningful way (hopefully!)

"""


def resolve_additional_reactions(ori_reactions, adnl_reactions, cpds, rcts, growth_media, no_growth_media, biomass_eqn,
                                 minimum_tp=0, minimum_accuracy=0.5, verbose=False):
    """
    Iteratively resolve additional reactions that are required.

    :param cpds: Our compounds dictionary object
    :type cpds: dict
    :param ori_reactions: the set of original reactions that form the base of the model
    :type ori_reactions: set
    :param adnl_reactions: a list of tuples of how the reactions were suggested, and the set of additional reactions
    :type adnl_reactions: list of tuple
    :param rcts: our reactions object
    :type rcts: dict
    :param growth_media: A list of media where we should grow
    :type growth_media: list of Set of Compounds
    :param no_growth_media: A list of media where we should not grow
    :type no_growth_media: list of Set of Compounds
    :param biomass_eqn: our biomass object
    :type biomass_eqn: metabolism.Reaction
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
        new_essential = PyFBA.gapfill.minimize_by_accuracy(ori, new, cpds, rcts, growth_media, no_growth_media,
                                                           biomass_eqn, minimum_tp, minimum_accuracy, verbose=verbose)
        for new_r in new_essential:
            reactions[new_r].is_gapfilled = True
            reactions[new_r].gapfill_method = how
        reqd_additional.update(new_essential)

    return reqd_additional


if __name__ == '__main__':
    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Import a list of reactions and then iterate through our gapfilling'
                                                 ' steps to see when we get growth. You can specify multiple -p & -n')
    parser.add_argument('-r', help='reactions file', required=True)
    parser.add_argument('-p', help='media file(s) on which the organism can grow', required=True, action='append')
    parser.add_argument('-n', help='media file(s) on which the organism can NOT grow', default=[], action='append')
    parser.add_argument('-f', help='fraction of growth conditions on which we want growth', default=0.8, type=float)
    parser.add_argument('-c', help='close genomes reactions file')
    parser.add_argument('-g', help='other genera reactions file')
    parser.add_argument('-t', help='organism type for the model (currently allowed are {}). Default=gramnegative'.format(
        orgtypes), default='gramnegative')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.t not in orgtypes:
        sys.exit("Sorry, {} is not a valid organism type".format(args.t))

    # read the enzyme data
    compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes(args.t)

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
            if r in reactions:
                reactions2run.add(r)

    for r in reactions2run:
        reactionsource[r] = args.r

    growth_media = [PyFBA.parse.read_media_file(m) for m in args.p]
    no_growth_media = [PyFBA.parse.read_media_file(m) for m in args.n]

    biomass_eqtn = PyFBA.metabolism.biomass.biomass_equation('gramnegative')

    min_growth_conditions = args.f
    if args.f < 1:
        min_growth_conditions = 1.0 * args.f * len(growth_media)

    sys.stderr.write("Looking for growth on at least {} media\n".format(min_growth_conditions))

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions, reactions2run, biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("For the initial run we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        sys.exit("{} is bigger than {} : no need to gapfill!".format(pr['tp'], min_growth_conditions))

    added_reactions = []
    original_reactions = copy.copy(reactions2run)

    #############################################################################################
    #                                       Gapfilling                                          #
    #                                                                                           #
    #  We do this in the order:                                                                 #
    #     media: because you should be importing everything in the media                        #
    #     closely related organisms: because you should have roles your friends have            #
    #     essential reactions: because you need to have these, but it is stronger evidence if   #
    #              your friends have it too!                                                    #
    #     subsystems: to complete things you already have                                       #
    #     orphans: to make sure everything is produced/consumed                                 #
    #     probability: because there are other reactions we can add                             #
    #     reactions with proteins: to make sure you can at least grow on the media              #
    #                                                                                           #
    #############################################################################################

    #############################################################################################
    #                                       Media import reactions                              #
    #############################################################################################

    sys.stderr.write("Gap filling from MEDIA\n")
    media_reactions = set()
    for media in growth_media:
        media_reactions.update(PyFBA.gapfill.suggest_from_media(compounds, reactions, reactions2run, media, verbose=args.v))
    added_reactions.append(("media", media_reactions))
    reactions2run.update(media_reactions)

    for r in media_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'media_reactions'

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions, reactions2run,
                                                  biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("After adding media we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=args.v)
        for r in original_reactions.union(additions):
            if r not in reactionsource:
                reactionsource[r] = "UNKNOWN??"
            print("{}\t{}".format(r, reactionsource[r]))
        sys.exit(0)

    #############################################################################################
    #                                        Other genomes and organisms                        #
    #############################################################################################

    sys.stderr.write("Gap filling from CLOSE GENOMES\n")
    close_reactions = set()
    if args.c:
        # add reactions from roles in close genomes
        close_reactions = PyFBA.gapfill.suggest_from_roles(args.c, reactions, threshold=0, verbose=True)
        # find the new reactions
        close_reactions.difference_update(reactions2run)
        added_reactions.append(("close genomes ", close_reactions))
        reactions2run.update(close_reactions)
        for r in close_reactions:
            if r not in reactionsource:
                reactionsource[r] = 'close_reactions'

        pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions, reactions2run,
                                                      biomass_eqtn)
        acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

        sys.stderr.write("After adding reactions from {} we get: ".format(args.c))
        for t in ['tp', 'fp', 'tn', 'fn']:
            sys.stderr.write("{}: {} ".format(t, pr[t]))
        sys.stderr.write("Accuracy: {}\n".format(acc))

        if pr['tp'] > min_growth_conditions:
            additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                     growth_media, no_growth_media, biomass_eqtn,
                                                     minimum_tp=min_growth_conditions, verbose=args.v)
            for r in original_reactions.union(additions):
                if r not in reactionsource:
                    reactionsource[r] = "UNKNOWN??"
                print("{}\t{}".format(r, reactionsource[r]))
            sys.exit(0)

    genus_reactions = set()
    if args.g:
        # add reactions from roles in similar genera
        genus_reactions = PyFBA.gapfill.suggest_from_roles(args.g, reactions, threshold=0, verbose=True)
        genus_reactions.difference_update(reactions2run)
        added_reactions.append(("other genera", genus_reactions))
        reactions2run.update(genus_reactions)

        for r in genus_reactions:
            if r not in reactionsource:
                reactionsource[r] = 'genus_reactions'

        pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions, reactions2run,
                                                      biomass_eqtn)
        acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

        sys.stderr.write("After adding reactions from {} we get: ".format(args.g))
        for t in ['tp', 'fp', 'tn', 'fn']:
            sys.stderr.write("{}: {} ".format(t, pr[t]))
        sys.stderr.write("Accuracy: {}\n".format(acc))

        if pr['tp'] > min_growth_conditions:
            additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                     growth_media, no_growth_media, biomass_eqtn,
                                                     minimum_tp=min_growth_conditions, verbose=args.v)
            for r in original_reactions.union(additions):
                if r not in reactionsource:
                    reactionsource[r] = "UNKNOWN??"
                print("{}\t{}".format(r, reactionsource[r]))
            sys.exit(0)

    #############################################################################################
    #                                       ESSENTIAL PROTEINS                                  #
    #############################################################################################

    sys.stderr.write("Gap filling from ESSENTIAL PROTEINS\n")
    essential_reactions = PyFBA.gapfill.suggest_essential_reactions()
    # find only the new reactions
    essential_reactions.difference_update(reactions2run)
    added_reactions.append(("essential", essential_reactions))
    reactions2run.update(essential_reactions)

    for r in essential_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'essential_reactions'

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions,
                                                  reactions2run,
                                                  biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("After adding essential reactions we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=args.v)
        for r in original_reactions.union(additions):
            if r not in reactionsource:
                reactionsource[r] = "UNKNOWN??"
            print("{}\t{}".format(r, reactionsource[r]))
        sys.exit(0)

    #############################################################################################
    #                                        Subsystems                                         #
    #############################################################################################

    sys.stderr.write("Gap filling from SUBSYSTEMS\n")
    subsystem_reactions = PyFBA.gapfill.suggest_reactions_from_subsystems(reactions, reactions2run, threshold=0.5, verbose=args.v)
    added_reactions.append(("subsystems", subsystem_reactions))
    reactions2run.update(subsystem_reactions)
    for r in subsystem_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'subsystem_reactions'

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions,
                                                  reactions2run,
                                                  biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("After adding subystems we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=args.v)
        for r in original_reactions.union(additions):
            if r not in reactionsource:
                reactionsource[r] = "UNKNOWN??"
            print("{}\t{}".format(r, reactionsource[r]))
        sys.exit(0)

    #############################################################################################
    #                                        Orphan compounds                                   #
    #############################################################################################

    sys.stderr.write("Gap filling from ORPHAN COMPOUNDS\n")
    orphan_compounds = PyFBA.gapfill.suggest_by_compound(compounds, reactions, reactions2run, 1)
    added_reactions.append(("orphans", orphan_compounds))
    reactions2run.update(orphan_compounds)
    for r in orphan_compounds:
        if r not in reactionsource:
            reactionsource[r] = 'orphan_compounds'

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions,
                                                  reactions2run,
                                                  biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("After adding orphans we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=args.v)
        for r in original_reactions.union(additions):
            if r not in reactionsource:
                reactionsource[r] = "UNKNOWN??"
            print("{}\t{}".format(r, reactionsource[r]))
        sys.exit(0)


    #############################################################################################
    #                                        Probability of inclusion                           #
    #############################################################################################

    sys.stderr.write("Gap filling from PROBABILITY\n")
    # use reactions wtih pLR or pRL > cutoff
    prob_reactions = PyFBA.gapfill.compound_probability(reactions, reactions2run, 0, True, True)
    prob_reactions.difference_update(reactions2run)
    added_reactions.append(("probability", prob_reactions))
    reactions2run.update(prob_reactions)

    for r in prob_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'probable_reactions'

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions,
                                                  reactions2run,
                                                  biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("After adding probability we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=args.v)
        for r in original_reactions.union(additions):
            if r not in reactionsource:
                reactionsource[r] = "UNKNOWN??"
            print("{}\t{}".format(r, reactionsource[r]))
        sys.exit(0)

    #############################################################################################
    #                       Reactions that map to proteins                       #
    #############################################################################################

    sys.stderr.write("Gap filling from ALL OTHER REACTIONS WITH PROTEINS\n")
    # propose other reactions that we have proteins for
    with_p_reactions = PyFBA.gapfill.suggest_reactions_with_proteins(reactions, True)
    # find the new reactions
    with_p_reactions.difference_update(reactions2run)
    added_reactions.append(("With proteins", with_p_reactions))
    reactions2run.update(with_p_reactions)
    for r in with_p_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'reactions_with_proteins'

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions,
                                                  reactions2run,
                                                  biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("After adding reactions with proteins we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=args.v)
        for r in original_reactions.union(additions):
            if r not in reactionsource:
                reactionsource[r] = "UNKNOWN??"
            print("{}\t{}".format(r, reactionsource[r]))
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
    without_p_reactions = PyFBA.gapfill.suggest_reactions_without_proteins(reactions, True)
    # we have to limit this to things we have compounds in our reaction list, or we will not be able to solve the
    # FBA (We may not be able to anyway)
    without_p_reactions = PyFBA.gapfill.limit_reactions_by_compound(reactions, reactions2run, without_p_reactions)
    # find the new reactions
    without_p_reactions.difference_update(reactions2run)
    added_reactions.append(("Without proteins", without_p_reactions))
    reactions2run.update(without_p_reactions)

    for r in without_p_reactions:
        if r not in reactionsource:
            reactionsource[r] = 'reactions_without_proteins'

    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, compounds, reactions,
                                                  reactions2run, biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    sys.stderr.write("After adding all reactions we get: ")
    for t in ['tp', 'fp', 'tn', 'fn']:
        sys.stderr.write("{}: {} ".format(t, pr[t]))
    sys.stderr.write("Accuracy: {}\n".format(acc))

    if pr['tp'] > min_growth_conditions:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=args.v)
        for r in original_reactions.union(additions):
            if r not in reactionsource:
                reactionsource[r] = "UNKNOWN??"
            print("{}\t{}".format(r, reactionsource[r]))
        sys.exit(0)
    else:
        sys.stderr.write("Sorry, regardless of how many reactions, we can never get more than {}".format(pr['tp']))
        sys.stderr.write("true positives.\n")
        sys.exit(-1)