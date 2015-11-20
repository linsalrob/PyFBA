import argparse
import copy
import sys

import fba
import gapfill
import parse
from metabolism import biomass
from parse import model_seed

__author__ = 'Rob Edwards'

"""
This code is designed to exmplify some of the gap-filling approaches. If you start with an ungapfilled set of
reactions, we iteratively try to build on the model until it is complete, and then we use the bisection code
to trim out reactions that are not necessary, to end up with the smallest set of reactions.

"""

def resolve_additional_reactions(ori_reactions, adnl_reactions, cpds, rcts, mediaset, biomass_eqn):
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
    :param mediaset: our media object
    :type mediaset: set
    :param biomass_eqn: our biomass object
    :type biomass_eqn: metabolism.Reaction
    :return: set of additional reactions from all of the added_reactions
    :rtype: set
    '''
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
        new_essential = gapfill.minimize_additional_reactions(ori, new, cpds, rcts, mediaset, biomass_eqn,
                                                              verbose=True)
        for new_r in new_essential:
            reactions[new_r].is_gapfilled = True
            reactions[new_r].gapfill_method = how
        reqd_additional.update(new_essential)

    return reqd_additional


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Import a list of reactions and then iterate through our gapfilling'
                                                 ' steps to see when we get growth')
    parser.add_argument('-r', help='reactions file', required=True)
    parser.add_argument('-m', help='media file', required=True)
    parser.add_argument('-c', help='close genomes reactions file')
    parser.add_argument('-g', help='other genera reactions file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    # read the enzyme data
    compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes('gramnegative')

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

    media = parse.read_media_file(args.m)
    biomass_eqtn = biomass.biomass_equation('gramnegative')

    status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass_eqtn, verbose=args.v)
    print("Initial run has " + str(value) + " --> Growth: " + str(growth))
    if growth:
        sys.exit("No need to gapfill!")

    added_reactions = []
    original_reactions = copy.copy(reactions2run)

    # get the suggestions
    essential_reactions = gapfill.suggest_essential_reactions(True)
    # find only the new reactions
    essential_reactions.difference_update(reactions2run)
    added_reactions.append(("essential", essential_reactions))
    reactions2run.update(essential_reactions)
    status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass_eqtn)
    print("After adding predefined reactions we get " + str(value) + " (growth is " + str(growth) + ")\n\n")

    # if this grows then we want to find the minimal set of reactions
    # that we need to add for growth and call it good.
    if growth:
        additional = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                  media, biomass_eqtn)
        print('reactions' + " : " + str(original_reactions.union(additional)))
        sys.exit(0)

    close_reactions = set()
    if args.c:
        # add reactions from roles in close genomes
        close_reactions = gapfill.suggest_from_roles(args.c, reactions, True)
        # find the new reactions
        close_reactions.difference_update(reactions2run)
        added_reactions.append(("close genomes ", close_reactions))
        reactions2run.update(close_reactions)

        sys.stderr.write("Reactions to run has " + str(len(reactions2run)) + " reactions\n")
        status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass_eqtn)
        sys.stderr.write("After adding reactions in " + args.c + " we get " +
                         str(value) + " (growth is " + str(growth) + ")\n\n")

        # if this grows then we want to find the minimal set of reactions
        # that we need to add for growth and call it good.
        if growth:
            additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                     media, biomass_eqtn)
            # print("Additional reactions required: " + str(additions) + "\n")
            print('reactions' + " : " + str(original_reactions.union(additions)))
            sys.exit(0)

    genus_reactions = set()
    if args.g:
        # add reactions from roles in similar genera
        genus_reactions = gapfill.suggest_from_roles(args.g, reactions, True)
        # find the new reactions
        genus_reactions.difference_update(reactions2run)
        added_reactions.append(("other genera", genus_reactions))
        reactions2run.update(genus_reactions)

        sys.stderr.write("Reactions to run has " + str(len(reactions2run)) + " reactions\n")
        status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass_eqtn)
        sys.stderr.write("After adding reactions in " + args.g + " we get " +
                         str(value) + " (growth is " + str(growth) + ")\n\n")

        # if this grows then we want to find the minimal set of reactions
        # that we need to add for growth and call it good.
        if growth:
            additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                     media, biomass_eqtn)
            # print("Additional reactions required: " + str(additions) + "\n")
            print('reactions' + " : " + str(original_reactions.union(additions)))
            sys.exit(0)

    # use reactions wtih pLR or pRL > cutoff
    prob_reactions = gapfill.compound_probability(reactions, reactions2run, 0, True, True)
    prob_reactions.difference_update(reactions2run)
    added_reactions.append(("probability", prob_reactions))
    reactions2run.update(prob_reactions)
    status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass_eqtn)
    sys.stderr.write("After adding by probability we get " +
                     str(value) + " (growth is " + str(growth) + ")\n\n")
    # if this grows then we want to find the minimal set of reactions
    # that we need to add for growth and call it good.
    if growth:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 media, biomass_eqtn)
        # print("Additional reactions required: " + str(additions) + "\n")
        print('reactions' + " : " + str(original_reactions.union(additions)))
        sys.exit(0)

    # propose other reactions that we have proteins for
    with_p_reactions = gapfill.suggest_reactions_with_proteins(reactions, True)
    # find the new reactions
    with_p_reactions.difference_update(reactions2run)
    added_reactions.append(("With proteins", with_p_reactions))
    reactions2run.update(with_p_reactions)
    status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass_eqtn)
    sys.stderr.write("After adding all reactions with proteins we get: " +
                     str(value) + " (growth is " + str(growth) + ")\n\n")

    # if this grows then we want to find the minimal set of reactions
    # that we need to add for growth and call it good.
    if growth:
        additions = resolve_additional_reactions(original_reactions, added_reactions, compounds, reactions,
                                                 media, biomass_eqtn)
        # print("Additional reactions required: " + str(additions) + "\n")
        print('reactions' + " : " + str(original_reactions.union(additions)))
        sys.exit(0)
