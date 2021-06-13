"""
Test which reactions are essential for growth on a media.
"""
import argparse
import copy
import sys

import PyFBA


def test_growth(reactions_to_delete, reactions_to_run, modeldata, media, biomass_eqn, verbose):
    """
    Test growth of reactions_to_run after we have deleted reactions_to_delete. Returns True on growth, False on no growth

    :param reactions_to_delete: reactions to delete from reactions to run
    :type reactions_to_delete: set
    :param reactions_to_run: Reactions to run for the model
    :type reactions_to_run: set
    :param modeldata: The modeldata object
    :param PyFBA.model_seed.ModelData
    :param media: The media object
    :type media: set
    :param biomass_eqn: The biomass equation
    :type biomass_eqn: PyFBA.metabolism.reaction.Reaction
    :param verbose: Print more output
    :type verbose: bool
    :return: Whether the remaining reactions result in growth
    :rtype: bool
    """

    new_r2r = set([x for x in reactions_to_run if x not in reactions_to_delete])
    PyFBA.fba.remove_uptake_and_secretion_reactions(reactions)

    status, value, growth = PyFBA.fba.run_fba(modeldata, new_r2r, media, biomass_eqn)

    if verbose:
        sys.stderr.write("Deleted {} rxns. Use {}. Growth: {}\n".format(len(reactions_to_delete), len(new_r2r), growth))

    return growth


def not_essential_reactions(reactions_to_delete, reactions_to_run, modeldata, media, biomass_eqn, verbose):
    """
    Iterate through the reactions and return the minimal set that are/are  not essential

    We split reactions_to_delete in half, and test. While we don't get growth, we keep splitting. If we get to a single
    element and we don't get growth we don't return it.

    It is most likely that when we split a model in half, we won't get growth. Therefore, looking for the halves that
    give us growth.

    :param reactions_to_delete: reactions to delete from reactions to run
    :type reactions_to_delete: set
    :param reactions_to_run: Reactions to run for the model
    :type reactions_to_run: set
    :param modeldata: The modeldata object
    :param PyFBA.model_seed.ModelData
    :param media: The media object
    :type media: set
    :param biomass_eqn: The biomass equation
    :type biomass_eqn: PyFBA.metabolism.reaction.Reaction
    :param verbose: Print more output
    :type verbose: bool
    :return: Whether the remaining reactions result in growth
    :rtype: Set[str]
    """

    # if we have a one element list, we need to test it, and either return it if there is growth or return an empty
    # set if there is not growth
    if len(reactions_to_delete) == 1:
        if test_growth(reactions_to_delete, reactions_to_run, modeldata, media, biomass_eqn, verbose):
            return set(reactions_to_delete)
        else:
            return set()

    left, right = PyFBA.gapfill.bisect(list(reactions_to_delete))
    # test left to see if every element is redundant
    redundant_elements = set()
    if test_growth(left, reactions_to_run, modeldata, media, biomass_eqn, verbose):
        # we get growth
        redundant_elements.update(left)
    else:
        # test the left half again
        redundant_elements.update(
            not_essential_reactions(left, reactions_to_run, modeldata, media, biomass_eqn, verbose)
        )

    # now test the right half
    if test_growth(right, reactions_to_run, modeldata, media, biomass_eqn, verbose):
        # we get growth
        redundant_elements.update(right)
    else:
        # test the right half again
        redundant_elements.update(
            not_essential_reactions(right, reactions_to_run, modeldata, media, biomass_eqn, verbose)
        )
    return redundant_elements


def test_all_reactions(reactions_to_run, modeldata, media, biomass_eqn, verbose):
    """
    Test all the reactions and print out which are required and which are redundant

    :param reactions_to_run: Reactions to run for the model
    :type reactions_to_run: set
    :param modeldata: The modeldata object
    :param PyFBA.model_seed.ModelData
    :param media: The media object
    :type media: set
    :param biomass_eqn: The biomass equation
    :type biomass_eqn: PyFBA.metabolism.reaction.Reaction
    :param verbose: Print more output
    :type verbose: bool
    :return: Whether the remaining reactions result in growth
    :rtype: bool
    """

    todelete = copy.deepcopy(reactions_to_run)
    # prevent inadvertent edits to reactions_to_run
    reactions_to_run = frozenset(reactions_to_run)
    redundant = not_essential_reactions(todelete, reactions_to_run, modeldata, media, biomass_eqn, verbose)

    for r in reactions_to_run:
        if r in redundant:
            print("{}\tREDUNDANT".format(r))
        else:
            print("{}\tESSENTIAL".format(r))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Test all reactions in a model")
    parser.add_argument('-r', help='reactions file', required=True)
    parser.add_argument('-m', help='media file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    # read the enzyme data
    modeldata = PyFBA.parse.model_seed.parse_model_seed_data('gramnegative', verbose=args.v)
    reactions_to_run = set()
    with open(args.r, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue
            if "biomass_equation" in l:
                if args.v:
                    sys.stderr.write("Biomass reaction was skipped from the list as it is imported\n")
                continue
            r = l.strip()
            if r in modeldata.reactions:
                reactions_to_run.add(r)

    media = PyFBA.parse.read_media_file(args.m)
    biomass_eqn = PyFBA.metabolism.biomass_equation('gramnegative')

    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions_to_run, media, biomass_eqn, verbose=args.v)
    sys.stderr.write("Before we test components, FBA has " + str(value) + " --> Growth: " + str(growth))
    if not growth:
        sys.exit("Since the complete model does not grow, we can't parse out the important parts!")

    test_all_reactions(reactions_to_run, modeldata, media, biomass_eqn, args.v)