"""
Given two media conditions, one where we grow and one where we do not,
and a set of reactions, test all the reactions to see which is required
for growth
"""

import os
import sys
import argparse
import copy

import PyFBA
from PyFBA import log_and_message

model_data = PyFBA.model_seed.ModelData()

__author__ = 'Rob Edwards'


def read_reactions(reaction_file, verbose=False):
    """
    Read the reactions file and return a set of reactions
    :param reaction_file: the reactions file
    :param verbose: more output
    :return: a set of reactions
    :rtype: set[str]
    """

    global model_data
    reactions2run = set()
    with open(reaction_file, 'r') as f:
        for li in f:
            if li.startswith('#'):
                continue
            if "biomass" in li.lower():
                if verbose:
                    sys.stderr.write("Biomass reaction was skipped from the list as it is auto-imported\n")
                continue
            r = li.strip()
            if r in model_data.reactions:
                reactions2run.add(r)
    return reactions2run

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


def compare_media(reactions, positive, negative, model_data, orgtype='gramnegative',  verbose=False):
    """
    Iteratively remove reactions until we don't have any more to test
    :param reactions: list of reactins to test
    :type reactions: set[str]
    :param positive: The media that we can grow on
    :type positive: Set[PyFBA.metabolism.Compound]
    :param negative: The media that we can not grow on
    :type negative: Set[PyFBA.metabolism.Compound]
    :param model_data: the model seed data
    :type model_data: PyFBA.model_seed.ModelData
    :param orgtype: the type of organism
    :type orgtype: str
    :param verbose: more output
    :type verbose: bool
    :return: a set of reactions where we grow
    :rtype: set[str]
    """
    biomass_equation = PyFBA.metabolism.biomass_equation(orgtype)

    original_reactions = copy.deepcopy(reactions)
    value, growth = run_eqn("Initial", model_data, reactions, positive, biomass_equation, verbose=verbose)
    if not growth:
        log_and_message(f"The initial set of {len(reactions)} reactions doesn't grow on your media {positive} (growth: {growth})",
                        stderr = True)
        return reactions

    num = len(original_reactions)
    c=0
    both = set()
    nonly = set()
    ponly = set()
    neither = set()
    log_and_message("REACTION ID\tPositive Media\tNegative Media", stderr = True)
    for r in original_reactions:
        reactions.remove(r)
        c += 1
        pvalue, pgrowth = run_eqn("Initial", model_data, reactions, positive, biomass_equation, verbose=verbose)
        nvalue, ngrowth = run_eqn("Initial", model_data, reactions, negative, biomass_equation, verbose=verbose)
        log_and_message(f"{r}\t{pgrowth}\t{ngrowth}", stderr = True)
        if pgrowth and ngrowth:
            both.add(r)
        elif pgrowth:
            ponly.add(r)
        elif ngrowth:
            nonly.add(r)
        else:
            neither.add(r)
        reactions.add(r)

    newreactions = both.union(ponly)
    pvalue, pgrowth = run_eqn("Initial", model_data, newreactions, positive, biomass_equation, verbose=verbose)
    nvalue, ngrowth = run_eqn("Initial", model_data, newreactions, negative, biomass_equation, verbose=verbose)
    log_and_message(f"After refining the reactions to {len(newreactions)} combined reactions we end up with: " +
                    f"{positive}: growth {pgrowth} and {negative}: growth {ngrowth} ", stderr = True)
    return newreactions

def compare_two_media():
    """
    Compare two media and test every reaction for growth/no growth
    """
    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Import a list of reactions and then compare growth on two '
                                                 'media conditions to identify essential/non-essential media')
    parser.add_argument('-r', '--reactions', help='reactions file', required=True)
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-p', '--positive', help='media name where we should grow', required=True)
    parser.add_argument('-n', '--negative', help='media name where we should not grow', required=True)
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    if not os.path.exists(args.reactions):
        sys.stderr.write(f"FATAL: {args.reactions} does not exist. Please check your files\n")
        sys.exit(1)

    if args.type not in orgtypes:
        sys.exit("Sorry, {} is not a valid organism type".format(args.type))

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)

    # read the enzyme data
    # compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes(args.type)
    global model_data
    model_data = PyFBA.parse.model_seed.parse_model_seed_data(args.type, verbose=args.verbose)
    reactions = read_reactions(args.reactions, args.verbose)
    log_and_message(f"Found {len(reactions)} reactions", stderr=args.verbose)
    positive = PyFBA.parse.read_media.find_media_file(args.positive, model_data, args.verbose)
    negative = PyFBA.parse.read_media.find_media_file(args.negative, model_data, args.verbose)

    rxn = compare_media(reactions, positive, negative, model_data, args.type, args.verbose)
    with open(args.output, 'w') as out:
        for r in rxn:
            out.write(f"{r}\n")
