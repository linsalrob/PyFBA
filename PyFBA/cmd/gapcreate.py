"""
Read a set of reactions and a media and recursively remove reactions to find the minimum set
required for growth on this media
"""
import copy
import os
import sys
import argparse

import PyFBA
from PyFBA import log_and_message, initiate_logger

model_data = PyFBA.model_seed.ModelData()


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
    log_and_message(f"FBA run {why} has a biomass flux value of {value:.2f} --> Growth: {growth}", stderr=verbose)
    return value, growth

def gap_create(reactions, media, model_data, orgtype='gramnegative',  flux_fraction = 0, verbose=False):
    """
    Iteratively remove reactions until we don't have any more to test
    :param reactions: list of reactins to test
    :type reactions: set[str]
    :param media: The media to test growth on
    :type media: Set[PyFBA.metabolism.Compound]
    :param model_data: the model seed data
    :type model_data: PyFBA.model_seed.ModelData
    :param orgtype: the type of organism
    :type orgtype: str
    :param flux_fraction: The fraction of the initial flux that we still consider growth. Set to 0 to just use growth
    :type flux_fraction: float
    :param verbose: more output
    :type verbose: bool
    :return: a set of reactions where we grow
    :rtype: set[str]
    """
    biomass_equation = PyFBA.metabolism.biomass_equation(orgtype)

    original_reactions = copy.deepcopy(reactions)
    initial_value, initial_growth = run_eqn("Initial", model_data, reactions, media, biomass_equation, verbose=verbose)
    if not initial_growth:
        log_and_message(f"The initial set of {len(reactions)} reactions doesn't grow on your media (flux: {initial_value})",
                        stderr = True)
        return reactions

    num = len(original_reactions)
    c=0
    for r in original_reactions:
        reactions.remove(r)
        c += 1
        value, growth = run_eqn(f"Reaction {c}: {r} ", model_data, reactions, media, biomass_equation, verbose=verbose)
        if flux_fraction > 0:
            if value/initial_value >= flux_fraction:
                # this is growth
                log_and_message(f"Reaction {c}/{num} ({r}) Flux: {value:.2f} Flux fraction {value/initial_value:.3f} {r} NOT required", stderr=verbose)
            else:
                log_and_message(f"Reaction {c}/{num} ({r}) Flux: {value:.2f} Flux fraction {value/initial_value:.3f} {r} REQUIRED", stderr=verbose)
                reactions.add(r)
        else:
            if growth:
                log_and_message(f"Reaction {c}/{num} ({r}) Flux: {value:.2f} {r} not required", stderr=verbose)
            else:
                log_and_message(f"Reaction {c}/{num} ({r}) is required for growth", stderr=verbose)
                reactions.add(r)
    log_and_message(f"After testing all the reactions, {len(reactions)} are required", stderr=verbose)
    return reactions

def create_reaction_gaps():
    """
    Create gaps in reactions while we still grow
    """
    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Import a list of reactions and then iterate through testing each'
                                                 'reaction to see if the model can still grow')
    parser.add_argument('-r', '--reactions', help='reactions file', required=True)
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-m', '--media', help='media name', required=True)
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-f', '--flux_fraction', default=0, type=float,
                        help='Flux fraction to consider grwoth. By default we use any flux but you can set it to e.g. 0.75 of the initial flux')
    parser.add_argument('-l', '--log', help='log file to write the detailed output (optional)', type=str)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    if not os.path.exists(args.reactions):
        sys.stderr.write(f"FATAL: {args.reactions} does not exist. Please check your files\n")
        sys.exit(1)

    if args.type not in orgtypes:
        sys.exit("Sorry, {} is not a valid organism type".format(args.type))

    if args.log:
        initiate_logger(args.log)
    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)

    # read the enzyme data
    # compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes(args.type)
    global model_data
    model_data = PyFBA.parse.model_seed.parse_model_seed_data(args.type, verbose=args.verbose)
    reactions = read_reactions(args.reactions, args.verbose)
    log_and_message(f"Found {len(reactions)} reactions", stderr=args.verbose)
    media = PyFBA.parse.read_media.find_media_file(args.media, model_data, args.verbose)
    rxn = gap_create(reactions, media, model_data, orgtype=args.type,
                     flux_fraction=args.flux_fraction, verbose=args.verbose)
    with open(args.output, 'w') as out:
        for r in rxn:
            out.write(f"{r}\n")
