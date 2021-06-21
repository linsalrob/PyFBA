"""
Read a set of reactions and a media and recursively remove reactions to find the minimum set
required for growth on this media
"""
import copy
import os
import sys
import argparse

import PyFBA
from PyFBA import log_and_message

model_data = PyFBA.model_seed.ModelData()


def read_media(mediafile, modeldata, verbose=False):
    """
    Read the media file and return a set of compounds
    :param modeldata: the modeldata object
    :type modeldata: PyFBA.model_seed.ModelData
    :param mediafile: the media file to read
    :param verbose: more output
    :type verbose: bool
    :return: a set of media compounds
    :rtype: Set[PyFBA.metabolism.Compound]
    """

    if mediafile in PyFBA.parse.media_files():
        log_and_message(f"parsing media directly from {mediafile}", stderr=verbose)
        # pyfba media already corrects the names, so we can  just return it.
        return PyFBA.parse.pyfba_media(mediafile, modeldata)
    elif os.path.exists(mediafile):
        log_and_message(f"parsing media file {mediafile}", stderr=verbose)
        media = PyFBA.parse.read_media_file(mediafile)
    elif 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)):
        log_and_message(f"parsing media file {os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)}", stderr=verbose)
        media = PyFBA.parse.read_media_file(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile))
    else:
        log_and_message(f"Can't figure out how to parse media from {mediafile}", stderr=True, loglevel="CRITICAL")
        sys.exit(-1)
    return PyFBA.parse.correct_media_names(media, modeldata.compounds)


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

def gap_create(reactions, media, model_data, orgtype='gramnegative',  verbose=False):
    """
    Iteratively remove reactions until we don't have any more to test
    :param reactions:
    :type reactions:
    :param media:
    :type media:
    :param model_data:
    :type model_data:
    :param orgtype:
    :type orgtype:
    :param verbose:
    :type verbose:
    :return:
    :rtype:
    """
    biomass_equation = PyFBA.metabolism.biomass_equation(orgtype)

    original_reactions = copy.deepcopy(reactions)
    value, growth = run_eqn("Initial", model_data, reactions, media, biomass_equation, verbose=verbose)
    if not growth:
        log_and_message(f"The initial set of {len(reactions)} reactions doesn't grow on your media (growth: {growth})",
                        stderr = True)
        return reactions

    num = len(original_reactions)
    c=0
    for r in original_reactions:
        reactions.remove(r)
        c += 1
        value, growth = run_eqn("Initial", model_data, reactions, media, biomass_equation, verbose=verbose)
        if growth:
            log_and_message(f"Reaction {c}/{num} ({r}) is required for growth", stderr=verbose)
            reactions.add(r)
        else:
            log_and_message(f"Reaction {c}/{num} ({r}) {r} not required", stderr=verbose)
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
    media = read_media(args.media, model_data, args.verbose)
    rxn = gap_create(reactions, media, model_data, args.type, args.verbose)
    with open(args.output, 'w') as out:
        for r in rxn:
            out.write(f"{r}\n")
