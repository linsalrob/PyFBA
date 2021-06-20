"""
Measure the fluxes in a set of reactions
"""
import argparse
import os
import sys

import PyFBA
from PyFBA import log_and_message


def fluxes(reactions, modeldata, media, biomass_equation, verbose=False):
    """
    Run the FBA and return the fluxes through each reaction
    :param biomass_equation: the biomass equation
    :type biomass_equation: PyFBA.metabolism.Reaction
    :param reactions: a set of the reactions to run
    :type reactions: set[str]
    :param modeldata: the model data object
    :type modeldata: PyFBA.model_seed.ModelData
    :param media: the media object
    :type media: set[PyFBA.metabolism.Compound]
    :param verbose: more output
    :return: a dict of the reactions and their fluxes
    :rtype: dict[str, float]
    """

    todelete = set()
    for r in reactions:
        if r not in modeldata.reactions:
            log_and_message(f"WARNING: Reaction {r} not found in our reaction set", stderr=verbose)
            todelete.add(r)
    for r in todelete:
        reactions.remove(r)

    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions, media, biomass_equation, verbose=verbose)
    if not growth:
        msg = f'ERROR: The set of {len(reactions)} reactions that you provided did not result in growth. '
        msg += 'We can not report fluxes if there was no growth!'
        log_and_message(msg, stderr=verbose)
        return {}
    return PyFBA.fba.reaction_fluxes()


def measure_fluxes():
    """
    Parse the arguments and start measuring the fluxes.
    """

    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Run Flux Balance Analysis and calculate reaction fluxes')
    parser.add_argument('-r', '--reactions', help='A list of the reactions in this model, one per line', required=True)
    parser.add_argument('-o', '--output', help='file to save the fluxes list to', required=True)
    parser.add_argument('-m', '--media', help='media name', required=True)
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-b', '--biomass', help='biomass equation to use. Default is the same as --type option')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    if not os.path.exists(args.reactions):
        sys.stderr.write(f"FATAL: {args.reactions} does not exist. Please check your files\n")
        sys.exit(1)

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)

    if args.verbose:
        log_and_message(f"Reading reactions from {args.reactions}", stderr=args.verbose)
    rxns = set()
    with open(args.reactions, 'r') as f:
        for li in f:
            if li.startswith('rxn'):
                rxns.add(li.strip())
            else:
                log_and_message(f'Skipped reaction {li} from {args.reactions} as it is not a standard reaction',
                                stderr=args.verbose)

    modeldata = PyFBA.parse.model_seed.parse_model_seed_data(args.type, verbose=args.verbose)
    if args.biomass:
        biomass_equation = PyFBA.metabolism.biomass_equation(args.biomass)
    else:
        biomass_equation = PyFBA.metabolism.biomass_equation(args.type)

    media = PyFBA.parse.pyfba_media(args.media, modeldata, args.verbose)
    fl = fluxes(reactions=rxns, modeldata=modeldata, media=media, biomass_equation=biomass_equation,
                verbose=args.verbose)

    if fl:
        with open(args.output, 'w') as out:
            for rx in fl:
                out.write(f"{rx}\t{fl[rx]}\n")
