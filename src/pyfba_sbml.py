"""
Run PyFBA on an SBML file.

Much of this code is taken from the complimentary jupyter notebook, Using_an_SBML_model
"""
import os
import argparse
import sys

import PyFBA
from PyFBA import log_and_message


def run_fba_sbml(sbmlfile, medianame, verbose=False):
    """
    Run the FBA on the SBML file and return the output
    :param sbmlfile: The SBML file to parse
    :param medianame: The name of the media
    :param verbose: More output
    :return: the flux and whether the model grew
    """

    sbml = PyFBA.parse.parse_sbml_file(sbmlfile)

    # Get a dict of reactions.
    # The key is the reaction ID, and the value is a metabolism.reaction.Reaction object
    reactions = sbml.reactions
    reactions_to_run = set()
    uptake_secretion_reactions = {}
    biomass_equation = None
    for r in reactions:
        if 'biomass_equation' == r or 'Biomass' == reactions[r].readable_name:
            biomass_equation = reactions[r]
            print(f"Our biomass equation is {biomass_equation.readable_name}")
            continue
        is_boundary = False
        for c in reactions[r].all_compounds():
            if c.uptake_secretion:
                is_boundary = True
                break
        if is_boundary:
            reactions[r].is_uptake_secretion = True
            uptake_secretion_reactions[r] = reactions[r]
        else:
            reactions_to_run.add(r)

    all_compounds = sbml.compounds
    # Filter for compounds that are boundary compounds
    filtered_compounds = set()
    for c in all_compounds:
        if not c.uptake_secretion:
            filtered_compounds.add(c)

    ms = PyFBA.model_seed.ModelData(compounds=filtered_compounds, reactions=reactions)
    # Read the media file
    media = PyFBA.parse.pyfba_media(medianame, ms)

    # Adjust the lower bounds of uptake secretion reactions
    # for things that are not in the media
    for u in uptake_secretion_reactions:
        reactions[u].lower_bound = -1000.0
        uptake_secretion_reactions[u].lower_bound = -1000.0
        reactions[u].upper_bound = 1000.0
        uptake_secretion_reactions[u].upper_bound = 1000.0

        is_media_component = False
        for c in uptake_secretion_reactions[u].all_compounds():
            if c in media:
                is_media_component = True
        if not is_media_component:
            reactions[u].lower_bound = 0.0
            uptake_secretion_reactions[u].lower_bound = 0.0

    if verbose:
        log_and_message(f"The biomass equation is {biomass_equation}", stderr=verbose)
        log_and_message(f"There are {len(reactions)} reactions in the model", stderr=verbose)
        log_and_message(f"There are {len(uptake_secretion_reactions)} uptake/secretion reactions in the model",
                        stderr=verbose)
        log_and_message(f"There are {len(reactions_to_run)} reactions to be run in the model", stderr=verbose)
        log_and_message(f"There are {len(all_compounds)} total compounds in the model", stderr=verbose)
        log_and_message(f"There are {len(filtered_compounds)} compounds that are not involved in uptake and secretion",
                        stderr=verbose)
        log_and_message(f"There are {len(media)} media compounds in {medianame}", stderr=verbose)

    status, value, growth = PyFBA.fba.run_fba(ms, reactions_to_run, media, biomass_equation, uptake_secretion_reactions,
                                              verbose=verbose)
    print(f"The FBA on {medianame} completed with a flux value of {value} --> growth: {growth}")
    return value, growth


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse an SBML file and run flux balance analysis')
    parser.add_argument('-s', help='SBML input file', required=True)
    parser.add_argument('-m', help='media to choose [Default %(default)s]', default="ArgonneLB")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.s):
        sys.stderr.write(f"ERROR: {args.s} was not found. Please check the path and try again\n")
        sys.exit(-1)

    if args.m not in PyFBA.Biochemistry.media:
        sys.stderr.write(f"Error: {args.m} is not a recognised media format. Please choose from this list:\n")
        sys.stderr.write("\n".join(PyFBA.Biochemistry.media.keys()))
        sys.exit(-1)

    run_fba_sbml(args.s, args.m, args.v)
