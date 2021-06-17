"""
Read a set of assigned functions, and then run the Flux Balance Analysis on the file
"""

import sys
import argparse
import PyFBA


def read_af(assigned_functions):
    pass


def read_reactions(r2rfile, verbose=False):
    reactions_to_run = set()
    with open(r2rfile, 'r') as f:
        for li in f:
            if li.startswith('rxn'):
                reactions_to_run.add(li.strip())
    if verbose:
        print(f"There are {len(reactions_to_run)} reactions to run")
    return reactions_to_run


def run_fba(reactions_to_run, mediafile, verbose=False):
    modeldata = PyFBA.parse.model_seed.parse_model_seed_data('gramnegative', verbose=verbose)
    if verbose:
        print(f"There are {len(modeldata.compounds):,} compounds, {len(modeldata.reactions):,} reactions, "
              f"and {len(modeldata.enzymes):,} enzymes in total")
    biomass_equation = PyFBA.metabolism.biomass_equation('gramnegative', modeldata.compounds)
    # Read the media file and correct the names
    media = PyFBA.parse.pyfba_media(mediafile, modeldata, verbose=verbose)
    print(f"The media has {len(media)} compounds")

    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions_to_run, media, biomass_equation,
                                              {}, verbose=verbose)
    print(f"Running FBA on {mediafile} completed with a flux value of {value} --> growth: {growth}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-m', help='media to choose [Default %(default)s]', default="ArgonneLB")
    parser.add_argument('-r', help='reactions file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    r2r = None
    if args.r:
        r2r = read_reactions(args.r, args.v)
    else:
        sys.stderr.write("Errr, need reactions somehow\n")
        sys.exit(-1)

    run_fba(r2r, args.m, args.v)
