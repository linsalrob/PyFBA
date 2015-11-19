import argparse
import json
import os
import sys
import fba
from parse import model_seed, read_media_file
from metabolism import biomass_equation
import metabolism

__author__ = 'Rob Edwards'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Dump an FBA run to JSON, keeping only the relevant parts')
    parser.add_argument('-r', help='reactions file (required)', required=True)
    parser.add_argument('-m', help='media file (required)', required=True)
    parser.add_argument('-j', help='json output file (required)', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    # read the enzyme data
    compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes('gramnegative')
    reactions2run = set()
    with open(args.r, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue
            if "biomass_equation" in l:
                if args.v:
                    sys.stderr.write("Biomass reaction was skipped from the list as it is imported\n")
                continue
            r = l.strip()
            if r in reactions:
                reactions2run.add(r)

    media = read_media_file(args.m)
    biomass_equation = biomass_equation.biomass_equation('gramnegative')

    # trim the reactions to only those ones that are used in the model
    reactions = {r: reactions[r] for r in reactions2run}
    with open(args.j, 'w') as out:
        json.dump({'reactions': reactions, 'reactions_to_run': reactions2run, 'compounds': compounds,
                   'media': media, 'biomass_equation': biomass_equation}, out)

    status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass_equation, verbose=True)
    print("Initial run has " + str(value) + " --> Growth: " + str(growth))
