import argparse
import os
import sys

import fba
from parse import model_seed, read_media_file
from metabolism import biomass

import metabolism

__author__ = 'Rob Edwards'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Import a list of reactions and run the FBA')
    parser.add_argument('-r', help='reactions file', required=True)
    parser.add_argument('-m', help='media file', required=True)
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
                    sys.stderr.write("Biomass reaction was skipped from the list as it is imported\n")
                continue
            r = l.strip()
            if r in reactions:
                reactions2run.add(r)

    media = read_media_file(args.m)
    bme = biomass.biomass_equation('gramnegative')

    status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, bme, True)
    print("Initial run has " + str(value) + " --> Growth: " + str(growth))
