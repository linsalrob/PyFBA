"""
Test all the reactions one at a time. We need to run theta(number of reactions) complexity to do this.

See PyFBA.gapgeneration.test_reactions.py for a O(n)/omega(log n) approach
"""
import copy
import os
import sys
import argparse

import PyFBA

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

    status, value, growth = PyFBA.fba.run_fba(modeldata, reactions_to_run, media, biomass_eqn)
    print("Before we test components, FBA has " + str(value) + " --> Growth: " + str(growth))
    if not growth:
        sys.exit("Since the complete model does not grow, we can't parse out the important parts!")

    ori_reactions = copy.copy(reactions_to_run)
    for r in ori_reactions:
        reactions_to_run = copy.copy(ori_reactions)
        reactions_to_run.remove(r)
        PyFBA.fba.remove_uptake_and_secretion_reactions(modeldata.reactions)
        status, value, growth = PyFBA.fba.run_fba(modeldata, reactions_to_run, media, biomass_eqn)
        print("{}\t{}".format(r, growth))
