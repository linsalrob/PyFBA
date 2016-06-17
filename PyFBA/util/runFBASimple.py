#!/usr/bin/python2.7
from __future__ import print_function
import sys
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("reactions", help="Reactions file")
parser.add_argument("media", help="Media file")
parser.add_argument("-v", "--verbose", help="Verbose stderr output",
                    action="store_true")
parser.add_argument("--header", help="Header line in reaction file",
                    action="store_true")
parser.add_argument("--dev", help="Use PyFBA dev code",
                    action="store_true")
args = parser.parse_args()

vrb = args.verbose
if args.dev:
    # Import PyFBA from absolute path
    sys.path.insert(0, os.path.expanduser("~") + "/Projects/PyFBA/")
    sys.path.insert(0, os.path.expanduser("~") + "/PyFBA/")
    print("IN DEV MODE", file=sys.stderr)
import PyFBA

# Load ModelSEED database
compounds, reactions, enzymes =\
    PyFBA.parse.model_seed.compounds_reactions_enzymes('gramnegative')

# Read in reactions file
modelReactions = set()
with open(args.reactions, "r") as f:
    if args.header:
        next(f)  # Header line
    for l in f:
        l = l.strip()
        ll = l.split("\t")
        modelReactions.add(ll[0])

if vrb:
    print("Loaded", len(modelReactions), "reactions", file=sys.stderr)

# Remove reactions IDs that do not not have a reaction equation associated
tempset = set()
for r in modelReactions:
    if r in reactions:
        tempset.add(r)
    elif args.verbose:
        print("Reaction ID", r,
              "is not in our reactions list. Skipped",
              file=sys.stderr)
modelReactions = tempset
tempset = None

# Load our media
media = PyFBA.parse.read_media_file(args.media)
print("Our media has", len(media) ,"components", file=sys.stderr)

# Define a biomass equation
biomass_equation = PyFBA.metabolism.biomass_equation('gramnegative')

# Run FBA
status, value, growth = PyFBA.fba.run_fba(compounds, reactions,
                                          modelReactions,
                                          media, biomass_equation)
print("The biomass reaction has a flux of", value,
      "--> Growth:", growth, file=sys.stderr)

# Save flux values
for rxn, val in PyFBA.lp.col_primal_hash().items():
    print(rxn, val, sep="\t")
