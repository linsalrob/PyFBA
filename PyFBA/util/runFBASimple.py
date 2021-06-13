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
modeldata = PyFBA.parse.model_seed.parse_model_seed_data('gramnegative', verbose=True)

# Read in reactions file
modelReactions = set()
modelInfo = {}
with open(args.reactions, "r") as f:
    if args.header:
        next(f)  # Header line
    for l in f:
        l = l.strip()
        rxn, func, eqn, gf = l.split("\t")
        modelReactions.add(rxn)
        modelInfo[rxn] = (func, eqn, gf)

if vrb:
    print("Loaded", len(modelReactions), "reactions", file=sys.stderr)

# Remove reactions IDs that do not not have a reaction equation associated
tempset = set()
for r in modelReactions:
    if r in modeldata.reactions:
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
status, value, growth = PyFBA.fba.run_fba(modeldata,
                                          modelReactions,
                                          media, biomass_equation)
print("The biomass reaction has a flux of", value,
      "--> Growth:", growth, file=sys.stderr)

# Save flux values
print("reaction\tflux\tfunction\tequation")
for rxn, val in PyFBA.lp.col_primal_hash().items():
    func, eqn, gf = modelInfo[rxn]
    print(rxn, val, func, eqn, sep="\t")
