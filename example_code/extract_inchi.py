"""
This code is designed to extract InCHI keys for the relevant compounds in a metabolic model to be used in
metabolomic analysis.
"""

import argparse
import copy
import os
import sys
import json

import PyFBA

MODELSEED_DIR = ""
if 'ModelSEEDDatabase' in os.environ:
        MODELSEED_DIR = os.environ['ModelSEEDDatabase']
else:
    sys.stderr.write("Please ensure that you install the Model SEED Database somewhere, and set the environment " +
                     "variable ModelSEEDDatabase to point to that directory.\n" +
                     " See INSTALLATION.md for more information\n")
    sys.exit(-1)

if not MODELSEED_DIR:
    sys.stderr.write("The ModelSEEDDatabase environment variable is not set.\n")
    sys.stderr.write("Please install the ModelSEEDDatabase, set the variable, and try again")
    sys.exit(-1)

if not os.path.exists(MODELSEED_DIR):
    sys.stderr.write("The MODEL SEED directory: {} does not exist.\n".format(MODELSEED_DIR))
    sys.stderr.write("Please check your installation.\n")
    sys.exit(-1)

parser = argparse.ArgumentParser(description="Create a list of InCHI keys from compounds in an sbml/reaction file")
parser.add_argument('-s', help='sbml input file')
parser.add_argument('-r', help='reactions input file')
parser.add_argument('-o', help='output file')
parser.add_argument('-v', help='verbose output', action='store_true')
args = parser.parse_args()

# Checking that an input is specified
if not args.s and not args.r:
    sys.exit("Either -s (sbml file) or -r (reactions file) must be specified.")

# Checking for double input
if args.s and args.r:
    sys.exit("Only one input type at a time please.")

# Mapping reactions, compounds, and InChI keys
cnames = {}
ckeys = {}
rcpds = {}

try:
    with open(os.path.join(MODELSEED_DIR, 'Biochemistry/compounds.json'), 'r') as cfile:
        compounds = json.load(cfile)
        for cpd in compounds:
            cnames[cpd] = compounds[cpd]['name']
            ckeys[cpd] = compounds[cpd]['inchikey']

    with open(os.path.join(MODELSEED_DIR, 'Biochemistry/reactions.json'), 'r') as rfile:
        reactions = json.load(rfile)
        for rxn in reactions:
            cpds = reactions[rxn]['compound_ids']
            cpds = cpds.strip().split(';')
            rcpds[rxn] = cpds

except:
    print('Failed to parse the ModelSEED Biochemistry Database.')

# Attempt to parse an SBML
if args.s:
    try:
        sbml = PyFBA.parse.parse_sbml_file(args.s, False)

        if args.o:
            out = args.o
        else:
            out = "{}.inchi.txt".format(args.s)

        with open(out, 'w') as outfile:
            for cpd in sbml.compounds_by_id:
                parts = cpd.split('_')
                c = parts[0]
                name = cnames.get(c, "none")
                key = ckeys.get(c, "none")
                outfile.write(f"{c}\t{name}\t{key}\n")

    except:
            sys.exit("There was an error parsing the SBML.")

# Attempt to parse a reactions file
if args.r:
    try:
        cpds = set()
        with open(args.r, 'r') as infile:
            for line in infile:
                rxn = line.strip()
                rxn_cpds = rcpds.get(rxn, '')
                for c in rxn_cpds:
                    cpds.add(c)

        if args.o:
            out = args.o
        else:
            out = "{}.inchi.txt".format(args.r)

        with open(out, 'w') as outfile:
            for cpd in cpds:
                name = cnames.get(cpd, "none")
                key = ckeys.get(cpd, "none")
                outfile.write(f"{cpd}\t{name}\t{key}\n")

    except:
            sys.exit("There was an error parsing the reactions file.")
