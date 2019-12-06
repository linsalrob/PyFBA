"""
This code is designed to extract InCHI keys for the relevant compounds in a metabolic model to be used in
metabolomic analysis. Currently only usable with the ModelSEED Database.
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
    sys.stderr.write(f"{bcolors.FAIL}FATAL: {bcolors.ENDC}")
    sys.stderr.write("Please ensure that you install the Model SEED Database somewhere, and set the environment " +
                     "variable ModelSEEDDatabase to point to that directory.\n" +
                     " See INSTALLATION.md for more information\n")
    sys.exit(-1)

if not MODELSEED_DIR:
    sys.stderr.write(f"{bcolors.FAIL}FATAL: {bcolors.ENDC}")
    sys.stderr.write("The ModelSEEDDatabase environment variable is not set.\n")
    sys.stderr.write("Please install the ModelSEEDDatabase, set the variable, and try again")
    sys.exit(-1)

if not os.path.exists(MODELSEED_DIR):
    sys.stderr.write(f"{bcolors.FAIL}FATAL: {bcolors.ENDC}")
    sys.stderr.write("The MODEL SEED directory: {} does not exist.\n".format(MODELSEED_DIR))
    sys.stderr.write("Please check your installation.\n")
    sys.exit(-1)

parser = argparse.ArgumentParser(description="Create a list of InCHI keys from compounds in an sbml file")
parser.add_argument('-i', help='sbml input file', required=True)
parser.add_argument('-o', help='output file')
parser.add_argument('-v', help='verbose output', action='store_true')
args = parser.parse_args()

if not os.path.exists(args.i):
    sys.exit("SBML file {} was not found".format(args.i))

out = "{}.inchi.tsv".format(args.i[:-5])
if args.o:
    out = args.o

sbml = PyFBA.parse.parse_sbml_file(args.i, False)
cids = {}
ckeys = {}

compounds_file = os.path.join(MODELSEED_DIR, 'Biochemistry/compounds.json')

# creating a dict of each compound ID to its InCHI key
try:
    with open(compounds_file, 'r') as f:
        mseed = json.load(f)

    for cpd in mseed:
        key = mseed[cpd]['inchikey']
        ckeys[cpd] = key
except:
    print('Could not read compounds.json; please set the MODELSEED_DIR path variable.')

# get the compounds from the SBML file. This is a dict of metabolite.Compound objects
cpdids = sbml.compounds_by_id
if args.v:
    sys.stderr.write("List of 'cpds': {}\n".format(len(cpds)))

# creating a dict of each compound in the sbml matching its name to its id
for cpd_id, cpd_obj in cpdids.items():
    try:
        newcm, cpdname, cpdloc = cpd_id.split("_")
    except:
        cpdname, cpdloc = cpd_id.split("_")
    cpd_name = cpd_obj.name
    cids[cpdname] = cpd_name

with open(out, 'w') as outfile:
    for id, name in cids.items():
        inkey = ckeys[id]
        try:
            if inkey:
                outfile.write("{}\t{}\t{}\n".format(id, name, inkey))
            else:
                outfile.write("{}\t{}\tnone\n".format(id, name))
        except:
            print("{}: {} has no InCHI key".format(id, name))


