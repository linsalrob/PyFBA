
"""
Parse the modelSEED json file, and create a table of all the compounds. This is largely a utility for our
friends using modelSEED to analyze MS/MS data.
"""

import argparse
import os
import sys
import json

# from PyFBA.tools import bcolors
from roblib import bcolors

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


def compounds(compounds_file=None, outfile="compounds.tsv", verbose=False):
    """
    Load the compounds mapping. This maps from cpd id to name (we use
    the name in our reactions, but use the cpd id to parse the model
    seed database to avoid ambiguities.

    Optionally, you can provide a compounds file. If not, the default
    in MODELSEED_DIR/Biochemistry/compounds.json will be used.

    Note that the compounds file must be in json format. See model_seed_tsv for a tab separated parser.

    :param compounds_file: An optional filename of a compounds file to parse
    :type compounds_file: str
    :param outfile: Output file to write
    :param verbose: more output
    :return: Nothing

    """

    if not compounds_file:
        compounds_file = os.path.join(MODELSEED_DIR, 'Biochemistry/compounds.json')

    if (os.path.exists(outfile)):
        sys.stderr.write(f"{bcolors.RED}FATAL: {outfile} already exists. Not overwriting{bcolors.ENDC}\n")
        sys.exit(1)

    out = open(outfile, 'w')

    adducts = {
        "M+H": +1.007825,
        "M+Na": +22.98977,
        "M+K": +38.96371,
        "M+2H20+H": +37.02896,
        "M+C2H3N+H": +42.03437,
        "M+C2H3N+Na": +64.01632,
        "M+2Na-H": +44.97172,
        "M+H+Na": +23.9976,
        "M-H": -1.007825,
        "M+Na-2H": +20.93812,
        "M+K-2H": +36.91206,
        "M+CH3COO(-)": +59.013853,
        "M+Cl(-)": +34.968304,
        "M-H2O+H": -17.002474,
    }

    alladducts = sorted(adducts.keys())



    out.write("Compound\tName\tFormula\tInchikey\tCharge\tpKa\tpKb\tSMILES\tMass")
    for a in alladducts:
        out.write(f"\t{a}")
    out.write('\t2M+H\tM+2H\tM-2H\n')

    try:
        with open(compounds_file, 'r') as f:
            data = json.load(f)
    except IOError as e:
        sys.stderr.write(f"{bcolors.FAIL}FATAL: {bcolors.ENDC} ")
        sys.exit("There was an error parsing " +
             compounds_file + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))

    for cpd in data:
        out.write("\t".join(map(str, [cpd, data[cpd]['name'], data[cpd]['formula'], data[cpd]['inchikey'],
                                      data[cpd]['charge'], data[cpd]["pka"],
                                      data[cpd]["pkb"], data[cpd]["smiles"], data[cpd]['mass']])))
        if data[cpd]['mass'] and 'null' != data[cpd]['mass']:
            mymass = float(data[cpd]['mass'])
            for a in alladducts:
                out.write("\t{}".format(mymass + adducts[a]))
            out.write("\t{}".format((2*mymass) + 1.007825)) # 2M + H
            out.write("\t{}".format((mymass + 2.05165)/2))  # M+2H
            out.write("\t{}".format((mymass - 2.05165)/2))  # M-2H
        out.write("\n")

    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a tsv of all compounds and their adducts')
    parser.add_argument('-o', help='output file to write', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    compounds(None, args.o, args.v)