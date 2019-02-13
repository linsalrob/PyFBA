
"""
A replacement parser for the Model SEED biochemistry modules that are available on Github
at https://github.com/ModelSEED/ModelSEEDDatabase.

The older version used the SOLR dumps to parse the data. This version uses the json
which is a lot cleaner.

We parse compounds from the compounds file in Biochemistry. Locations
are currently hardcoded because the ModelSeedDirectory does not contain
a mapping for compartments (the mapping files do not have the integers
used in the reactions file!).
"""

import copy
import os
import re
import sys
import io
import json

import PyFBA
from PyFBA.tools import bcolors

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


def template_reactions(modeltype='microbial'):
    """
    Load the template reactions to adjust the model. These are in the Templates directory, and just
    adjust some of the reactions to be specific for

    Returns a hash of some altered parameters for the model.
    :param modeltype: which type of model to load e.g. GramNegative, GramPositive, Microbial
    :type modeltype: str
    :return: A hash of the new model parameters that should be used to update the reactions object
    :rtype: dict
    """

    inputfile = ""
    if modeltype.lower() == 'microbial':
        inputfile = "Templates/Microbial/Reactions.tsv"
    elif modeltype.lower() == 'gramnegative' or modeltype.lower() == 'gram_negative':
        inputfile = "Templates/GramNegative/Reactions.tsv"
    elif modeltype.lower() == 'grampositive' or modeltype.lower() == 'gram_positive':
        inputfile = "Templates/GramPositive/Reactions.tsv"
    elif modeltype.lower() == 'mycobacteria':
        inputfile = "Templates/Mycobacteria/Reactions.tsv"
    elif modeltype.lower() == 'plant':
        inputfile = "Templates/Plant/Reactions.tsv"
    else:
        raise NotImplementedError("Parsing data for " + inputfile + " has not been implemented!")

    if not os.path.exists(os.path.join(MODELSEED_DIR, inputfile)):
        raise IOError(f"{bcolors.FAIL}FATAL: {bcolors.ENDC}" + os.path.join(MODELSEED_DIR, inputfile) +
                      " was not found. Please check your model SEED directory (" + MODELSEED_DIR + ")")

    new_enz = {}
    with open(os.path.join(MODELSEED_DIR, inputfile), 'r') as f:
        for l in f:
            if l.startswith('id'):
                continue
            p = l.strip().split("\t")
            new_enz[p[0]] = {}
            new_enz[p[0]]['direction'] = p[2]
            new_enz[p[0]]['enzymes'] = set(p[-1].split("|"))

    return new_enz

def compounds(compounds_file=None):
    """
    Load the compounds mapping. This maps from cpd id to name (we use
    the name in our reactions, but use the cpd id to parse the model
    seed database to avoid ambiguities.

    Optionally, you can provide a compounds file. If not, the default
    in MODELSEED_DIR/Biochemistry/compounds.json will be used.

    Note that the compounds file must be in json format. See model_seed_tsv for a tab separated parser.

    :param compounds_file: An optional filename of a compounds file to parse
    :type compounds_file: str
    :return: A hash of compounds with the str(compound) as the key and the compound object as the value
    :rtype: dict

    """

    cpds = {}

    if not compounds_file:
        compounds_file = os.path.join(MODELSEED_DIR, 'Biochemistry/compounds.json')

    try:
        with open(compounds_file, 'r') as f:
            data = json.load(f)

        for cpd in data:
            cc = PyFBA.metabolism.Compound(data[cpd]['name'], 'c')
            cc.model_seed_id = cpd
            ce = PyFBA.metabolism.Compound(p[2], 'e')
            ce.model_seed_id = cpd

            for k in data[cpd].keys():
                if data[cpd][k] == "null" or data[cpd][k] == "none":
                    data[cpd][k] = None

            for k in ["abbreviation", "abstract_compound", "aliases", "charge", "comprised_of", "deltag", "deltagerr",
                      "formula", "id", "inchikey", "is_cofactor", "is_core", "is_obsolete", "linked_compound", "mass",
                      "name", "pka", "pkb", "smiles", "source"]:
                cc.k = data[cpd][k]
                ce.k = data[cpd][k]

            # there are some compounds (like D-Glucose and Fe2+) that appear >1x in the table
            if str(cc) in cpds:
                cpds[str(cc)].alternate_seed_ids.add(cpd)
            else:
                cpds[str(cc)] = cc

            # there are some compounds (like D-Glucose and Fe2+) that appear >1x in the table
            if str(ce) in cpds:
                cpds[str(ce)].alternate_seed_ids.add(cpd)
            else:
                cpds[str(ce)] = ce
    except IOError as e:
        sys.stderr.write(f"{bcolors.FAIL}FATAL: {bcolors.ENDC} ")
        sys.exit("There was an error parsing " +
                 compounds_file + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))

    return cpds

def location():
    """Parse or return the codes for the locations. The ModelSEEDDatabase
    uses codes, and has a compartments file but they do not match up.

    This is currently hardcoded, but is put here so we can rewrite it as
    if the compartments file is updated.

    Note: It looks like the locations in the Model SEED files are specific to different organisms

    :return: A dict of location numeric IDs and string IDs
    :rtype: dict
    """

    # 0: cytoplasmic, 1: extracellular, 2: chloroplast

    global all_locations
    all_locations = {'0': 'c', '1': 'e', '2': 'h'}
    return all_locations


def reactions(organism_type="", rctf='Biochemistry/reactions.master.tsv', verbose=False):
    """
    Parse the reaction information in Biochemistry/reactions.json

    One reaction ID is associated with one equation and thus many
    compounds and parts.

    If the boolean verbose is set we will print out error/debugging
    messages.

    You can supply an alternative reactions file (rctf) if you
    don't like the default but this must be in json format. See model_seed_tsv.py
    to parse a tab separated values file.

    :param organism_type: The type of organism, eg. microbial, gram_negative, gram_positive
    :type organism_type: str
    :param rctf: The optional reaction file to provide
    :type rctf: str
    :param verbose: Print more output
    :type verbose: bool
    :return: Two components, a dict of the reactions and a dict of all the compounds used in the reactions.
    :rtype: dict, dict

    """

    locations = location()
    cpds = compounds()
    # cpds_by_id = {cpds[c].model_seed_id: cpds[c] for c in cpds}
    cpds_by_id = {"e": {}, "c": {}, "h": {}}
    for c in cpds:
        cpds_by_id[cpds[c].location][cpds[c].model_seed_id] = cpds[c]
        for asi in cpds[c].alternate_seed_ids:
            cpds_by_id[cpds[c].location][asi] = cpds[c]

    all_reactions = {}

    try:
        with open(os.path.join(MODELSEED_DIR, rctf), 'r') as rxnf:
            data = json.load(rxnf)

        for rid in data:
            rxn = data[rid]['equation']

            for k in data[rid].keys():
                if data[rid][k] == "null" or data[rid][k] == "none":
                    data[rid][k] = None

                # @TODO: check this:
                rxn = pieces[6]

                if data[rid]['deltag']:
                    deltaG = float(data[rid]['deltag'])
                else:
                    deltaG = 0.0
                if data[rid]['deltagerr']:
                    deltaG_error = float(data[rid]['deltagerr'])
                else:
                    deltaG_error = 0.0

                # we need to split the reaction, but different reactions
                # have different splits!

                separator = ""
                for separator in [" <=> ", " => ", " <= ", " = ", " < ", " > ", "Not found"]:
                    if separator in rxn:
                        break
                if separator == "Not found":
                    if verbose:
                        sys.stderr.write("WARNING: Could not find a seperator in " + rxn +
                                         ". This reaction was skipped. Please check it\n")
                    continue

                left, right = rxn.split(separator)

                # create a new reaction object to hold all the information ...

                r = PyFBA.metabolism.Reaction(rid)

                r.deltaG = deltaG
                r.deltaG_error = deltaG_error
                if data[rid]['is_transport'] != '0':
                    r.is_transport = True
                all_reactions[rid] = r

                r.direction = data[rid]['direction']

                # we have to rewrite the equation to accomodate
                # the proper locations
                newleft = []
                newright = []

                # deal with the compounds on the left side of the equation
                m = re.findall('\(([\d\.e-]+)\)\s+(.*?)\[(\d+)\]', left)
                if m == [] and verbose:
                    sys.stderr.write("ERROR: Could not parse the compounds" + " on the left side of the reaction " +
                                     rid + ": " + rxn + "\n")

                for p in m:
                    (q, cmpd, locval) = p

                    if locval in locations:
                        loc = locations[locval]
                    else:
                        if verbose:
                            sys.stderr.write("WARNING: Could not get a location " + " for " + locval + "\n")
                        loc = locval

                    # we first look up to see whether we have the compound
                    # and then we need to create a new compound with the
                    # appropriate location

                    if cmpd in cpds_by_id[loc]:
                        nc = cpds_by_id[loc][cmpd]
                    else:
                        if verbose:
                            sys.stderr.write("ERROR: Did not find " + cmpd + " in the compounds file.\n")
                        nc = PyFBA.metabolism.Compound(cmpd, loc)

                    ncstr = str(nc)
                    nc.add_reactions({rid})
                    cpds[ncstr] = nc

                    r.add_left_compounds({nc})
                    r.set_left_compound_abundance(nc, float(q))

                    newleft.append("(" + str(q) + ") " + nc.name + "[" + loc + "]")

                # deal with the right side of the equation
                m = re.findall('\(([\d\.e-]+)\)\s+(.*?)\[(\d+)\]', right)
                if m == [] and verbose:
                        sys.stderr.write("ERROR: Could not parse the compounds on the right side of the reaction " +
                                         rid + ": " + rxn + " >>" + right + "<<\n")

                for p in m:
                    (q, cmpd, locval) = p

                    if locval in locations:
                        loc = locations[locval]
                    else:
                        if verbose:
                            sys.stderr.write("WARNING: Could not get a location " + " for " + locval + "\n")
                        loc = locval

                    # we first look up to see whether we have the compound
                    # and then we need to create a new compound with the
                    # appropriate location

                    if cmpd in cpds_by_id[loc]:
                         nc = cpds_by_id[loc][cmpd]
                    else:
                        if verbose:
                            sys.stderr.write("ERROR: Did not find " + cmpd + " in the compounds file.\n")
                        nc = PyFBA.metabolism.Compound(cmpd, loc)

                    ncstr = str(nc)

                    nc.add_reactions({rid})
                    cpds[ncstr] = nc

                    r.add_right_compounds({nc})
                    r.set_right_compound_abundance(nc, float(q))

                    newright.append("(" + str(q) + ") " + nc.name + "[" + loc + "]")

                r.equation = " + ".join(newleft) + " <=> " + " + ".join(newright)

                # @TODO figure out which elements should be arrays. Need to do this for the compounds, too.
                r.aliases = data[rid]['aliases'].split(";")

                all_reactions[rid] = r
    except IOError as e:
        sys.exit("There was an error parsing " + rctf + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))

    # finally, if we need to adjust the organism type based on Template reactions, we shall
    if organism_type:
        new_rcts = template_reactions(organism_type)
        for r in new_rcts:
            all_reactions[r].direction = new_rcts[r]['direction']
            all_reactions[r].enzymes = new_rcts[r]['enzymes']

    return cpds, all_reactions

