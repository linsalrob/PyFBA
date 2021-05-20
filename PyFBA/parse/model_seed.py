"""
A parser for the SEED biochemistry modules that are available on Github
at https://github.com/ModelSEED/ModelSEEDDatabase. We have also included
them in our repo as a submodule.

We parse compounds from the compounds file in Biochemistry. Locations
are currently hardcoded because the ModelSeedDirectory does not contain
a mapping for compartments (the mapping files do not have the integers
used in the reactions file!).

PLEASE NOTE: This version uses the ModelSeed JSON files, please revert
to model_seed_solr.py if you want to use the old SOLR datadumps.
You should be able to do that by changing the functions in __init__.py

"""

import copy
import os
import re
import sys
import json
import PyFBA

from .config import MODELSEED_DIR


def template_reactions(modeltype='Microbial'):
    """
    Load the template reactions to adjust the model. Returns a hash of some altered parameters for the model
    :param modeltype: which type of model to load e.g. GramNegative, GramPositive, Microbial
    :type modeltype: str
    :return: A hash of the new model parameters that should be used to update the reactions object
    :rtype: dict
    """

    inputfile = ""
    if modeltype.lower() == 'core':
        inputfile = "Templates/Core/Reactions.tsv"
    elif modeltype.lower() == 'fungi':
        inputfile = "Templates/Fungi/Reactions.tsv"
    elif modeltype.lower() == 'gramnegative' or modeltype.lower() == 'gram_negative':
        inputfile = "Templates/GramNegative/Reactions.tsv"
    elif modeltype.lower() == 'grampositive' or modeltype.lower() == 'gram_positive':
        inputfile = "Templates/GramPositive/Reactions.tsv"
    elif modeltype.lower() == 'human':
        inputfile = "Templates/Human/Reactions.tsv"
    elif modeltype.lower() == 'microbial':
        inputfile = "Templates/Microbial/Reactions.tsv"
    elif modeltype.lower() == 'mycobacteria':
        inputfile = "Templates/Mycobacteria/Reactions.tsv"
    elif modeltype.lower() == 'plant':
        inputfile = "Templates/Plant/Reactions.tsv"
    else:
        raise NotImplementedError("Parsing data for " + inputfile + " has not been implemented!")

    if not os.path.exists(os.path.join(MODELSEED_DIR, inputfile)):
        raise FileExistsError(f"{inputfile} was not found. Please check your model SEED directory {MODELSEED_DIR}")

    new_enz = {}
    with open(os.path.join(os.path.join(MODELSEED_DIR, inputfile)), 'r') as f:
        for l in f:
            if not l.startswith('id'):
                p = l.strip().split("\t")
                new_enz[p[0]] = {}
                new_enz[p[0]]['direction'] = p[2]
                new_enz[p[0]]['enzymes'] = set(p[-1].split("|"))

    return new_enz


def compounds(compounds_file=None):
    """
    Load the compounds mapping. This maps from cpd id to name (we use
    the name in our reactions, but use the cpd to parse the model
    seed database to avoid ambiguities.

    Optionally, you can provide a compounds file. If not, the default
    in MODELSEED_DIR/Biochemistry/compounds.master.tsv will be used.

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
            for jc in json.load(f):
                c = PyFBA.metabolism.Compound(jc['name'], '')
                c.model_seed_id = jc['id']
                c.mw = jc['mass']
                # this should be all the keys (except name and ID)
                for ck in ["abbreviation", "abstract_compound", 
                           "aliases", "charge", "comprised_of", "deltag",
                           "deltagerr", "formula", "id", "inchikey",
                           "is_cofactor", "is_core", "is_obsolete",
                           "linked_compound", "mass", "name", "notes",
                           "pka", "pkb", "smiles", "source"]:
                    if ck in jc:
                        c.add_attribute(ck, jc[ck])

                # there are some compounds (like D-Glucose and Fe2+) that appear >1x in the table
                if str(c) in cpds:
                    cpds[str(c)].alternate_seed_ids.add(jc['id'])
                else:
                    cpds[str(c)] = c
    except IOError as e:
        sys.exit("There was an error parsing " +
                 compounds_file + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))

    return cpds


def location():
    """Parse or return the codes for the locations. The ModelSEEDDatabase
    uses codes, and has a compartments file but they do not match up.

    This is currently hardcoded, but is put here so we can rewrite it as
    if the compartments file is updated

    :return: A dict of location numeric IDs and string IDs
    :rtype: dict
    """

    # 0: cytoplasmic, 1: extracellular, 2: chloroplast

    all_locations = {'0': 'c', '1': 'e', '2': 'h'}
    return all_locations


def reactions(organism_type="", rctf='Biochemistry/reactions.json', verbose=False):
    """
    Parse the reaction information in Biochemistry/reactions.master.tsv

    One reaction ID is associated with one equation and thus many
    compounds and parts.

    If the boolean verbose is set we will print out error/debugging
    messages.

    You can supply an alternative reactions file (rctf) if you
    don't like the default.

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
    cpds_by_id = {}
    for c in cpds:
        cpds_by_id[cpds[c].model_seed_id] = cpds[c]
        for asi in cpds[c].alternate_seed_ids:
            cpds_by_id[asi] = cpds[c]

    all_reactions = {}  # type Dict[Any, Reaction]

    with open(os.path.join(MODELSEED_DIR, rctf), 'r') as rxnf:
        for rxn in json.load(rxnf):
            r = PyFBA.metabolism.Reaction(rxn['id'])
            if 'name' in rxn and rxn['name']:
                r.readable_name = rxn['name']
            r.model_seed_id = rxn['id']
            # convert a few 0/1 to True/False
            if rxn['is_transport']:
                r.is_transport = True
            if rxn['is_obsolete']:
                r.is_obsolete = True

            for rxnkey in ["abbreviation", "abstract_reaction", "aliases", "code", "compound_ids", "definition",
                           "deltag", "deltagerr", "direction", "ec_numbers",
                           "linked_reaction", "notes", "pathways", "reversibility",
                           "source", "status", "stoichiometry"]:
                if rxnkey in rxn:
                    r.add_attribute(rxnkey, rxn[rxnkey])

                for separator in [" <=> ", " => ", " <= ", " = ", " < ", " > ", "Not found"]:
                    if separator in rxn['equation']:
                        break
                if separator == "Not found":
                    if verbose:
                        sys.stderr.write("WARNING: Could not find a seperator in " + rxn['equation'] +
                                         ". This reaction was skipped. Please check it\n")
                    continue

                left, right = rxn['equation'].split(separator)

                # check and see we have a valid equation
                left = left.strip()
                right = right.strip()

                # we have to rewrite the equation to accomodate
                # the proper locations
                newleft = []
                newright = []

                # deal with the compounds on the left side of the equation
                m = re.findall(r'\(([\d.e-]+)\)\s+(.*?)\[(\d+)]', left)
                if m == [] and verbose:
                    sys.stderr.write("ERROR: Could not parse the compounds" + " on the left side of the reaction " +
                                     r.id + ": " + rxn['equation'] + "\n")

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

                    if cmpd in cpds_by_id:
                        nc = PyFBA.metabolism.Compound(cpds_by_id[cmpd].name, loc)
                    else:
                        if verbose:
                            sys.stderr.write("ERROR: Did not find " + cmpd + " in the compounds file.\n")
                        nc = PyFBA.metabolism.Compound(cmpd, loc)

                    ncstr = str(nc)

                    if ncstr in cpds:
                        nc = copy.copy(cpds[ncstr])
                    nc.add_reactions({r.id})
                    cpds[ncstr] = nc

                    r.add_left_compounds({nc})
                    r.set_left_compound_abundance(nc, float(q))

                    newleft.append("(" + str(q) + ") " + nc.name + "[" + loc + "]")

                # deal with the right side of the equation
                m = re.findall(r'\(([\d.e-]+)\)\s+(.*?)\[(\d+)]', right)
                if m == [] and verbose:
                    sys.stderr.write("ERROR: Could not parse the compounds on the right side of the reaction " +
                                     r.id + ": " + rxn['equation'] + " >>" + right + "<<\n")

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

                    if cmpd in cpds_by_id:
                        nc = PyFBA.metabolism.Compound(cpds_by_id[cmpd].name, loc)
                    else:
                        if verbose:
                            sys.stderr.write("ERROR: Did not find " + cmpd + " in the compounds file.\n")
                        nc = PyFBA.metabolism.Compound(cmpd, loc)

                    ncstr = str(nc)
                    if ncstr in cpds:
                        nc = copy.copy(cpds[ncstr])
                    nc.add_reactions({r.id})
                    cpds[ncstr] = nc

                    r.add_right_compounds({nc})
                    r.set_right_compound_abundance(nc, float(q))

                    newright.append("(" + str(q) + ") " + nc.name + "[" + loc + "]")

                r.equation = " + ".join(newleft) + " <=> " + " + ".join(newright)

                all_reactions[r.id] = r

    # finally, if we need to adjust the organism type based on Template reactions, we shall
    if not organism_type:
        sys.stderr.write("ERROR: A model type was not specified, and so your Enzyme Complexes are just microbial core")
        organism_type = "Core"

    new_rcts = template_reactions(organism_type)
    for r in new_rcts:
        all_reactions[r].direction = new_rcts[r]['direction']
        all_reactions[r].enzymes = new_rcts[r]['enzymes']

    return cpds, all_reactions


def ftr_to_roles(rf="Annotations/Roles.tsv"):
    """
    Read the roles file and create a dictionary of feature_id->role
    :param rf: the Roles file
    :return: a dict of feature_ids and roles
    """

    ftr2role = {}
    with open(rf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            ftr2role[p[0]] = p[1]
    return ftr2role


def complex_to_ftr(cf="Annotations/Complexes.tsv"):
    """
    Read the complexes file, and create a dict of complex->feature_ids
    :param cf: the Complexes file
    :return: a dict of complexes -> feature_ids
    """

    cpx2ftr = {}
    with open(cf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            cpx2ftr[p[0]] = set()
            if 'null' == p[5]:
                continue
            for ftr in p[5].split('|'):
                cpx2ftr[p[0]].add(ftr.split(';')[0])

    return cpx2ftr


def compounds_reactions_enzymes(organism_type='', verbose=False):
    """
    Convert each of the roles and complexes into a set of enzymes, and connect them to reactions.

    We return three dicts, the compounds, the reactions, and the enzymes. Compounds and reactions are as above,
    enzymes are from this method.

    :param organism_type: The type of organism, eg. Microbial, Gram_positive, Gram_negative
    :type organism_type:str
    :param verbose:Print more output
    :type verbose:bool
    :return: The compounds, the reactions, and the enzymes in that order
    :rtype: dict of Compound, dict of Reaction, dict of Enzyme
    """

    # The complex (Enzyme) is our key data structure as it connects reactions and roles.
    # So we need to go from complexes->roles and complexes->reactions
    # The complex->roles is through the two annotation files
    # The complex->reactions is through the Templates file

    cpds, rcts = reactions(organism_type, verbose=verbose)  # type
    enzs = {}

    # Set up enzymes with complexes and reactions
    c2f = complex_to_ftr()
    f2r = ftr_to_roles()
    for cmplx in c2f:
        if cmplx in enzs:
            if verbose:
                sys.stderr.write(f"Warning: have duplicate {cmplx} complexes that maybe in more than once. " +
                                 "Skipped later incantations\n")
            continue
        enzs[cmplx] = PyFBA.metabolism.Enzyme(cmplx)
        for ft in c2f[cmplx]:
            enzs[cmplx].add_roles({f2r[ft]})
            for ecno in re.findall(r'[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+', f2r[ft]):
                enzs[cmplx].add_ec(ecno)
    for r in rcts:
        for c in r.enzymes:
            if c in enzs:
                enzs[c].add_reaction(r.id)
