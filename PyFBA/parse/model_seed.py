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

import os
import re
import sys
import json
from importlib.resources import open_text

from typing import Dict, Set

import PyFBA

from PyFBA import MODELSEED_DIR
from PyFBA import Biochemistry

from PyFBA.model_seed import ModelSeed
from PyFBA import log_and_message

modelseedstore = ModelSeed()


def template_reactions(modeltype):
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
        for li in f:
            if not li.startswith('id'):
                p = li.strip().split("\t")
                new_enz[p[0]] = {}
                new_enz[p[0]]['direction'] = p[2]
                new_enz[p[0]]['enzymes'] = set(p[-1].split("|"))

    return new_enz


def compounds(compounds_file='Biochemistry/compounds.json', verbose=False) -> Set[PyFBA.metabolism.Compound]:
    """
    Load the compounds mapping that connects ID to compound objects

    Optionally, you can provide a compounds file. If not, the default
    in MODELSEED_DIR/Biochemistry/compounds.master.tsv will be used.

    :param compounds_file: An optional filename of a compounds file to parse
    :type compounds_file: str
    :parma verbose: more output
    :type verbose: bool
    :return: A set of Compound objects.
    :rtype: Set[PyFBA.metabolism.Compound]

    """

    global modelseedstore

    if modelseedstore.compounds:
        return modelseedstore.compounds

    modelseedstore.compounds = set()


    """
    try:
        if verbose:
            sys.stderr.write(f"Parsing compounds in {compounds_file}\n")
        with open(os.path.join(MODELSEED_DIR, compounds_file), 'r') as f:
            for jc in json.load(f):
                c = PyFBA.metabolism.Compound(jc['id'], jc['name'])
                c.model_seed_id = jc['id']
                c.mw = jc['mass']

                # parse the aliases. At the time of writing, the aliases are f'd up
                # and look like this:
                # 'Name: Manganese; Manganese(2+); Mn(II); Mn++; Mn+2; Mn2+; manganese (II) ion; manganese ion'
                # i.e. they have converted the hash but not encoded it
                if 'aliases' in jc and jc['aliases']:
                    allals = {}
                    for al in jc['aliases']:
                        parts = al.split(': ')
                        keys = parts[0]
                        vals = parts[1].split('; ')
                        allals[keys] = vals
                    c.add_attribute('aliases', allals)


                # this should be all the keys (except name and ID)
                for ck in ["abbreviation", "abstract_compound", 
                           "charge", "comprised_of", "deltag",
                           "deltagerr", "formula", "id", "inchikey",
                           "is_cofactor", "is_core", "is_obsolete",
                           "linked_compound", "mass", "name", "notes",
                           "pka", "pkb", "smiles", "source"]:
                    if ck in jc:
                        c.add_attribute(ck, jc[ck])

                modelseedstore.compounds.add(c)
    except IOError as e:
        sys.exit("There was an error parsing " +
                 compounds_file + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))
    """

    f = open_text("Biochemistry", 'ModelSEEDDatabase.Biochemistry.compounds.json')
    for jc in json.load(f):
        c = PyFBA.metabolism.Compound(jc['id'], jc['name'])
        c.model_seed_id = jc['id']
        c.mw = jc['mass']

        # parse the aliases. At the time of writing, the aliases are f'd up
        # and look like this:
        # 'Name: Manganese; Manganese(2+); Mn(II); Mn++; Mn+2; Mn2+; manganese (II) ion; manganese ion'
        # i.e. they have converted the hash but not encoded it
        if 'aliases' in jc and jc['aliases']:
            allals = {}
            for al in jc['aliases']:
                parts = al.split(': ')
                keys = parts[0]
                vals = parts[1].split('; ')
                allals[keys] = vals
            c.add_attribute('aliases', allals)

        # this should be all the keys (except name and ID)
        for ck in ["abbreviation", "abstract_compound",
                   "charge", "comprised_of", "deltag",
                   "deltagerr", "formula", "id", "inchikey",
                   "is_cofactor", "is_core", "is_obsolete",
                   "linked_compound", "mass", "name", "notes",
                   "pka", "pkb", "smiles", "source"]:
            if ck in jc:
                c.add_attribute(ck, jc[ck])

        modelseedstore.compounds.add(c)

    return modelseedstore.compounds


def location() -> Dict[str, str]:
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


def reactions(organism_type=None, rctf='Biochemistry/reactions.json', verbose=False) \
        -> Dict[str, PyFBA.metabolism.Reaction]:
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
    :return: a dict of the reactions
    :rtype: Dict[str, PyFBA.metabolism.Reaction]
    """

    if not organism_type:
        if verbose:
            sys.stderr.write("ERROR: A model type was not specified, and so using microbial core")
        organism_type = "Core"

    global modelseedstore


    if modelseedstore.organism_type and modelseedstore.organism_type == organism_type and modelseedstore.reactions:
        return modelseedstore.reactions

    modelseedstore.organism_type = organism_type
    cpds = compounds(verbose=verbose)
    locations = location()

    modelseedstore.reactions = {}  # type Dict[Any, Reaction]

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

                # we parse out the two sides and compare them.
                # we do left and store in new[0] and right and store in new[1] and then rejoin
                new = [[],[]]
                for i, side in enumerate([left, right]):
                    side = side.strip()
                    if not side:
                        continue

                    # deal with the compounds on the left side of the equation
                    m = re.findall(r'\(([\d.e-]+)\)\s+(.*?)\[(\d+)]', side)
                    if m == [] and verbose:
                        sys.stderr.write("ERROR: Could not parse the compounds from {side}\n")

                    for p in m:
                        (q, cmpd, locval) = p

                        if locval in locations:
                            loc = locations[locval]
                        else:
                            if verbose:
                                sys.stderr.write("WARNING: Could not get a location for {locval}\n")
                            loc = locval

                        # we first look up to see whether we have the compound
                        # and then we need to create a new compound with the
                        # appropriate location
                        msg = f"Looking for {cmpd}: "
                        cpdbyid = modelseedstore.get_compound_by_id(cmpd)
                        if cpdbyid:
                            nc = PyFBA.metabolism.CompoundWithLocation.from_compound(cpdbyid, loc)
                            msg += " found by ID"
                        else:
                            cpdbyname = modelseedstore.get_compound_by_name(cmpd)
                            if cpdbyname:
                                nc = PyFBA.metabolism.CompoundWithLocation.from_compound(cpdbyname, loc)
                                msg += " found by name"
                            else:
                                nc = PyFBA.metabolism.CompoundWithLocation(cmpd, cmpd, loc)
                                msg += " not found"
                        if verbose:
                            log_and_message(msg, stderr=True)
                        nc.add_reactions({r.id})

                        if i == 0:
                            r.add_left_compounds({nc})
                            r.set_left_compound_abundance(nc, float(q))
                        else:
                            r.add_right_compounds({nc})
                            r.set_right_compound_abundance(nc, float(q))

                        new[i].append("(" + str(q) + ") " + nc.name + "[" + loc + "]")

                r.equation = " + ".join(new[0]) + " <=> " + " + ".join(new[1])

                modelseedstore.reactions[r.id] = r

    # finally, if we need to adjust the organism type based on Template reactions, we shall
    new_rcts = template_reactions(organism_type)
    for r in new_rcts:
        modelseedstore.reactions[r].direction = new_rcts[r]['direction']
        modelseedstore.reactions[r].enzymes = new_rcts[r]['enzymes']

    return modelseedstore.reactions


def ftr_to_roles(rf="Annotations/Roles.tsv") -> Dict[str, str]:
    """
    Read the roles file and create a dictionary of feature_id->role
    :param rf: the Roles file
    :return: a dict of feature_ids and roles
    """

    ftr2role = {}
    with open(os.path.join(MODELSEED_DIR, rf), 'r') as f:
        for li in f:
            if li.startswith('id'):
                continue
            p = li.strip().split("\t")
            ftr2role[p[0]] = p[1]
    return ftr2role


def complex_to_ftr(cf="Annotations/Complexes.tsv") -> Dict[str, set]:
    """
    Read the complexes file, and create a dict of complex->feature_ids
    :param cf: the Complexes file
    :return: a dict of complexes -> feature_ids
    """

    cpx2ftr = {}
    with open(os.path.join(MODELSEED_DIR, cf), 'r') as f:
        for li in f:
            if li.startswith('id'):
                continue
            p = li.strip().split("\t")
            cpx2ftr[p[0]] = set()
            if 'null' == p[5]:
                continue
            for ftr in p[5].split('|'):
                cpx2ftr[p[0]].add(ftr.split(';')[0])

    return cpx2ftr


def enzymes(organism_type="", verbose=False) -> Dict[str, PyFBA.metabolism.Enzyme]:
    """
    Convert each of the roles and complexes into a set of enzymes, and connect them to reactions.

    We return three dicts, the compounds, the reactions, and the enzymes. Compounds and reactions are as above,
    enzymes are from this method.

    :param organism_type: The type of organism, eg. Microbial, Gram_positive, Gram_negative
    :type organism_type:str
    :param verbose:Print more output
    :type verbose:bool
    :return: The compounds, the reactions, and the enzymes in that order
    :rtype: dict of Enzymes
    """

    # The complex (Enzyme) is our key data structure as it connects reactions and roles.
    # So we need to go from complexes->roles and complexes->reactions
    # The complex->roles is through the two annotation files
    # The complex->reactions is through the Templates file

    global modelseedstore

    if modelseedstore.enzymes:
        return modelseedstore.enzymes

    rcts = reactions(organism_type, verbose=verbose)

    modelseedstore.enzymes = {}

    # Set up enzymes with complexes and reactions
    c2f = complex_to_ftr()
    f2r = ftr_to_roles()
    for cmplx in c2f:
        if cmplx in modelseedstore.enzymes:
            if verbose:
                sys.stderr.write(f"Warning: have duplicate {cmplx} complexes that maybe in more than once. " +
                                 "Skipped later incantations\n")
            continue
        modelseedstore.enzymes[cmplx] = PyFBA.metabolism.Enzyme(cmplx)
        for ft in c2f[cmplx]:
            if ft in f2r:
                modelseedstore.enzymes[cmplx].add_roles({f2r[ft]})
            else:
                sys.stderr.write(f"Warning: No functional role for {ft}\n")
            for ecno in re.findall(r'[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+', f2r[ft]):
                modelseedstore.enzymes[cmplx].add_ec(ecno)
    for r in rcts:
        for c in rcts[r].enzymes:
            if c in modelseedstore.enzymes:
                modelseedstore.enzymes[c].add_reaction(rcts[r].id)

    return modelseedstore.enzymes


def compounds_reactions_enzymes(organism_type='', verbose=False) -> (Dict[str, PyFBA.metabolism.Compound],
                                                                     Dict[str, PyFBA.metabolism.Reaction],
                                                                     Dict[str, PyFBA.metabolism.Enzyme]):
    """
    Return three dicts, compounds, reactions, functions
    :param organism_type: The type of organism (e.g. Microbial, GramNegative)
    :param verbose: more output
    :return: dict, dict, dict
    """
    return (
        compounds(verbose=verbose),
        reactions(organism_type=organism_type, verbose=verbose),
        enzymes(organism_type=organism_type, verbose=verbose)
    )


def complexes(organism_type='', verbose=False) -> Dict[str, Set[PyFBA.metabolism.Reaction]]:
    """
    Generate a list of the complexes in the SEED data. Connection between complexes and reactions. A complex can be
    involved in many reactions.
    :return: a dict with key is complex and value is all reactions
    """

    global modelseedstore
    if not modelseedstore.complexes:
        enz = enzymes(organism_type=organism_type, verbose=verbose)
        modelseedstore.complexes = {e: enz[e].reactions for e in enz}

    return modelseedstore.complexes


def roles(organism_type='', verbose=False) -> Dict[str, Set[str]]:
    """
    Return a hash of the roles where the id is the role name and the value is the set of complex IDs that the role is
    inolved in
    :param organism_type: limit to a type of organism
    :param verbose: more output
    :return: a dict of role->complexes
    """

    if not modelseedstore.roles:
        enz = enzymes(organism_type=organism_type, verbose=verbose)
        modelseedstore.roles = {}
        for e in enz:
            for r in enz[e].roles:
                if r not in modelseedstore.roles:
                    modelseedstore.roles[r] = set()
                modelseedstore.roles[r].add(e)
    return modelseedstore.roles

def parse_model_seed_data(organism_type='', verbose=False):
    """
    Parse the model seed data and return a ModelSeed class that contains
    all the data
    :param organism_type: limit to a type of organism
    :param verbose: more output
    :return: a ModelSeed class
    :rtype: PyFBA.model_seed.ModelSeed
    """

    global modelseedstore
    
    if not modelseedstore.compounds:
        c = compounds(verbose=verbose),
    if not modelseedstore.reactions:
        r = reactions(organism_type=organism_type, verbose=verbose),
    if not modelseedstore.enzymes:
        e = enzymes(organism_type=organism_type, verbose=verbose)
    if not modelseedstore.complexes:
        c = complexes(organism_type=organism_type, verbose=verbose)
    if not modelseedstore.roles:
        r = roles(organism_type=organism_type, verbose=verbose)
    return modelseedstore

def reset_cache():
    """
    Reset the cache of modelseed data to force reparsing it
    """

    global modelseedstore
    modelseedstore = ModelSeed()
    modelseedstore.reset()

