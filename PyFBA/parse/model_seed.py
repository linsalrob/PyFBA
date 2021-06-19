"""
A parser for the SEED biochemistry modules that are available on Github
at https://github.com/ModelSEED/ModelSEEDDatabase. We have also included
them in our repo as a submodule.

We parse compounds from the compounds file in Biochemistry. Locations
are currently hardcoded because the ModelSeedDirectory does not contain
a mapping for compartments (the mapping files do not have the integers
used in the reactions file!).

PLEASE NOTE: This version uses the ModelData JSON files, please revert
to model_seed_solr.py if you want to use the old SOLR datadumps.
You should be able to do that by changing the functions in __init__.py

"""

import sys
import re
import json
try:
    from importlib.resources import open_text
except ImportError:
    # this is for python<3.7
    from importlib_resources import open_text
from typing import Dict, Set, Any, List

import PyFBA
from PyFBA.model_seed import ModelData
from PyFBA import log_and_message

modelseedstore = ModelData()


def template_reactions(modeltype):
    """
    Load the template reactions to adjust the model. Returns a hash of some altered parameters for the model
    :param modeltype: which type of model to load e.g. GramNegative, GramPositive, Microbial
    :type modeltype: str
    :return: A hash of the new model parameters that should be used to update the reactions object
    :rtype: dict
    """

    if modeltype.lower() == 'core':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.Core"
    elif modeltype.lower() == 'fungi':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.Fungi"
    elif modeltype.lower() == 'gramnegative' or modeltype.lower() == 'gram_negative':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.GramNegative"
    elif modeltype.lower() == 'grampositive' or modeltype.lower() == 'gram_positive':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.GramPositive"
    elif modeltype.lower() == 'human':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.Human"
    elif modeltype.lower() == 'microbial':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.Microbial"
    elif modeltype.lower() == 'mycobacteria':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.Mycobacteria"
    elif modeltype.lower() == 'plant':
        inputmodule = "PyFBA.Biochemistry.ModelSEEDDatabase.Templates.Plant"
    else:
        raise NotImplementedError(f"Parsing data for {modeltype} has not been implemented!")

    reactionsf = open_text(inputmodule, "Reactions.tsv")
    new_enz = {}
    for li in reactionsf:
        if not li.startswith('id'):
            p = li.strip().split("\t")
            new_enz[p[0]] = {}
            new_enz[p[0]]['direction'] = p[2]
            new_enz[p[0]]['gfdirection'] = p[3]
            new_enz[p[0]]['enzymes'] = set(p[-1].split("|"))
    reactionsf.close()

    return new_enz


def compounds(compounds_file='compounds.json', verbose=False) -> Set[PyFBA.metabolism.Compound]:
    """
    Load the compounds mapping that connects ID to compound objects

    Optionally, you can provide a compounds file. If not, the default
    in MODELSEED_DIR/Biochemistry/compounds.master.tsv will be used.

    Note there are some compounds that apppear more than once in the ModelSEED json file.

    For example, "mecillinam" appears as both "cpd30724" and "cpd35839". The former has "source: Orphan" while the
    latter has "source: Primary Database". We prefer the latter, and keep the other ID(s) as alternate IDs.

    We use name as a primary key, although we should consider smiles or something else, these are often NULL in
    non-primary entries.

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

    log_and_message(f"Reading compounds from PyFBA.Biochemistry.ModelSEEDDatabase.Biochemistry.{compounds_file}",
                    stderr=verbose)
    compf = open_text("PyFBA.Biochemistry.ModelSEEDDatabase.Biochemistry", compounds_file)

    primary_compounds: Dict[str, PyFBA.metabolism.Compound] = {}
    secondary_compounds: Dict[str, List[PyFBA.metabolism.Compound]] = {}

    for jc in json.load(compf):
        # If we are not the primary source, and this compound already has a primary source,
        # we just append the ID and move along. Otherwise we need to make a new compound
        # in case we don't have a Primary Database source for this compound.

        # Initially we checked if this was not a Primary Database source, but there are some compounds
        # like D-Glucose that are in the database twice as primary compounds, e.g cpd00027 and cpd26821
        # So now we just take the first instance!
        # if jc['source'] != "Primary Database" and jc['name'] in primary_compounds:
        if jc['name'] in primary_compounds:
            primary_compounds[jc['name']].alternate_seed_ids.add(jc['id'])
            continue

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
                   "linked_compound", "mass", "notes",
                   "pka", "pkb", "smiles", "source"]:
            if ck in jc:
                c.add_attribute(ck, jc[ck])

        if jc['source'] == "Primary Database":
            if jc['name'] in secondary_compounds:
                for s in secondary_compounds[jc['name']]:
                    c.alternate_seed_ids.add(s.id)
                del secondary_compounds[jc['name']]
            primary_compounds[jc['name']] = c
        else:
            if jc['name'] not in secondary_compounds:
                secondary_compounds[jc['name']] = []
            secondary_compounds[jc['name']].append(c)

    # now just flatten secondary compounds
    if len(secondary_compounds) > 0:
        log_and_message(f"We found {len(secondary_compounds)} compounds that do not have a Primary Database equivalent")
        # at the moment we just use the last one as the exemplar
        for n in secondary_compounds:
            exc = secondary_compounds[n].pop()
            # check if we have additional information
            for extra in secondary_compounds[n]:
                exc.alternate_seed_ids.add(extra.id)
                for ck in ["abbreviation", "abstract_compound",
                           "charge", "comprised_of", "deltag",
                           "deltagerr", "formula", "id", "inchikey",
                           "is_cofactor", "is_core", "is_obsolete",
                           "linked_compound", "mass", "notes",
                           "pka", "pkb", "smiles", "source"]:
                    if not exc.get_attribute(ck) and  extra.get_attribute(ck):
                        exc.add_attribute(ck, extra.get_attribute(ck))
            primary_compounds[n] = exc

    modelseedstore.compounds = set(primary_compounds.values())
    modelseedstore.rebuild_indices()


    compf.close()
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


def reactions(organism_type=None, rctf='reactions.json', verbose=False) \
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
            log_and_message("ERROR: A model type was not specified, and so using microbial core", stderr=verbose)
        organism_type = "Core"

    global modelseedstore

    if modelseedstore.organism_type and modelseedstore.organism_type == organism_type and modelseedstore.reactions:
        return modelseedstore.reactions

    modelseedstore.organism_type = organism_type
    compounds(verbose=verbose)
    locations = location()

    modelseedstore.reactions = {}  # type Dict[Any, Reaction]

    log_and_message(f"Reading reactions from PyFBA.Biochemistry.ModelSEEDDatabase.Biochemistry.{rctf}",
                    stderr=verbose)
    rxnf = open_text("PyFBA.Biochemistry.ModelSEEDDatabase.Biochemistry", rctf)
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

        r.direction = rxn['direction']
        r.ntdirection = rxn['direction']
        if rxn['ec_numbers']:
            r.ec_numbers = rxn['ec_numbers']

        for rxnkey in ["abbreviation", "abstract_reaction", "aliases", "code", "compound_ids", "definition",
                       "deltag", "deltagerr", "linked_reaction", "notes", "pathways", "reversibility",
                       "source", "status", "stoichiometry"]:
            if rxnkey in rxn:
                r.add_attribute(rxnkey, rxn[rxnkey])

        for separator in [" <=> ", " => ", " <= ", " = ", " < ", " > ", "Not found"]:
            if separator in rxn['equation']:
                break
        if separator == "Not found":
            if verbose:
                log_and_message("WARNING: Could not find a seperator in {rxn['equation']} "
                                "This reaction was skipped. Please check it", stderr=verbose)
            continue

        left, right = rxn['equation'].split(separator)

        # we parse out the two sides and compare them.
        # we do left and store in new[0] and right and store in new[1] and then rejoin
        new = [[], []]
        for i, side in enumerate([left, right]):
            side = side.strip()
            if not side:
                continue

            # deal with the compounds on the left side of the equation
            m = re.findall(r'\(([\d.e-]+)\)\s+(.*?)\[(\d+)]', side)
            if m == [] and verbose:
                log_and_message("ERROR: Could not parse the compounds from {side}", stderr=verbose)

            for p in m:
                (q, cmpd, locval) = p

                if locval in locations:
                    loc = locations[locval]
                else:
                    if verbose:
                        log_and_message(f"WARNING: Could not get a location for {locval}", stderr=verbose)
                    loc = locval

                # we first look up to see whether we have the compound
                # and then we need to create a new compound with the
                # appropriate location
                msg = f"Looking for {cmpd}: "
                cpdbyid = modelseedstore.get_compound_by_id(cmpd)
                if cpdbyid:
                    nc = PyFBA.metabolism.CompoundWithLocation.from_compound(cpdbyid, loc)
                    msg = None # no real interest in those that we just find!
                else:
                    cpdbyname = modelseedstore.get_compound_by_name(cmpd)
                    if cpdbyname:
                        nc = PyFBA.metabolism.CompoundWithLocation.from_compound(cpdbyname, loc)
                        msg += " found by name"
                    else:
                        nc = PyFBA.metabolism.CompoundWithLocation(cmpd, cmpd, loc)
                        msg += " not found"
                if msg and verbose:
                    log_and_message(msg, stderr=verbose)
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

    rxnf.close()
    # finally, if we need to adjust the organism type based on Template reactions, we shall
    new_rcts = template_reactions(organism_type)
    for r in new_rcts:
        modelseedstore.reactions[r].direction = new_rcts[r]['direction']
        modelseedstore.reactions[r].gfdirection = new_rcts[r]['gfdirection']
        modelseedstore.reactions[r].enzymes = new_rcts[r]['enzymes']

    return modelseedstore.reactions


def ftr_to_roles(rf="Roles.tsv", verbose=False) -> Dict[str, str]:
    """
    Read the roles file and create a dictionary of feature_id->role
    :param rf: the Roles file
    :return: a dict of feature_ids and roles
    """

    ftr2role = {}
    log_and_message(f"Reading roles from PyFBA.Biochemistry.ModelSEEDDatabase.Biochemistry.{rf}", stderr=verbose)
    rolesf = open_text("PyFBA.Biochemistry.ModelSEEDDatabase.Annotations", rf)
    for li in rolesf:
        if li.startswith('id'):
            continue
        p = li.strip().split("\t")
        ftr2role[p[0]] = p[1]
    rolesf.close()

    return ftr2role


def complex_to_ftr(cf="Complexes.tsv", verbose=False) -> Dict[str, set]:
    """
    Read the complexes file, and create a dict of complex->feature_ids
    :param cf: the Complexes file
    :return: a dict of complexes -> feature_ids
    """

    cpx2ftr = {}
    log_and_message(f"Reading complexes from PyFBA.Biochemistry.ModelSEEDDatabase.Biochemistry.{cf}", stderr=verbose)
    complexf = open_text("PyFBA.Biochemistry.ModelSEEDDatabase.Annotations", cf)
    for li in complexf:
        if li.startswith('id'):
            continue
        p = li.strip().split("\t")
        cpx2ftr[p[0]] = set()
        if 'null' == p[5]:
            continue
        for ftr in p[5].split('|'):
            cpx2ftr[p[0]].add(ftr.split(';')[0])

    complexf.close()
    return cpx2ftr


def enzymes(organism_type=None, verbose=False) -> Dict[str, PyFBA.metabolism.Enzyme]:
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

    log_and_message(f"Creating enzymes with complexes and reactions", stderr=verbose)
    # Set up enzymes with complexes and reactions
    c2f = complex_to_ftr()
    f2r = ftr_to_roles()
    for cmplx in c2f:
        if cmplx in modelseedstore.enzymes:
            if verbose:
                log_and_message(f"Warning: have duplicate {cmplx} complexes that maybe in more than once. "
                                f"Skipped later incantations", stderr=verbose)
            continue
        modelseedstore.enzymes[cmplx] = PyFBA.metabolism.Enzyme(cmplx)
        for ft in c2f[cmplx]:
            if ft in f2r:
                modelseedstore.enzymes[cmplx].add_roles({f2r[ft]})
            else:
                log_and_message(f"Warning: No functional role for {ft}", stderr=verbose)
            for ecno in re.findall(r'[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+', f2r[ft]):
                modelseedstore.enzymes[cmplx].add_ec(ecno)
    for r in rcts:
        for c in rcts[r].enzymes:
            if c in modelseedstore.enzymes:
                modelseedstore.enzymes[c].add_reaction(rcts[r].id)

    return modelseedstore.enzymes


def compounds_reactions_enzymes(organism_type=None, verbose=False) -> (Dict[str, PyFBA.metabolism.Compound],
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


def complexes(organism_type=None, verbose=False) -> Dict[str, Set[str]]:
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


def roles(organism_type=None, verbose=False) -> Dict[str, Set[str]]:
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


def parse_model_seed_data(organism_type=None, verbose=False):
    """
    Parse the model seed data and return a ModelSeed class that contains
    all the data
    :param organism_type: limit to a type of organism
    :param verbose: more output
    :return: a ModelSeed class
    :rtype: PyFBA.model_seed.ModelData
    """

    global modelseedstore
    
    if not modelseedstore.compounds:
        compounds(verbose=verbose),
    if not modelseedstore.reactions:
        reactions(organism_type=organism_type, verbose=verbose),
    if not modelseedstore.enzymes:
        enzymes(organism_type=organism_type, verbose=verbose)
    if not modelseedstore.complexes:
        complexes(organism_type=organism_type, verbose=verbose)
    if not modelseedstore.roles:
        roles(organism_type=organism_type, verbose=verbose)
    return modelseedstore


def reset_cache():
    """
    Reset the cache of modelseed data to force reparsing it
    """

    global modelseedstore
    modelseedstore = ModelData()
    modelseedstore.reset()
