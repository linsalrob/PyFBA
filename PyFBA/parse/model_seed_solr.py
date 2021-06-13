
"""

Please note that this is the historic version, and now we have moved to parsing the JSON
files directly. See model_data.py

A parser for the SEED biochemistry modules that are available on Github
at https://github.com/ModelSEED/ModelSEEDDatabase. We have also included
them in our repo as a submodule.

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

import PyFBA
from .model_seed import location

def template_reactions(modeltype='microbial'):
    """
    Load the template reactions to adjust the model. Returns a hash of some altered parameters for the model
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
        raise IOError(os.path.join(MODELSEED_DIR, inputfile) +
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
        compounds_file = os.path.join(MODELSEED_DIR, 'Biochemistry/compounds.master.tsv')

    try:
        with open(compounds_file, 'r') as f:
            for li, l in enumerate(f):
                if li == 0:
                    # skip the header line
                    continue
                p = l.strip().split("\t")
                c = PyFBA.metabolism.Compound(p[2], '')
                c.model_seed_id = p[0]
                c.abbreviation = p[1]
                c.formula = p[3]
                c.mw = p[4]
                # there are some compounds (like D-Glucose and Fe2+) that appear >1x in the table
                if str(c) in cpds:
                    cpds[str(c)].alternate_seed_ids.add(p[0])
                else:
                    cpds[str(c)] = c
    except IOError as e:
        sys.exit("There was an error parsing " +
                 compounds_file + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))

    return cpds


def reactions(organism_type="", rctf='Biochemistry/reactions.master.tsv', verbose=False):
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

    all_reactions = {}

    try:
        with open(os.path.join(MODELSEED_DIR, rctf), 'r') as rxnf:
            for l in rxnf:
                if l.startswith('id'):
                    # ignore the header line
                    continue
                if l.startswith("#"):
                    # ignore any comment lines
                    continue

                pieces = l.strip().split("\t")
                if len(pieces) < 20:
                    sys.stderr.write("ERROR PARSING REACTION INFO: " + l)
                    continue

                rid = pieces[0]

                rxn = pieces[6]
                for i in range(len(pieces)):
                    if pieces[i] == "none" or pieces[i] == "null":
                        pieces[i] = None

                if pieces[14]:
                    deltaG = float(pieces[14])
                else:
                    deltaG = 0.0
                if pieces[15]:
                    deltaG_error = float(pieces[15])
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

                # check and see we have a valid equation
                left = left.strip()
                right = right.strip()
                if False:
                    if left == "" or right == "":
                        if verbose:
                            sys.stderr.write("One side missing for " + rxn + " ignored\n")
                        continue

                # create a new reaction object to hold all the information ...

                r = PyFBA.metabolism.Reaction(rid)

                r.deltaG = deltaG
                r.deltaG_error = deltaG_error
                if pieces[5] != '0':
                    r.is_transport = True
                all_reactions[rid] = r

                r.direction = pieces[9]

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

                    if cmpd in cpds_by_id:
                        nc = PyFBA.metabolism.Compound(cpds_by_id[cmpd].name, loc)
                    else:
                        if verbose:
                            sys.stderr.write("ERROR: Did not find " + cmpd + " in the compounds file.\n")
                        nc = PyFBA.metabolism.Compound(cmpd, loc)

                    ncstr = str(nc)

                    if ncstr in cpds:
                        nc = copy.copy(cpds[ncstr])
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

                    if cmpd in cpds_by_id:
                        nc = PyFBA.metabolism.Compound(cpds_by_id[cmpd].name, loc)
                    else:
                        if verbose:
                            sys.stderr.write("ERROR: Did not find " + cmpd + " in the compounds file.\n")
                        nc = PyFBA.metabolism.Compound(cmpd, loc)

                    ncstr = str(nc)
                    if ncstr in cpds:
                        nc = copy.copy(cpds[ncstr])
                    nc.add_reactions({rid})
                    cpds[ncstr] = nc

                    r.add_right_compounds({nc})
                    r.set_right_compound_abundance(nc, float(q))

                    newright.append("(" + str(q) + ") " + nc.name + "[" + loc + "]")

                r.equation = " + ".join(newleft) + " <=> " + " + ".join(newright)

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


def complexes(cf="SOLRDump/TemplateReactions.tsv", verbose=False):
    """
    Connection between complexes and reactions. A complex can be
    involved in many reactions.

    In addition, many complexes are involved in one reaction, so we have
    a many:many relationship here

    Read the complex file and return a hash of the complexes where
    key is the complex id and the value is a set of reactions that the
    complex is involved in.

    You can provide an optional complexes file (cf) if you don't like
    the default!

    :param cf: An optional complexes file name
    :type cf: str
    :param verbose: Print more output
    :type verbose: bool
    :return A dict of the complexes where the key is the complex id and the value is the set of reactions
    :rtype: dict
    """

    cplxes = {}
    try:
        # io.open() to enable the encoding and errors arguments when using Python2
        # io.open() will read lines as unicode objects instead of str objects
        # In Python2, unicode objects are equivalent to Python3 str objects
        with io.open(os.path.join(MODELSEED_DIR, cf), 'r', encoding='utf-8', errors='replace') as rin:
            for l in rin:
                # If using Python2, must convert unicode object to str object
                if sys.version_info.major == 2:
                    l = l.encode('utf-8', 'replace')
                if l.startswith("#") or l.startswith('id'):
                    # ignore any comment lines
                    continue

                p = l.strip().split("\t")
                if len(p) < 30:
                    if verbose:
                        sys.stderr.write("WARNING: Malformed line in " + cf + ": " + l + "\n")
                    continue
                if p[28] == "":
                    continue
                for cmplx in p[28].split(';'):
                    if cmplx not in cplxes:
                        cplxes[cmplx] = set()
                    cplxes[cmplx].add(p[1])
    except IOError as e:
        sys.stderr.write("There was an error parsing {}\n".format(os.path.join(MODELSEED_DIR, cf)))
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
        sys.exit(-1)

    return cplxes


def roles_ec(rf="SOLRDump/ComplexRoles.tsv"):
    """
    Read the roles and EC and return a hash of the roles and EC where the id
    is the role name or EC number and the value is the set of complex IDs that
    the role is inolved in.

    One role or EC can be involved in many complexes.

    You can provide an alternate roles file (rf) if you don't like the
    default.

    :param rf: an alternate roles file
    :type rf: str
    :return: A dict of role name and complex ids that the roles is involved with
    :rtype: dict

    """
    rles_ec = {}
    try:
        with open(os.path.join(MODELSEED_DIR, rf), 'r') as rin:
            for l in rin:
                if l.startswith("#") or l.startswith('complex_id'):
                    # ignore any comment lines
                    continue
                p = l.strip().split("\t")
                if p[5] not in rles_ec:
                    rles_ec[p[5]] = set()
                rles_ec[p[5]].add(p[0])

                # Try to add EC number if it exists in role name
                for ecno in re.findall('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+', l):
                    if ecno not in rles_ec:
                        rles_ec[ecno] = set()
                    rles_ec[ecno].add(p[0])
    except IOError as e:
        sys.exit("There was an error parsing " + rf + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))

    return rles_ec


def roles(rf="SOLRDump/ComplexRoles.tsv"):
    """
    Read the roles and return a hash of the roles where the id is the
    role name and the value is the set of complex IDs that the role is
    inolved in.

    One role can be involved in many complexes.

    You can provide an alternate roles file (rf) if you don't like the
    default.

    :param rf: an alternate roles file
    :type rf: str
    :return: A dict of role name and complex ids that the roles is involved with
    :rtype: dict

    """
    rles = {}
    try:
        with open(os.path.join(MODELSEED_DIR, rf), 'r') as rin:
            for l in rin:
                if l.startswith("#") or l.startswith('complex_id'):
                    # ignore any comment lines
                    continue
                p = l.strip().split("\t")
                if p[5] not in rles:
                    rles[p[5]] = set()
                rles[p[5]].add(p[0])
    except IOError as e:
        sys.exit("There was an error parsing " + rf + "\n" + "I/O error({0}): {1}".format(e.errno, e.strerror))

    return rles


def enzymes(verbose=False):
    """
    Convert each of the roles and complexes into a set of enzymes, and
    connect them to reactions.

    Return just the enzyme objects.

    You probably want to use compounds_reactions_enzymes, this is partly here
    as a test case to make sure that enzymes and complexes play well
    together

    :param verbose: Print more output
    :type verbose: bool
    :return: A dict of with complex id as key and reaction id as value
    """

    roleset = roles()
    cmplxset = complexes()
    enzs = {}
    cpds, rcts = reactions()

    # for roles the key is the role name and the value is the complex it
    # is in
    for rolename in roleset:
        # what complex is this involved in
        for complexid in roleset[rolename]:
            if complexid not in cmplxset:
                if verbose:
                    sys.stderr.write("WARNING: " + complexid + " is not in the complexes\n")
                continue

            if complexid not in enzs:
                enzs[complexid] = PyFBA.metabolism.Enzyme(complexid)
            enzs[complexid].add_roles({rolename})
            for ecno in re.findall('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+', rolename):
                enzs[complexid].add_ec(ecno)

    for complexid in cmplxset:
        if complexid not in enzs:
            if verbose:
                sys.stderr.write("WARNING: No roles found that are part" + " of complex " + complexid + "\n")
            continue
        for reactid in cmplxset[complexid]:
            if reactid in rcts:
                enzs[complexid].add_reaction(reactid)
                rcts[reactid].add_enzymes({complexid})

    return enzs


def compounds_reactions_enzymes(organism_type='', verbose=False):
    """
    Convert each of the roles and complexes into a set of enzymes, and
    connect them to reactions.

    We return three dicts, the compounds, the enzymes, and the reactions. See the individual methods for the dicts
    that we return!

    :param organism_type: The type of organism, eg. Microbial, Gram_positive, Gram_negative
    :type organism_type:str
    :param verbose:Print more output
    :type verbose:bool
    :return: The compounds, the reactions, and the enzymes in that order
    :rtype: dict of Compound, dict of Reaction, dict of Enzyme

    """

    roleset = roles()
    cmplxset = complexes()
    cpds, rcts = reactions(organism_type, verbose=verbose)
    enzs = {}

    # for roles the key is the role name and the value is the complex it
    # is in
    for rolename in roleset:
        # what complex is this involved in
        for complexid in roleset[rolename]:
            if complexid not in cmplxset:
                if verbose:
                    sys.stderr.write("WARNING: " + complexid + " is not in the complexes\n")
                continue

            if complexid not in enzs:
                enzs[complexid] = PyFBA.metabolism.Enzyme(complexid)
            enzs[complexid].add_roles({rolename})
            for ecno in re.findall('[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+', rolename):
                enzs[complexid].add_ec(ecno)

    for complexid in cmplxset:
        if complexid not in enzs:
            if verbose:
                sys.stderr.write("WARNING: No roles found that are part" + " of complex " + complexid + "\n")
            continue
        for reactid in cmplxset[complexid]:
            if reactid in rcts:
                enzs[complexid].add_reaction(reactid)
                rcts[reactid].add_enzymes({complexid})

    return cpds, rcts, enzs




