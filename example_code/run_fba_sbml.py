"""
Read a model from a SBML file and run FBA using the
Gnu Linear Programming Toolkit GLPK.

This script requires that you have libsbml installed to read SBML files.
This can be installed using pip install python-libsbml-experimental
for OSX and Linux.

This is intended as a demonstration of the pieces and parts that you need
in place to get FBA running with GLPK.

To get a sbml file, I suggest that you build a model at http://kbase.us,
and run it there, and then you can download the model and use this code
to run it. We have also provided sample data sets with this code in the
models/Citrobacter folder.

Any comments or questions: Rob Edwards and Daniel Cuevas

This is public domain software. It has no warranty.

Copyright Rob Edwards, Daniel Cuevas 2015

"""

import sys
import libsbml
import PyFBA.lp as lp
import argparse


def parse_sbml_file(sbmlf):
    """Parse a sbml file and return a sbml object"""
    doc = libsbml.SBMLReader().readSBML(sbmlf)
    if doc.getNumErrors() > 0:
        sys.stderr.write("Errors occurred when reading the document\n")
        sys.exit(1)

    return doc


def parse_media_file(mediaf):
    """
    Get the present media compounds.

    Parameters:
        mediaf: filepath to the tab-delimited media file.

    Returns:
        media: dictionary with first level of keys as compound IDs. The values
               are another dictionary of compound name, formula, and charge.
    """
    media = {}
    with open(mediaf, "r") as f:
        for li, l in enumerate(f):
            # Skip header line
            if li == 0:
                continue
            l = l.rstrip("\n")
            cpdId, name, formula, charge = l.split("\t")
            media[cpdId] = {"name": name,
                            "formula": formula,
                            "charge": charge}
    return media


def parse_stoichiometry(model):
    """
    Get the stoichiometry from the models and load it into the
    FBA object.

    Parameters:
        model: the parsed SBML object of the model

    Returns:
        cpIds: a list of all compound ids
        rxnIds: a list of all reaction ids
        sm: stoichiometric dictionary
            {cpdId: {rxnId: coefficient}}
    """

    sm = {}  # Compound to reactions
    cpds = []  # Compound objects
    cpdIds = []  # Compound IDs
    rxns = []  # Reaction objects
    rxnIds = []  # Reaction IDs
    objFunc = []  # Objective function

    # Grab all reaction and compound objects
    rxns = model.getListOfReactions()

    # Only get compounds with boundary condition set to false
    # These are not to be used in the stoichiometric matrix
    cpds = [s for s in model.getListOfSpecies()
            if not s.getBoundaryCondition()]
    cpdIds = [s.getId() for s in cpds]
    sm = {x: {} for x in cpdIds}

    if verbose:
        print "Length of 'rxns' list: %d" % len(rxns)
        print "Length of 'cpds' list: %d" % len(cpds)

    # Create the stoichiometric Dictionary
    sm, objFunc = addStoichiometry(sm, rxns, rxnIds, objFunc)

    if verbose:
        print "Number of compounds in sm: %d" % len(sm)

    return cpdIds, rxnIds, sm, objFunc


def addStoichiometry(sm, rxns, rxnIds, objFunc):
    """Add stoichiometry information to matrix"""
    for r in rxns:
        rid = r.getId()
        rxnIds.append(rid)
        for sp in r.getListOfReactants():
            spName = sp.getSpecies()
            # Be sure to skip boundary compounds
            if spName not in sm:
                continue
            sm[spName][rid] = -float(sp.getStoichiometry())
        for sp in r.getListOfProducts():
            spName = sp.getSpecies()
            # Be sure to skip boundary compounds
            if spName not in sm:
                continue
            sm[spName][rid] = float(sp.getStoichiometry())

        # Add objective function value
        rk = r.getKineticLaw()
        coeff = float(rk.getParameter("OBJECTIVE_COEFFICIENT").getValue())
        objFunc.append(coeff)

    return sm, objFunc


def createSMatrix(cpdIds, rxnIds, sm, objFunc):
    """Create the SMatrix with the fba module"""
    # Use stoichiometric Dictionary to create S matrix
    # Preallocate size of array
    SMat = [[0.0] * len(rxnIds) for i in cpdIds]

    if verbose:
        print "Number of total compounds: %d" % len(cpdIds)
        print "Number of total reactions: %d" % len(rxnIds)
        print "SMat dimensions before fill: %d x %d" % (len(SMat),
                                                        len(SMat[0]))

    # Fill in S matrix
    for cIdx, c in enumerate(cpdIds):
        for rIdx, r in enumerate(rxnIds):
            try:
                val = sm[c][r]
            except KeyError:
                val = 0.0
            SMat[cIdx][rIdx] = val

    if verbose:
        print "SMat dimensions after fill: %d x %d" % (len(SMat), len(SMat[0]))

    # Load the data
    lp.load(SMat, cpdIds, rxnIds)

    # Objective function
    lp.objective_coefficients(objFunc)

    return cpdIds, rxnIds


def reaction_bounds(model, rxnIds, media):
    """
    Extract bounds for each reaction that are present in the SBML file.

    Parameters:
        model: SBML model
        rxnIds: List of reaction IDs from the model
        media: Media dictionary of compound IDs and other info. Value is None
               if a media file was not given.
    """
    # Generate list of external media compounds
    if media:
        mediaList = [m + "_e0" for m in media.keys()]
    rbounds = []
    for rIdx in rxnIds:
        r = model.getReaction(rIdx)
        rk = r.getKineticLaw()
        # Check if this is an exchange reaction
        # Need to make sure compound is present in media
        if "EX" in rIdx and media and rIdx.split("EX_")[-1] not in mediaList:
            lb = 0.0
        else:
            lb = float(rk.getParameter("LOWER_BOUND").getValue())
        ub = float(rk.getParameter("UPPER_BOUND").getValue())
        rbounds.append((lb, ub))

    if verbose:
        print "Number of reaction bounds: %d" % len(rbounds)
    lp.col_bounds(rbounds)


def compound_bounds(cpdIds):
    """
    This is the zero flux vector.

    Parameters:
        cpdIds: List of compound IDs from the model
    """
    cbounds = [(0.0, 0.0) for i in cpdIds]
    if verbose:
        print "Number of compound bounds: %d" % len(cbounds)
    lp.row_bounds(cbounds)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("sbml", help='SBML Model file')
    parser.add_argument('-m', '--media', help="Media input file")
    parser.add_argument('-o', '--out', help="Output file")
    parser.add_argument('-v', '--verbose', action="store_true",
                        help="Output status messages")
    args = parser.parse_args()
    sbmlf = args.sbml
    output = args.out if args.out else None
    mediaf = args.media if args.media else None
    verbose = args.verbose

    doc = parse_sbml_file(sbmlf)
    model = doc.getModel()
    if mediaf:
        media = parse_media_file(mediaf)
    else:
        media = None

    cpdIds, rxnIds, sm, objFunc = parse_stoichiometry(model)
    cpdIds, rxnIds = createSMatrix(cpdIds, rxnIds, sm, objFunc)

    reaction_bounds(model, rxnIds, media)
    compound_bounds(cpdIds)

    status, value = lp.solve()
    print("Solved the FBA with status: {}".format(status))
    print("Objective value: {}".format(value))
