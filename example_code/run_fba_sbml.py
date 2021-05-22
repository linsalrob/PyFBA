"""
Read a model from a SBML file and run FBA using the
Gnu Linear Programming Toolkit GLPK.

This script requires that you have libsbml installed to read SBML files.
This can be installed using pip install python-libsbml-experimental
for OSX and Linux.

This is intended as a demonstration of the pieces and parts that you need
in place to get FBA running with GLPK.

To get a sbml file, I suggest that you build a model at https://kbase.us,
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


def parse_sbml_file(sbmlfile):
    """
    Parse a sbml file and return a sbml object
    :param sbmlfile: The file to parse
    :return: the parsed lbsbml object
    """
    sbml_doc = libsbml.SBMLReader().readSBML(sbmlfile)
    if sbml_doc.getNumErrors() > 0:
        sys.stderr.write("Errors occurred when reading the document\n")
        sys.exit(1)
    return sbml_doc


def parse_media_file(mediafile):
    """
    Get the media compounds present
    :param mediafile: filepath to the tab-delimited media file.
    :return: dictionary with first level of keys as compound IDs. The values are another dictionary of
    compound name, formula, and charge.
    """

    readmedia = {}
    with open(mediafile, "r") as f:
        for li, ls in enumerate(f):
            # Skip header line
            if li == 0:
                continue
            ls = ls.rstrip("\n")
            cpid, name, formula, charge = ls.split("\t")
            readmedia[cpid] = {"name": name, "formula": formula, "charge": charge}
    return readmedia


def parse_stoichiometry(mdl):
    """
    Get the stoichiometry from the models and load it into the FBA object.
    :param mdl: the parsed SBML object of the model
    :return: cpIds: a list of all compound ids; rctnids: a list of all reaction ids; sm: stoichiometric matrix;
    objective function
    """

    rctnids = []  # Reaction IDs
    obfunc = []  # Objective function

    # Grab all reaction and compound objects
    rxns = mdl.getListOfReactions()

    # Only get compounds with boundary condition set to false
    # These are not to be used in the stoichiometric matrix
    cpds = [s for s in mdl.getListOfSpecies()
            if not s.getBoundaryCondition()]
    cpids = [s.getId() for s in cpds]
    stm = {x: {} for x in cpids}

    if verbose:
        print(f"Length of 'rxns' list: {len(rxns)}")
        print(f"Length of 'cpds' list: {len(cpds)}")

    # Create the stoichiometric Dictionary
    stm, obfunc = add_stoichiometry(stm, rxns, rctnids, obfunc)

    if verbose:
        print(f"Number of compounds in stiochometric matrix: {len(stm)}")

    return cpids, rctnids, stm, obfunc


def add_stoichiometry(stm, rxns, rctnids, obfunc):
    """Add stoichiometry information to matrix"""
    for r in rxns:
        rid = r.getId()
        rctnids.append(rid)
        for sp in r.getListOfReactants():
            species_name = sp.getSpecies()
            # Be sure to skip boundary compounds
            if species_name not in stm:
                continue
            stm[species_name][rid] = -float(sp.getStoichiometry())
        for sp in r.getListOfProducts():
            species_name = sp.getSpecies()
            # Be sure to skip boundary compounds
            if species_name not in stm:
                continue
            stm[species_name][rid] = float(sp.getStoichiometry())

        # Add objective function value
        rk = r.getKineticLaw()
        coeff = float(rk.getParameter("OBJECTIVE_COEFFICIENT").getValue())
        obfunc.append(coeff)

    return stm, obfunc


def create_smatrix(cpids, rctnids, stm, obfunc):
    """Create the stiochiomatrix with the fba module"""
    # Use stoichiometric Dictionary to create S matrix
    # Preallocate size of array
    st_mat = [[0.0] * len(rctnids) for i in cpids]

    if verbose:
        print(f"Number of total compounds: {len(cpids)}")
        print(f"Number of total reactions: {len(rctnids)}")
        print(f"Stiochiomatrix dimensions before fill: {len(st_mat)} x {len(st_mat[0])}")

    # Fill in S matrix
    for cidx, c in enumerate(cpids):
        for ridx, r in enumerate(rctnids):
            try:
                val = stm[c][r]
            except KeyError:
                val = 0.0
            st_mat[cidx][ridx] = val

    if verbose:
        print(f"Stiochiomatrix dimensions after fill: fill: {len(st_mat)} x {len(st_mat[0])}")

    # Load the data
    lp.load(st_mat, cpids, rctnids)

    # Objective function
    lp.objective_coefficients(obfunc)

    return cpids, rctnids


def reaction_bounds(mdl, rctnids, rmedia):
    """
    Extract bounds for each reaction that are present in the SBML file.

    Parameters:
        mdl: SBML model
        rctnids: List of reaction IDs from the model
        rmedia: Media dictionary of compound IDs and other info. Value is None
               if a media file was not given.
    """
    # Generate list of external media compounds
    media_list = []
    if rmedia:
        media_list = [m + "_e0" for m in rmedia.keys()]
    rbounds = []
    for ridx in rctnids:
        r = mdl.getReaction(ridx)
        rk = r.getKineticLaw()
        # Check if this is an exchange reaction
        # Need to make sure compound is present in rmedia
        if "EX" in ridx and rmedia and ridx.split("EX_")[-1] not in media_list:
            lb = 0.0
        else:
            lb = float(rk.getParameter("LOWER_BOUND").getValue())
        ub = float(rk.getParameter("UPPER_BOUND").getValue())
        rbounds.append((lb, ub))

    if verbose:
        print(f"Number of reaction bounds: {len(rbounds)}")
    lp.col_bounds(rbounds)


def compound_bounds(cpids):
    """
    This is the zero flux vector.

    Parameters:
        cpids: List of compound IDs from the model
    """
    cbounds = [(0.0, 0.0) for i in cpids]
    if verbose:
        print(f"Number of compound bounds: {len(cbounds)}")
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

    compound_ids, reaction_ids, stiochio_m, objFunc = parse_stoichiometry(model)
    compound_ids, reaction_ids = create_smatrix(compound_ids, reaction_ids, stiochio_m, objFunc)

    reaction_bounds(model, reaction_ids, media)
    compound_bounds(compound_ids)

    status, value = lp.solve()
    print("Solved the FBA with status: {}".format(status))
    print("Objective value: {}".format(value))
