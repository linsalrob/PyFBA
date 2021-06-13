#!/usr/bin/python2.7
from __future__ import print_function
import sys
import os
import copy
import argparse


def outputReactions(fn, sfx, rxns, db_rxns, db_enz, gf=None, v=False):
    fn += "_"
    fh = open(fn + sfx, "w")
    if gf:
        fh.write("reaction\tcomplex\tfunction\tgapfilled\tgfstep\tequation\n")
    else:
        fh.write("reaction\tcomplex\tfunction\tequation\n")
    for r in rxns:
        myEnz = db_rxns[r].enzymes
        eqn = db_rxns[r].equation
        if not gf:
            currGF = ""
            gfstep = ""
        elif r in gf:
            currGF = "yes\t"
            gfstep = db_rxns[r].gapfill_method + "\t"
        else:
            currGF = "no\t"
            gfstep = "\t"
        if len(myEnz) == 0:
            if v:
                print("No enzymes found for reaction", r, file=sys.stderr)
            fh.write("{}\tnull\tnull\t{}{}{}\n".format(r, currGF, gfstep, eqn))
            continue
        for e in myEnz:
            if e not in db_enz:
                if v:
                    print(e, "does not exist in 'enzymes'", file=sys.stderr)
                fh.write("{}\t{}\tnull\t{}{}{}\n".format(r, e, currGF,
                                                         gfstep, eqn))
                continue
            myRoles = db_enz[e].roles
            if len(myRoles) == 0:
                if v:
                    print("No roles found for enzyme", e, file=sys.stderr)
                fh.write("{}\t{}\tnull\t{}{}{}\n".format(r, e, currGF,
                                                         gfstep, eqn))
                continue
            for role in myRoles:
                fh.write("{}\t{}\t{}\t{}{}{}\n".format(r, e, role, currGF,
                                                       gfstep, eqn))
    fh.close()


def outputFlux(fn, sfx):
    fh = open(fn + "_reactions_" + sfx, "w")
    for rxn, val in PyFBA.lp.col_primal_hash().items():
        fh.write(rxn + "\t" + str(val) + "\n")
    fh.close()


parser = argparse.ArgumentParser(description="Build model from roles then gap-fill model")
parser.add_argument("functions", help="Assigned functions file")
parser.add_argument("cgfunctions", help="Closest genomes functions file")
parser.add_argument("media", help="Media file")
parser.add_argument("-o", "--outsuffix", help="Suffix for output files")
parser.add_argument("--draft", help="Output starting reactions",
                    action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose stderr output",
                    action="store_true")
parser.add_argument("--dev", help="Use PyFBA dev code",
                    action="store_true")
args = parser.parse_args()

outsfx = args.outsuffix if args.outsuffix else "out"

if args.dev:
    # Import PyFBA from absoulte path
    sys.path.insert(0, os.path.expanduser("~") + "/Projects/PyFBA/")
    sys.path.insert(0, os.path.expanduser("~") + "/PyFBA/")
    print("IN DEV MODE", file=sys.stderr)
import PyFBA

# Load ModelSEED database
modeldata = PyFBA.parse.model_seed.parse_model_seed_data('gramnegative', verbose=True)

# Read in assigned functions file
assigned_functions = PyFBA.parse.read_assigned_functions(args.functions)
roles = set([i[0] for i in [list(j) for j in assigned_functions.values()]])
print("There are {} unique roles in this genome".format(len(roles)),
      file=sys.stderr)

# Obtain dictionary of roles and their reactions
#roles_to_reactions = PyFBA.filters.roles_to_reactions(roles)
#reactions_to_run = set()
#for role in roles_to_reactions:
#    reactions_to_run.update(roles_to_reactions[role])
#print("There are {}".format(len(reactions_to_run)),
#      "unique reactions associated with this genome", file=sys.stderr)

# Obtain enzyme complexes from roles
complexes = PyFBA.filters.roles_to_complexes(roles)
if args.verbose:
    print("There are", len(complexes["complete"]), "complete and",
          len(complexes["incomplete"]), "incomplete enzyme complexes",
          file=sys.stderr)
# Get reactions from only completed complexes
reactions_to_run = set()
for c in complexes["complete"]:
    reactions_to_run.update(modeldata.enzymes[c].reactions)
print("There are {}".format(len(reactions_to_run)),
      "unique reactions associated with this genome", file=sys.stderr)

# Remove reactions IDs that do not not have a reaction equation associated
tempset = set()
for r in reactions_to_run:
    if r in modeldata.reactions:
        tempset.add(r)
    elif args.verbose:
        print("Reaction ID {}".format(r),
              "is not in our reactions list. Skipped",
              file=sys.stderr)
reactions_to_run = tempset
if args.draft:
    outputReactions("origreactions", outsfx,
                    reactions_to_run, modeldata.reactions, modeldata.enzymes, gf=None, v=args.verbose)

# Load our media
media = PyFBA.parse.read_media_file(args.media)
print("Our media has {} components".format(len(media)), file=sys.stderr)

# Define a biomass equation
biomass_equation = PyFBA.metabolism.biomass_equation('gramnegative')


# Run FBA on our media
status, value, growth = PyFBA.fba.run_fba(modeldata,
                                          reactions_to_run,
                                          media,
                                          biomass_equation)
print("Initial run has a biomass flux value of",
      "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# Check to see if model needs any gap-filling
if growth:
    print("Model grew without gap-filling", file=sys.stderr)
    sys.exit()

# Gap-fill the model
added_reactions = []
original_reactions_to_run = copy.copy(reactions_to_run)

# Media import reactions
if not growth:
    print("Adding media import reactions", file=sys.stderr)
    media_reactions = PyFBA.gapfill.suggest_from_media(modeldata,
                                                       reactions_to_run, media)
    added_reactions.append(("media", media_reactions))
    reactions_to_run.update(media_reactions)
    print("Attempting to add", len(media_reactions), "reacitons",
          file=sys.stderr)
    print("Total # reactions:", len(reactions_to_run), file=sys.stderr)

    status, value, growth = PyFBA.fba.run_fba(modeldata,
                                              reactions_to_run,
                                              media,
                                              biomass_equation)
    print("The biomass reaction has a flux of",
          "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# Essential reactions
if not growth:
    print("Adding essential reactions", file=sys.stderr)
    essential_reactions = PyFBA.gapfill.suggest_essential_reactions()
    added_reactions.append(("essential", essential_reactions))
    reactions_to_run.update(essential_reactions)
    print("Attempting to add", len(essential_reactions), "reacitons",
          file=sys.stderr)
    print("Total # reactions:", len(reactions_to_run), file=sys.stderr)

    status, value, growth = PyFBA.fba.run_fba(modeldata,
                                              reactions_to_run,
                                              media,
                                              biomass_equation)
    print("The biomass reaction has a flux of",
          "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# Reactions from closely related organisms
if not growth:
    print("Adding close organisms reactions", file=sys.stderr)
    reactions_from_other_orgs =\
        PyFBA.gapfill.suggest_from_roles(args.cgfunctions, modeldata.reactions)
    added_reactions.append(("close genomes", reactions_from_other_orgs))
    reactions_to_run.update(reactions_from_other_orgs)
    print("Attempting to add", len(reactions_from_other_orgs), "reacitons",
          file=sys.stderr)
    print("Total # reactions:", len(reactions_to_run), file=sys.stderr)

    status, value, growth = PyFBA.fba.run_fba(modeldata,
                                              reactions_to_run,
                                              media,
                                              biomass_equation)
    print("The biomass reaction has a flux of",
          "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# Subsystems
if not growth:
    print("Adding subsystem reactions", file=sys.stderr)
    subsystem_reactions =\
        PyFBA.gapfill.suggest_reactions_from_subsystems(modeldata.reactions,
                                                        reactions_to_run,
                                                        threshold=0.5)
    added_reactions.append(("subsystems", subsystem_reactions))
    reactions_to_run.update(subsystem_reactions)
    print("Attempting to add", len(subsystem_reactions), "reacitons",
          file=sys.stderr)
    print("Total # reactions:", len(reactions_to_run), file=sys.stderr)

    status, value, growth = PyFBA.fba.run_fba(modeldata,
                                              reactions_to_run,
                                              media,
                                              biomass_equation)
    print("The biomass reaction has a flux of",
          "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# EC-related reactions
if not growth:
    print("Adding EC-related reactions", file=sys.stderr)
    ec_reactions = PyFBA.gapfill.suggest_reactions_using_ec(roles,
                                                            modeldata.reactions,
                                                            reactions_to_run)
    added_reactions.append(("ec", ec_reactions))
    reactions_to_run.update(ec_reactions)
    print("Attempting to add", len(ec_reactions), "reacitons",
          file=sys.stderr)
    print("Total # reactions:", len(reactions_to_run), file=sys.stderr)

    status, value, growth = PyFBA.fba.run_fba(modeldata,
                                              reactions_to_run,
                                              media,
                                              biomass_equation)
    print("The biomass reaction has a flux of",
          "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# Compound-probability-based reactions
if not growth:
    print("Adding compound-probability-based reactions", file=sys.stderr)
    probable_reactions = PyFBA.gapfill.compound_probability(modeldata.reactions,
                                                            reactions_to_run,
                                                            cutoff=0,
                                                            rxn_with_proteins=True)
    added_reactions.append(("probable", probable_reactions))
    reactions_to_run.update(probable_reactions)
    print("Attempting to add", len(probable_reactions), "reacitons",
          file=sys.stderr)
    print("Total # reactions:", len(reactions_to_run), file=sys.stderr)

    status, value, growth = PyFBA.fba.run_fba(modeldata,
                                              reactions_to_run,
                                              media,
                                              biomass_equation)
    print("The biomass reaction has a flux of",
          "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# Orphan compounds
if not growth:
    print("Adding orphan-compound reactions", file=sys.stderr)
    orphan_reactions =\
        PyFBA.gapfill.suggest_by_compound(modeldata,
                                          reactions_to_run,
                                          max_reactions=1)
    added_reactions.append(("orphans", orphan_reactions))
    reactions_to_run.update(orphan_reactions)
    print("Attempting to add", len(orphan_reactions), "reacitons",
          file=sys.stderr)
    print("Total # reactions:", len(reactions_to_run), file=sys.stderr)

    status, value, growth = PyFBA.fba.run_fba(modeldata,
                                              reactions_to_run,
                                              media,
                                              biomass_equation)
    print("The biomass reaction has a flux of",
          "{} --> Growth: {}".format(value, growth), file=sys.stderr)

if not growth:
    print("UNABLE TO GAP-FILL MODEL", file=sys.stderr)
    sys.exit()

# Trimming the model
reqd_additional = set()

# Begin loop through all gap-filled reactions
while added_reactions:
    ori = copy.copy(original_reactions_to_run)
    ori.update(reqd_additional)
    # Test next set of gap-filled reactions
    # Each set is based on a method described above
    how, new = added_reactions.pop()

    # Get all the other gap-filled reactions we need to add
    for tple in added_reactions:
        ori.update(tple[1])

    # Use minimization function to determine the minimal
    # set of gap-filled reactions from the current method
    new_essential =\
        PyFBA.gapfill.minimize_additional_reactions(ori,
                                                    new,
                                                    modeldata,
                                                    media,
                                                    biomass_equation)
    # Record the method used to determine
    # how the reaction was gap-filled
    for new_r in new_essential:
        modeldata.reactions[new_r].is_gapfilled = True
        modeldata.reactions[new_r].gapfill_method = how
    reqd_additional.update(new_essential)

# Combine old and new reactions
all_reactions = original_reactions_to_run.union(reqd_additional)


status, value, growth = PyFBA.fba.run_fba(modeldata, all_reactions,
                                          media, biomass_equation)
print("The biomass reaction has a flux of",
      "{} --> Growth: {}".format(value, growth), file=sys.stderr)

# Save flux values
outputFlux("flux", outsfx)

# Save all reactions
outputReactions("allreactions", outsfx, all_reactions, modeldata.reactions, modeldata.enzymes, reqd_additional,
                args.verbose)