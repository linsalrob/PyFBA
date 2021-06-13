import argparse
import os
import sys

import PyFBA

parser = argparse.ArgumentParser(description="Print the equations of reactions in an sbml file")
parser.add_argument('-s', help='sbml file', required=True)
parser.add_argument('-m', help='media file', required=True)
parser.add_argument('-v', help='verbose output', action='store_true')
args = parser.parse_args()

if not os.path.exists(args.s):
    sys.exit("SBML file {} was not found".format(args.s))

sbml = PyFBA.parse.parse_sbml_file(args.s, False)
# create modeldata objects
modeldata = PyFBA.model_seed.ModelData(compounds=sbml.compounds, reactions=sbml.reactions)
modeldata.rebuild_indices()

# get the reactions from the SBML file. This is a dict of metabolite.Reaction objects
rxn = sbml.reactions
if args.v:
    sys.stderr.write("List of 'rxns' list: {}\n".format(len(rxn)))

# filter out reactions where one of the components is an uptake/secretion reaction (aka boundary reaction)
# The uptake and secretion compounds typically have reaction bounds that allow them to be consumed
# (i.e. diffuse away from the cell) but not produced. However, our media components can also increase in
# concentration (i.e. diffuse to the cell) and thus the bounds are set higher. Whenever you change the growth media,
# you also need to adjust the reaction bounds to ensure that the media can be consumed!
reactions_to_run = set()
uptake_secretion_reactions = {}
biomass_equation = None
for r in rxn:
    if 'biomass_equation' in rxn[r].name.lower():
        biomass_equation = rxn[r]
        continue
    is_boundary = False
    for c in rxn[r].all_compounds():
        if c.uptake_secretion:
            is_boundary = True
    if is_boundary:
        rxn[r].is_uptake_secretion = True
        uptake_secretion_reactions[r] = rxn[r]
    else:
        reactions_to_run.add(r)

# get all the compounds from the SBML file. This is a dict of metabolite.Compound objects
cps = sbml.compounds
# filter for compounds that are boundary compounds
cpds = {}
for c in cps:
    if not cps[c].uptake_secretion:
        cpds[c] = cps[c]

if args.v:
    sys.stderr.write("List of 'cpds' list: {}\n".format(len(cpds)))

# read the media file
media = PyFBA.parse.read_media_file(args.m)
# correct the names
media = PyFBA.parse.correct_media_names(media, cps)

# adjust the lower bounds of uptake secretion reactions for things that are not in the media
for u in uptake_secretion_reactions:
    is_media_component = False
    for c in uptake_secretion_reactions[u].all_compounds():
        if c in media:
            is_media_component = True
    if not is_media_component:
        rxn[u].lower_bound = 0.0
        uptake_secretion_reactions[u].lower_bound = 0.0

# run the fba
status, value, growth = PyFBA.fba.run_fba(modeldata, reactions_to_run, media, biomass_equation,
                                          uptake_secretion_reactions, args.v)
sys.stdout.write("The FBA completed with value; {} and growth: {}\n".format(value, growth))
