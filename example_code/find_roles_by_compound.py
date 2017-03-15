import os
import sys

import PyFBA

import argparse

parser = argparse.ArgumentParser(description='Find roles associated with a compound')
parser.add_argument('-c', help='compound name', required=True)
args = parser.parse_args()

compounds, reactions, enzymes = PyFBA.parse.compounds_reactions_enzymes('gramnegative')

generic = PyFBA.metabolism.Compound(args.c, '')
if str(generic) in compounds:
    print("Roles associated with generic compound {}:".format(args.c))
    for r in PyFBA.filters.reactions_to_roles(compounds[str(generic)].reactions):
        print(r)
else:
    sys.stderr.write("No compound like {} found\n".format(generic))

intracellular = PyFBA.metabolism.Compound(args.c, 'c')
if str(intracellular) in compounds:
    print("Roles associated with the intracellular compound {}:".format(args.c))
    roles = PyFBA.filters.reactions_to_roles(compounds[str(intracellular)].reactions)
    for r in roles:
        print(r + "\t" + str(roles[r]))
else:
    sys.stderr.write("No compound like {} found\n".format(intracellular))

extracellular = PyFBA.metabolism.Compound(args.c, 'e')
if str(extracellular) in compounds:
    print("Roles associated with the extracellular compound {}:".format(args.c))
    for r in PyFBA.filters.reactions_to_roles(compounds[str(extracellular)].reactions):
        print(r)
else:
    sys.stderr.write("No compound like {} found\n".format(extracellular))


