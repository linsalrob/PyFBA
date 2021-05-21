import os
import sys


import PyFBA

import argparse

parser = argparse.ArgumentParser(description='Find reactions associated with a compound')
parser.add_argument('-c', help='compound name', required=True)
args = parser.parse_args()

compounds, reactions, enzymes = PyFBA.parse.compounds_reactions_enzymes('gramnegative')

for c in compounds:
    if c.name == args.c:
        wanted = c
        break

if str(wanted) in compounds:
    print("Reactions associated with generic compound {}:".format(args.c))
    for r in compounds[str(wanted)].reactions:
        print(r)
else:
    sys.stderr.write("No compound like {} found\n".format(wanted))

intracellular = wanted
intracellular.location = 'c'
if str(intracellular) in  compounds:
    print("Reactions associated with the intracellular compound {}:".format(args.c))
    for r in compounds[str(intracellular)].reactions:
        print(r)
else:
    sys.stderr.write("No compound like {} found\n".format(intracellular))

extracellular = wanted
extracellular.location = 'e'
if str(extracellular) in compounds:
    print("Reactions associated with the extracellular compound {}:".format(args.c))
    for r in compounds[str(extracellular)].reactions:
        print(r)
else:
    sys.stderr.write("No compound like {} found\n".format(extracellular))


