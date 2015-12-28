import argparse
import os
import sys
import PyFBA


parser = argparse.ArgumentParser(description='Print information about a reaction')
parser.add_argument('-r', help='reaction(s) to print information about', action='append', required=True)
args = parser.parse_args()

compounds, reactions, enzymes = PyFBA.parse.compounds_reactions_enzymes('gramnegative')

for r in args.r:
    roles = PyFBA.filters.reactions_to_roles(r)
    if roles:
        rolestr = "; ".join(roles[r])
    else:
        rolestr = 'No roles for this reaction'
    print("{}\t{}\t{}".format(r, reactions[r].equation, rolestr))
