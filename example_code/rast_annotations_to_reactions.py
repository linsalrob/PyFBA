import PyFBA

__author__ = 'Rob Edwards'

import argparse
import os
import sys

parser=argparse.ArgumentParser(description='Convert a RAST spreadsheet file to a list of roles')
parser.add_argument('-a', help='RAST spreadsheet (tab-separated text) file')
parser.add_argument('-v', help='verbose', action='store_true')
args = parser.parse_args()

if args.a.endswith('xls'):
    sys.exit("Please download the RAST spreadsheet file as a text file, not as an excel file.")

af = PyFBA.parse.read_assigned_functions(args.a)
roles = set()
[roles.update(i) for i in af.values()]
rc = PyFBA.filters.roles_to_reactions(roles)
reactions = set()
for r in rc:
    reactions.update(rc[r])
print("\n".join(reactions))
