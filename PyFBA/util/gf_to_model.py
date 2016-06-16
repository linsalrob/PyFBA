#!/usr/bin/python2.7
from __future__ import print_function
import argparse
import sys
from os import listdir
from os.path import isfile, isdir, join

parser = argparse.ArgumentParser()
parser.add_argument("gfs", help="Gap-fill results directory")
parser.add_argument("-v", "--verbose", help="Verbose stderr output",
                    action="store_true")

args = parser.parse_args()
gfdir = args.gfs

# Check if directory exists
if not isdir(gfdir):
    print("Directory '", gfdir, "' does not exist", sep="", file=sys.stderr)
    sys.exit(1)

files =  [fi for fi in listdir(gfdir) if isfile(join(gfdir, fi))]
if args.verbose:
    print(gfdir, "contains", len(files), "files", file=sys.stderr)

reactions = set()
# Iterate through gap-filled files
for f in [fi for fi in listdir(gfdir) if isfile(join(gfdir, fi))]:
    with open(join(gfdir, f), "r") as fh:
        # Header line
        next(fh)
        for l in fh:
            l = l.strip()
            ll = l.split("\t")
            # Store reaction ID
            reactions.add(ll[0])

if args.verbose:
    print(len(reactions), "unique reactions found", file=sys.stderr)

# Iterate through reaction ID set
for r in reactions:
    print(r)
