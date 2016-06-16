#!/usr/bin/python2.7
from __future__ import print_function
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("params_file", help="Parameters file from PMAnalyzer")
parser.add_argument("-m", "--media_map",
                    default="/Users/dcuevas/bin/pm1_to_kbase.txt",
                    help=("Mapping file from FBA to PMAnalyzer media. "
                          "[default: bin/pm1_to_kbase.txt]"))
parser.add_argument("-t", "--threshold", type=float, default=0.40,
                    help="Growth threshold to use. [default: 0.40]")
parser.add_argument("--growth_class", action="store_true",
                    help="Use class column to identify growth")

args = parser.parse_args()
gclass = args.growth_class

media = {}
# Read in mapping file and store in dictionary
with open(args.media_map, "r") as f:
    next(f)  # Header line
    for l in f:
        l = l.strip()
        m, ms, c, w = l.split("\t")
        media[w] = m + ".txt"

# Read in growth curve parameters file
with open(args.params_file, "r") as f:
    next(f)  # Header line
    for l in f:
        l = l.strip()
        ll = l.split("\t")
        gl = float(ll[8])
        w = ll[1]
        try:
            med = media[w]
        except:
            print(w, ll[2], ll[3],
                  "doesn't have a media mapping",
                  file=sys.stderr)
            med = "unknown"
        # Check if the growth level meets threshold
        if ((gclass and "+" in ll[-1]) or
                (not gclass and gl >= args.threshold)):
            print(med)
