import argparse
import os

parser = argparse.ArgumentParser(description="Convert a .gff file from rast into its assigned functions.")
parser.add_argument('-i', help='gff input file', required=True)
parser.add_argument('-o', help='output file')
parser.add_argument('-v', help='verbose output', action='store_true')
args = parser.parse_args()

if not os.path.exists(args.i):
    sys.exit(".gff file {} was not found".format(args.i))

out = "{}.assigned_functions".format(args.i[:-4])
if args.o:
    out = args.o

try:
    with open(args.i, 'r') as infile:
        with open(out, 'w') as outfile:
            header = infile.readline()
            for line in infile:
                line = line.strip().split("\t")
                data = line[8].split(";")
                fid, name = data[0], data[1]
                outfile.write("{}\t{}\n".format(fid[3:], name[5:]))

except:
    print('Could not read .gff file; please make sure your file is correct.')
