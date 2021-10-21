"""
Convert a set of reactions to a list of roles
"""

import os
import sys
import argparse

import PyFBA
from PyFBA import log_and_message


def convert_reactions_to_roles():
    """
    Parse the arguments and start the gapfilling.
    """

    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Get the roles associated with a file of reactions')
    parser.add_argument('-r', '--reactions', help='A list of the reactions you have, one per line')
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)
    model_data = PyFBA.parse.model_seed.parse_model_seed_data(args.type)
    rcts = set()
    with open(args.reactions, 'r') as fin:
        for li in fin:
            rcts.add(li.strip())

    roles = PyFBA.filters.reactions_to_roles(rcts, args.type, verbose=args.verbose)

    with open(args.output, 'w') as out:
        for rid in roles:
            for rl in roles[rid]:
                out.write(f"{rid}\t{rl}\n")


def convert_reactions_to_aliases():
    """
    Convert the reactions to a list of aliases
    """

    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Get the roles associated with a file of reactions')
    parser.add_argument('-r', '--reactions', help='A list of the reactions you have, one per line')
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)
    model_data = PyFBA.parse.model_seed.parse_model_seed_data(args.type)
    rcts = set()
    with open(args.reactions, 'r') as fin:
        for li in fin:
            rcts.add(li.strip())

    aliases = {}
    allall  = set()

    for r in rcts:
        if r in model_data.reactions:
            aliases[r] = {}
            if model_data.reactions[r].aliases:
                for a in model_data.reactions[r].aliases:
                    log_and_message(f"Splitting {a}")
                    k, v = a.split(": ", 1)
                    aliases[r][k] = v
                    allall.add(k)

    allaliases = sorted(allall)

    with open(args.output, 'w') as out:
        out.write("\t".join(["reaction id"] + allaliases))
        out.write("\n")
        for r in rcts:
            out.write(r)
            for a in allaliases:
                if a in aliases[r]:
                    out.write(f"\t{aliases[r][a]}")
                else:
                    out.write("\t")
            out.write("\n")



if __name__ == "__main__":
    convert_reactions_to_roles()
