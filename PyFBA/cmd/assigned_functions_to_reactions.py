import os

import PyFBA
import argparse
import sys
from PyFBA import log_and_message


def roles_to_reactions_to_run(roles, orgtype='gramnegative', verbose=False):
    roles_to_reactions = PyFBA.filters.roles_to_reactions(roles, organism_type=orgtype, verbose=verbose)
    reactions_to_run = set()
    for role in roles_to_reactions:
        reactions_to_run.update(roles_to_reactions[role])
    log_and_message(f"There are {len(reactions_to_run)} unique reactions associated with this genome", stderr=verbose)
    return reactions_to_run


def to_reactions():
    """
    Parse the arguments.
    """

    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Convert a set of functions or roles to a list of reactions')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--roles', help='A list of functional roles in this genome, one per line')
    group.add_argument('-a', '--assigned_functions',
                       help='RAST assigned functions role (tab separated PEG/Functional Role)')
    group.add_argument('-f', '--features', help='PATRIC features.txt file (with 5 columns)')
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)
    PyFBA.parse.model_seed.parse_model_seed_data(args.type)
    if args.roles:
        if not os.path.exists(args.roles):
            sys.stderr.write(f"FATAL: {args.roles} does not exist. Please check your files\n")
            sys.exit(1)
        log_and_message(f"Getting the roles from {args.roles}", stderr=args.verbose)
        roles = PyFBA.parse.read_functional_roles(args.roles, args.verbose)
    elif args.assigned_functions:
        if not os.path.exists(args.assigned_functions):
            sys.stderr.write(f"FATAL: {args.assigned_functions} does not exist. Please check your files\n")
            sys.exit(1)
        log_and_message(f"Getting the roles from {args.assigned_functions}", stderr=args.verbose)
        roles = PyFBA.parse.assigned_functions_set(args.assigned_functions)
    elif args.features:
        if not os.path.exists(args.features):
            sys.stderr.write(f"FATAL: {args.features} does not exist. Please check your files\n")
            sys.exit(1)
        log_and_message(f"Getting the roles from {args.features}", stderr=args.verbose)
        roles = PyFBA.parse.read_features_file(args.features, args.verbose)
    else:
        sys.stderr.write("FATAL. Either a roles or functions file must be provided")
        sys.exit(1)

    reactions_to_run = roles_to_reactions_to_run(roles, args.type, args.verbose)
    with open(args.output, 'w') as out:
        for r in reactions_to_run:
            out.write(f"{r}\n")
