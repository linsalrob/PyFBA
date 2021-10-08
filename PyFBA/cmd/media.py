"""
Some things to get information about the media
"""

import PyFBA
import argparse
import sys


def list_media():
    for m in PyFBA.Biochemistry.media:
        print(m)


def media_compounds():
    """
    Print the compounds in the media
    :param m: the media
    :param model_data: a model data object (not necessary)
    :param verbose: more output
    """

    parser = argparse.ArgumentParser(description='List the compounds in a media formulation')
    parser.add_argument('-m', '--media', help='the name of the media', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    model_data = PyFBA.parse.model_seed.parse_model_seed_data(verbose=args.verbose)

    media = PyFBA.parse.read_media.find_media_file(args.media, model_data, args.verbose)
    for c in media:
        print(c)