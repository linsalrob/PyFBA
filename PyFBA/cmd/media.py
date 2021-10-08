"""
Some things to get information about the media
"""

import PyFBA
import argparse



def list_media():
    for m in PyFBA.Biochemistry.media:
        print(m)


def compounds_in_media(m, model_data=None, verbose=False):
    """
    Print the compounds in the media
    :param m: the media
    :param model_data: a model data object (not necessary)
    :param verbose: more output
    """

    if not model_data:
        model_data = PyFBA.parse.model_seed.parse_model_seed_data(verbose=verbose)

    media = PyFBA.parse.read_media.find_media_file(m, model_data, verbose)
    for c in media:
        print(c)