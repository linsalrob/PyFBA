import os
import sys

from PyFBA import log_and_message

try:
    from importlib.resources import open_text, path
except ImportError:
    # this is for python<3.7
    from importlib_resources import open_text, path

import PyFBA

def media_files():
    """
    Get a list of the media names that we know about
    :return: a dict of the media names and the data repositories associated wtih those names
    """

    all_media = {}
    with path('PyFBA.Biochemistry.media', 'ArgonneLB.txt') as media_path:
        for m in os.listdir(media_path):
            if 'README.md' == m or '__init__.py' == m:
                continue
            n = m.replace('.txt', '')
            all_media[n] = m
    return all_media


def pyfba_media(media_name, verbose=False):
    """
    Parse a media file that we have provided with PyFBA. This is probably the most common method to access media.
    If you want alternative media included, just let us know!

    Returns a set of compounds that are in the media.

    :param mediaf: The name of the media to read
    :type mediaf: str
    :return: A set of media components
    :rtype: set of metabolism.Compound
    """

    all_media = media_files()
    if media_name not in media_files():
        am = "\n".join(map(str, all_media.keys()))
        log_and_message(f"Sorry, we don't know what media {media_name} is. It should be one of {am}", stderr=True)
        return set()

    media = set()
    log_and_message(f"Reading media from {all_media[media_name]}", stderr=verbose)
    with open_text("PyFBA.Biochemistry.media", all_media[media_name]) as f:
        for li, l in enumerate(f):
            # skip the header line
            if li == 0:
                continue
            p = l.strip().split("\t")
            if len(p) < 2:
                sys.stderr.write("Skipped line {} as it does not have enough columns\n".format(l.strip()))
                continue
            c = PyFBA.metabolism.CompoundWithLocation(f"Media{li:03}", p[1], 'e')
            media.add(c)

    return media


def read_media_file(mediaf):
    """
    Read a media file and return a set of compounds with the media added. If the environment variable PYFBA_MEDIA_DIR
    is set, we will look in there for mediaf if we can not find it.
        
    Returns a set of compounds that are in the media.

    :param mediaf: The file to read
    :type mediaf: str
    :return: A set of media components
    :rtype: set of metabolism.Compound
    """

    media = set()

    if not os.path.exists(mediaf):
        if 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediaf)):
            mediaf = os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediaf)
        else:
            raise IOError(f"Media file {mediaf} can not be found\nSet the environment variable PYFBA_MEDIA_DIR" +
                          " to point to a directory with all the media files")

    with open(mediaf, 'r') as f:
        for li, l in enumerate(f):
            # skip the header line
            if li == 0:
                continue
            p = l.strip().split("\t")
            if len(p) < 2:
                sys.stderr.write("Skipped line {} as it does not have enough columns\n".format(l.strip()))
                continue
            c = PyFBA.metabolism.CompoundWithLocation(f"Media{li:03}", p[1], 'e')
            media.add(c)
    
    return media
