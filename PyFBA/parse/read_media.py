import os
import sys

import PyFBA
from PyFBA import log_and_message

try:
    from importlib.resources import open_text, path
except ImportError:
    # this is for python<3.7
    from importlib_resources import open_text, path


def media_files():
    """
    Get a list of the media names that we know about
    :return: a dict of the media names and the data repositories associated wtih those names
    """
    
    return PyFBA.Biochemistry.media


def raw_media(media_name, verbose=False):
    """
    Parse a media file that we have provided with PyFBA.

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


def correct_media_names(media, cpds, verbose=False):
    """
    Correct the names in media files so they match names in the SBML files. Basically replacing '-' with '_'
    or '+' with ' '

    :param cpds: A set of compounds that are in the model
    :type cpds: set
    :param media: A set of compounds that define the media
    :type media: set
    :param verbose: more output
    :type verbose: bool
    :return: A new media set with corrected names
    :rtype: set
    """

    # correct some of the media names so that they match the compounds in the
    # SBML file. This is why we should use compound IDs and not names!
    newmedia = set()
    compounds_by_name = {c.name: c for c in cpds}
    warned_compounds = False
    for m in media:
        if m.name in compounds_by_name:
            media_component = PyFBA.metabolism.CompoundWithLocation.from_compound(compounds_by_name[m.name], 'e')
            newmedia.add(media_component)
            # log_and_message(f"Found media component by name {media_component}\n", "GREEN", stdout=True)
            continue

        testname = m.name.replace('-', '_')
        if testname in compounds_by_name:
            media_component = PyFBA.metabolism.CompoundWithLocation.from_compound(compounds_by_name[testname], 'e')
            newmedia.add(media_component)
            log_and_message(f"Found media component by corrected name (-:_) {media_component}\n", "GREEN", stderr=verbose)
            continue

        testname = m.name.replace('+', '')
        if testname in compounds_by_name:
            media_component = PyFBA.metabolism.CompoundWithLocation.from_compound(compounds_by_name[testname], 'e')
            newmedia.add(media_component)
            log_and_message(f"Found media component by corrected name (+:'') {media_component}\n", "GREEN", stderr=verbose)
            continue

        log_and_message(f"Checking media compounds: Our compounds do not include  {m.name}", stderr=verbose)
        warned_compounds = True
        newmedia.add(m)

    if warned_compounds:
        log_and_message("""
Please note: The warnings about media not being found in compounds are not fatal.
It just means that we did not find that compound anywhere in the reactions, and so it is unlikely to be
needed or used. We typically see a few of these in rich media. 
        """, stderr=verbose)
    return newmedia


def pyfba_media(media_name, modeldata, verbose=False):
    """
    Read the media file, and correct the media for the compound names in modeldata.compounds
    This is probably the most common method to access media.

    If you want alternative media included, just let us know!

    :param media_name: the media file to read
    :type media_name: str
    :param modeldata: the modeldata object
    :type mediaf: PyFBA.model_seed.ModelData
    :param verbose: more output
    :type verbose: bool
    :return: a set of media normalized to the compounds in modeldata
    """

    media = raw_media(media_name, verbose=verbose)
    return correct_media_names(media, modeldata.compounds, verbose=verbose)


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

def find_media_file(mediafile, modeldata, verbose=False):
    """
    Find a media file, read it, and return a set of compounds
    :param modeldata: the modeldata object
    :type modeldata: PyFBA.model_seed.ModelData
    :param mediafile: the media file to read
    :param verbose: more output
    :type verbose: bool
    :return: a set of media compounds
    :rtype: Set[PyFBA.metabolism.Compound]
    """

    if mediafile in media_files():
        log_and_message(f"parsing media directly from {mediafile}", stderr=verbose)
        # pyfba media already corrects the names, so we can  just return it.
        return pyfba_media(mediafile, modeldata)
    elif os.path.exists(mediafile):
        log_and_message(f"parsing media file {mediafile}", stderr=verbose)
        media = read_media_file(mediafile)
    elif 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)):
        log_and_message(f"parsing media file {os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)}",
                        stderr=verbose)
        media = read_media_file(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile))
    else:
        log_and_message(f"Can't figure out how to parse media from {mediafile}", stderr=True, loglevel="CRITICAL")
        sys.exit(-1)
    return correct_media_names(media, modeldata.compounds)
