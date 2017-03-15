import os
import sys

import PyFBA


def read_media_file(mediaf):
    """
    Read a media file and return a set with the media added. If the environment variable PYFBA_MEDIA_DIR
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
            raise IOError("Media file {} can not be found\nPlease set the environment variable PYFBA_MEDIA_DIR to point to a directory with all the media files".format(mediaf))

    with open(mediaf, 'r') as f:
        for li, l in enumerate(f):
            # skip the header line
            if li == 0:
                continue
            p = l.strip().split("\t")
            if len(p) < 2:
                sys.stderr.write("Skipped line {} as it does not have enough columns\n".format(l.strip()))
                continue
            c = PyFBA.metabolism.Compound(p[1], 'e')
            media.add(c)
    
    return media


