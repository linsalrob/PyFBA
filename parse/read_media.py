
from metabolism import Compound


def read_media_file(mediaf):
    """
    Read a media file and return a set with the media added.
        
    Returns a set of compounds that are in the media.

    """

    media = set()
    with open(mediaf, 'r') as f:
        for li, l in enumerate(f):
            # skip the header line
            if li == 0:
                continue
            p=l.strip().split("\t")
            c=Compound.Compound(p[1], 'e')
            media.add(c)
    
    return media


