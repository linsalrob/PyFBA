import PyFBA

def read_media_file(mediaf):
    """
    Read a media file and return a set with the media added.
        
    Returns a set of compounds that are in the media.

    :param mediaf: The file to read
    :type mediaf: str
    :return: A set of media components
    :rtype: set of metabolism.Compound
    """

    media = set()
    with open(mediaf, 'r') as f:
        for li, l in enumerate(f):
            # skip the header line
            if li == 0:
                continue
            p = l.strip().split("\t")
            c = PyFBA.metabolism.Compound(p[1], 'e')
            media.add(c)
    
    return media


