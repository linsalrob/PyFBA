from __future__ import print_function, absolute_import, division


def plate_to_media(main_source, compound):
    """
    Convert plate growth condition to a PyFBA media file

    :param main_source: The compound source being targeted (e.g. Carbon)
    :type main_source: str
    :param compound: The compound used in the well
    :type compound: str
    :return: PyFBA media file name
    :rtype: str
    """

    with open("pmanalyzer_to_pyfba.tsv", "r") as f:
        for l in f:
            ms, c, name = l.rstrip().split("\t")
            if (ms.lower() == main_source.lower() and
                    c.lower() == compound.lower()):
                return name
    return None


def plate_to_media_map():
    """
    Get mapping from plate growth conditions to PyFBA media file names
    :return: Mapping of (main source, compound) as key and value is media file
    :rtype: dict
    """

    media_map = dict()
    with open("pmanalyzer_to_pyfba.tsv", "r") as f:
        for l in f:
            ms, c, name = l.rstrip().split("\t")
            media_map[(ms, c)] = name

    return media_map
