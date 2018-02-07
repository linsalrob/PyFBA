from __future__ import print_function, absolute_import, division
from os.path import isfile, join, dirname
import sys


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

    map_file = join(dirname(__file__), "pmanalyzer_to_pyfba.tsv")
    with open(map_file, "r") as f:
        for l in f:
            ms, c, name = l.rstrip().split("\t")
            if (ms.lower() == main_source.lower() and
                    c.lower() == compound.lower()):
                return name
    raise ValueError("Media not found")


def plate_to_media_map():
    """
    Get mapping from plate growth conditions to PyFBA media file names

    :return: Mapping of (main source, compound) as key and value is media file
    :rtype: dict
    """

    media_map = dict()
    map_file = join(dirname(__file__), "pmanalyzer_to_pyfba.tsv")
    with open(map_file, "r") as f:
        for l in f:
            ms, c, name = l.rstrip().split("\t")
            media_map[(ms, c)] = name

    return media_map


def growth_plate_results(pm_file):
    """
    Build mapping of growth results for each media condition

    :param pm_file: PMAnalyzer parameter results file
    :return: Growth results for each media condition for each sample-rep
    :rtype: dict
    """
    # Check that file exists
    if not isfile(pm_file):
        raise FileNotFoundError

    results = dict()
    media_map = plate_to_media_map()  # PyFBA media mapping from PM plates
    #  Iterate through file
    with open(pm_file, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        # Determine if we have a replicate or not
        rep = "rep" in header

        # Find growth class column
        try:
            growth_col = header.index("growthclass")
            ms_col = header.index("mainsource")
            cp_col = header.index("compound")
        except ValueError:
            raise ValueError("growthclass column not found in file")

        # Build data set
        for l in f:
            contents = l.rstrip("\n").split("\t")

            # Name is a tuple with either just the sample name or with rep ID
            name = (contents[0], contents[1]) if rep else (contents[0],)

            # Get media info
            mainsource = contents[ms_col]
            compound = contents[cp_col]
            growth = "+" in contents[growth_col]

            try:
                pyfba_media = media_map[(mainsource, compound)]
            except KeyError:
                print("Media not found: ({}, {}). Skipping.".format(
                    mainsource, compound))
                continue

            # Store growth data
            if name not in results:
                results[name] = {}
            results[name][pyfba_media] = growth

    return results
