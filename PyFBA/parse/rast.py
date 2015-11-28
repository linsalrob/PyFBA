import os
import sys

__author__ = 'Rob Edwards'


def read_downloaded_data(spreadsheet_file):
    """
    Read data downloaded from RAST as a 'spreadsheet (tab-separated text format)' and return a hash of the
    proteins and functional roles in the genome.
    :param spreadsheet_file:The file downloaded from RAST
    :type spreadsheet_file: str
    :return: Dictionary of protein encoding gene identifiers and functional roles of those genes
    :rtype:dict of str and str
    """
    if not os.path.exists(spreadsheet_file):
        raise IOError("ERROR: {} does not exist".format(spreadsheet_file))

    function = {}
    with open(spreadsheet_file, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            function[p[1]] = p[7]
    return function


def read_assigned_functions(assigned_functions_file):
    """
    Read the assigned functions file from RAST and return a data structure

    :param assigned_functions_file: The assigned functions file downloaded from RAST
    :type assigned_functions_file: str
    :return: A hash of peg and function
    :rtype: dict

    """

    if not os.path.exists(assigned_functions_file):
        raise IOError("ERROR: {} does not exist".format(assigned_functions_file))

    function = {}
    with open(assigned_functions_file, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            function[p[0]] = p[1]
    return function
