import os
import sys

import re

__author__ = 'Rob Edwards'


def roles_of_function(role):
    """
    Separate a function into a set of roles.

    In the SEED a function associated with a PEG can have more than one role. This is denoted by joining roles
    with either ' / ' (a multi-domain, multifunctional gene), ' ; ' (an ambiguous function), or ' @ '
    (a Single domain playing multiple roles). For more information see
    http://www.nmpdr.org/FIG/Html/SEED_functions.html

    This method just takes a function and returns a set of the role(s). If there is a single role, it is still a set.

    :param role: The functional role
    :type role: str
    :return: A set of the roles
    :rtype: set
    """

    # remove comments from functions and split multiple functions
    func = re.sub('\s+[#!]\s.*$', '', role)
    return set(re.split('\s*;\s+|\s+[;/@]\s+', func))


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
            function[p[1]] = roles_of_function(p[7])
    return function


def read_assigned_functions(assigned_functions_file):
    """
    Read the assigned functions file from RAST and return a data structure. Note that a peg can be associated with >1
    function (e.g. a bifunctional protein).

    :param assigned_functions_file: The assigned functions file downloaded from RAST
    :type assigned_functions_file: str
    :return: A hash of peg and a set of the function(s) associated with that peg
    :rtype: dict of sets

    """

    if not os.path.exists(assigned_functions_file):
        raise IOError("ERROR: {} does not exist".format(assigned_functions_file))

    function = {}
    with open(assigned_functions_file, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            function[p[0]] = roles_of_function(p[1])
    return function
