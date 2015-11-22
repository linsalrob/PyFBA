import os
import sys

__author__ = 'Rob Edwards'


def read_assigned_functions(assigned_functions_file):
    """
    Read the assigned functions file from RAST and return a data structure

    :param assigned_functions_file: The assigned functions file downloaded from RAST
    :type assigned_functions_file: str
    :return: A hash of peg and function
    :rtype: dict

    """

    if not os.path.exists(assigned_functions_file):
        raise IOError("ERROR: " + assigned_functions_file + " does not exist")

    function = {}
    with open(assigned_functions_file, 'r') as ain:
        for l in ain:
            p = l.strip().split("\t")
            function[p[0]] = p[1]
    return function
