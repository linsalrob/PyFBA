import os
import sys

__author__ = 'Rob Edwards'
head, tail = os.path.split(__file__)
print("{} AND {}".format(head, tail))