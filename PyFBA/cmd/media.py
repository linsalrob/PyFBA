"""
Some things to get information about the media
"""

import PyFBA
import argparse



def list_media():
    for m in PyFBA.Biochemistry.media:
        print(m)


