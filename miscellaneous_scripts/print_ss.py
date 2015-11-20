import os
import sys

from gapfill import subsystem
from random import shuffle
from parse import model_seed

__author__ = 'Rob Edwards'

compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()
#r2r = reactions.keys()
#shuffle(r2r)
# r2r = r2r[0:100]

r2r = {'rxn00006', 'rxn00189', 'rxn00194', 'rxn00206', 'rxn00322', 'rxn00405', 'rxn05216', 'rxn10136',
       'rxn17196', 'rxn20595', 'rxn26755'}
suggestions = subsystem.suggest_reactions_from_subsystems(reactions, r2r, threshold=1, verbose=True)
