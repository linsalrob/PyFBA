"""
Generate some statistics about the data. You can use this to adjust the
numbers in the test harness.
"""


import os
import sys
import PyFBA
from timeit import default_timer as timer

for i in range(5):
    print(f"\nIteration {i}")
    start = timer()
    compounds = PyFBA.parse.model_seed.compounds()
    end = timer()
    print(f"Initially there are {len(compounds)} compounds.", end=" ")
    print(f"This search took {end-start} seconds")
    
    start = timer()
    reactions = PyFBA.parse.model_seed.reactions()
    end = timer()
    print(f"Initially there are {len(reactions)} reactions.", end=" ")
    print(f"This search took {end-start} seconds")
    compounds = PyFBA.parse.model_seed.compounds()
    print(f"After adding the reactions, there are {len(compounds)} compounds.")

    start = timer()
    enzymes = PyFBA.parse.model_seed.enzymes()
    end = timer()
    print(f"Initially there are {len(enzymes)} enzymes.", end=" ")
    print(f"This search took {end-start} seconds")
    reactions = PyFBA.parse.model_seed.reactions()
    compounds = PyFBA.parse.model_seed.compounds()
    print(f"After adding the enzymes, there are {len(compounds)} compounds.")
    print(f"After adding the enzymes, there are {len(reactions)} reactions.")

    start = timer()
    compounds, reactions, enzymes =  PyFBA.parse.model_seed.compounds_reactions_enzymes()
    end = timer()
    print(f"Getting all three gave {len(compounds)} compoounds,", end =" ")
    print(f"{len(reactions)} reactions, and {len(enzymes)} enzymes.")
    print(f"Getting all three  took {end-start} seconds")
    PyFBA.parse.model_seed.reset_cache()

