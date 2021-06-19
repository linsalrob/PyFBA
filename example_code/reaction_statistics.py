"""
Print statistics about the reactions. Note that this is used in test_model_seed_parsing.py, and so you
can run this code to update those numbers!
"""

import PyFBA
from timeit import default_timer as timer

start = timer()
modeldata = PyFBA.parse.model_seed.parse_model_seed_data("gramnegative", verbose=False)
end = timer()

print(f"self.assertGreaterEqual(len(compounds), {len(modeldata.compounds)})")
print(f"self.assertGreaterEqual(len(reactions), {len(modeldata.reactions)})")
is_transport = 0
direction = {}
for r in modeldata.reactions:
    if modeldata.reactions[r].is_transport:
        is_transport += 1
    direction[modeldata.reactions[r].direction] = direction.get(modeldata.reactions[r].direction, 0) + 1


print(f"self.assertGreaterEqual(is_transport, {is_transport})")
for k in (direction):
    print(f"self.assertGreaterEqual(direction['{k}'], {direction[k]})")


print(f"Initially there are {len(modeldata.compounds)} compounds.")
print(f"Initially there are {len(modeldata.reactions)} reactions.")
print(f"Initially there are {len(modeldata.enzymes)} enzymes.")
print(f"This search took {end-start} seconds")
reactions = PyFBA.parse.model_seed.reactions()
compounds = PyFBA.parse.model_seed.compounds()

PyFBA.parse.model_seed.reset_cache()

