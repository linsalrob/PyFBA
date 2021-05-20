"""
Print statistics about the reactions. Note that this is used in test_model_seed_parsing.py, and so you
can run this code to update those numbers!
"""

import os
import sys
import PyFBA

import PyFBA
MODELSEED_DIR = ""
if 'ModelSEEDDatabase' in os.environ:
        MODELSEED_DIR = os.environ['ModelSEEDDatabase']
else:
    sys.stderr.write("Please ensure that you install the Model SEED Database somewhere, and set the environment " +
                     "variable ModelSEEDDatabase to point to that directory.\n" +
                     " See INSTALLATION.md for more information\n")
    sys.exit(-1)
if not MODELSEED_DIR:
    sys.stderr.write("The ModelSEEDDatabase environment variable is not set.\n")
    sys.stderr.write("Please install the ModelSEEDDatabase, set the variable, and try again")
    sys.exit(-1)
if not os.path.exists(MODELSEED_DIR):
    sys.stderr.write("The MODEL SEED directory: {} does not exist.\n".format(MODELSEED_DIR))
    sys.stderr.write("Please check your installation.\n")
    sys.exit(-1)


compounds, reactions = PyFBA.parse.model_seed.reactions()

print(f"self.assertGreaterEqual(len(compounds), {len(compounds)})")
print(f"self.assertGreaterEqual(len(reactions), {len(reactions)})")
is_transport = 0
direction = {}
for r in reactions:
    if reactions[r].is_transport:
        is_transport += 1
    direction[reactions[r].direction] = direction.get(reactions[r].direction, 0) + 1


print(f"self.assertGreaterEqual(is_transport, {is_transport})")
for k in (direction):
    print(f"self.assertGreaterEqual(direction['{k}'], {direction[k]})")