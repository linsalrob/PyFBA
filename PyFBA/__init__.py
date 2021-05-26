import logging
import os
import sys
from datetime import datetime



from .version import __version__
from .logs import log_and_message, message

# define a logger
logger = logging.getLogger('PyFBA')
logger.setLevel(5)
loglocation = os.path.join(os.getcwd(), 'PyFBA.log')
message(f"We are logging to {loglocation}", "GREEN", 'stdout')
hdlr = logging.FileHandler(loglocation)
fmt = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
hdlr.setFormatter(fmt)
logger.addHandler(hdlr)

log_and_message(f"Welcome to PyFBA version {__version__} started at {datetime.now()}", stdout=True)

MODELSEED_DIR = ""
if os.path.exists('Biochemistry/ModelSEEDDatabase'):
    MODELSEED_DIR = 'Biochemistry/ModelSEEDDatabase'
elif 'ModelSEEDDatabase' in os.environ:
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

sys.stderr.write(f"We are using {MODELSEED_DIR} for our data\n")


import PyFBA.fba
import PyFBA.filters
import PyFBA.gapfill
import PyFBA.metabolism
import PyFBA.parse
import PyFBA.model
import PyFBA.model_seed



all = [
    '__version__', 'logger', 'log_and_message', 'message', 'MODELSEED_DIR'
    ]


