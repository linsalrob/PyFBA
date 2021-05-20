import os
import sys
import unittest

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

class TestModelSeedParsing(unittest.TestCase):

    def setUp(self):
        """
        This is run before everything else
        """
        self.assertTrue(os.path.exists(MODELSEED_DIR))

    def test_template_working(self):
        """Test the template parsing is correcting the orientation of reactions"""
        # In core rxn00148 has direction <
        # after Microbial template parsing rxn00148 has direction =

        compounds, reactions = PyFBA.parse.model_seed.reactions()
        self.assertEqual(reactions['rxn00148'].direction, '<')
        compounds, reactions = PyFBA.parse.model_seed.reactions('microbial')
        self.assertEqual(reactions['rxn00148'].direction, '=')

    def test_template_reactions(self):
        """
        Test parsing the template reactions in the model seed
        """

        enz = PyFBA.parse.model_seed.template_reactions('microbial')
        errstr = f"The microbial template enzymes are {len(enz)}"
        self.assertEqual(len(enz), 20351, errstr)
        allkeys = enz.keys()
        errstr = f"The allkeys is {len(allkeys)}"
        self.assertEqual(len(allkeys), 20351, errstr)

        self.assertIn('direction', enz[list(allkeys)[0]], "The model seed template data should contain direction")
        self.assertIn('enzymes', enz[list(allkeys)[0]], "The model seed template data should contain enzymes")

    def test_compounds(self):
        """
        Test the compounds() method in model seed parsing
        """

        cmps = PyFBA.parse.model_seed.compounds()
        self.assertEqual(len(cmps), 33845, 'The compounds list has changed. Most likely the model seed has been ' +
                         'updated and the test code is wrong!')

    def test_locations(self):
        """
        Test the location strings. These should be hard coded in the parser code
        """

        locs = PyFBA.parse.model_seed.location()
        self.assertEqual(len(locs), 3)
        self.assertEqual(locs['0'], 'c')
        self.assertEqual(locs['1'], 'e')
        self.assertEqual(locs['2'], 'h')

    def test_reactions(self):
        """Test parsing the reactions by parse.model_seed"""
        compounds, reactions = PyFBA.parse.model_seed.reactions()
        # in the current version of modelseeddatabase (May 2021)
        # we have the following data -
        #
        # Note that these numbers are occasionally updated, and so you may need to update the test values.
        # To mitigate this, we use >= in our comparison (in the hope none are deleted!)
        #
        # Please also see the code in example_code/reaction_statistics.py that will generate an
        # updated version of these numbers for you
        self.assertGreaterEqual(len(compounds), 56767)
        self.assertGreaterEqual(len(reactions), 43774)
        is_transport = 0
        direction = {}
        for r in reactions:
            if reactions[r].is_transport:
                is_transport += 1
            direction[reactions[r].direction] = direction.get(reactions[r].direction, 0) + 1

        self.assertGreaterEqual(is_transport, 5505)

        self.assertEqual(len(direction), 3)
        self.assertGreaterEqual(direction['='], 27676)
        self.assertGreaterEqual(direction['>'], 12767)
        self.assertGreaterEqual(direction['<'], 3331)


    def test_compounds_enzymes_and_reactions(self):
        """Test retrieving the compounds, enzymes, and reactions"""
        cpds, rcts, enzs = PyFBA.parse.model_seed.compounds_reactions_enzymes()
        self.assertEqual(len(enzs), 4067, f"THere are {len(enzs)} enzymes")
        self.assertGreaterEqual(len(rcts), 34696, f"There are {len(rcts)} reactions")
        self.assertGreaterEqual(len(cpds), 45616, f"There are {len(cpds)} compounds")

