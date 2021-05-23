import os
import sys
import unittest

import PyFBA

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

class TestModelSeedParsing(unittest.TestCase):

    def setUp(self):
        """
        This is run before everything else
        """
        self.assertTrue(os.path.exists(MODELSEED_DIR))

    def test_compounds(self):
        """
        Test the compounds() method in model seed parsing
        """
        PyFBA.parse.model_seed.reset_cache()
        cmps = PyFBA.parse.model_seed.compounds()
        self.assertEqual(len(cmps), 33992, 'The compounds list has changed. Most likely the model seed has been ' +
                         'updated and the test code is wrong!')

    def test_template_working(self):
        """Test the template parsing is correcting the orientation of reactions"""
        # In core rxn00148 has direction <
        # after Microbial template parsing rxn00148 has direction =

        PyFBA.parse.model_seed.reset_cache()
        reactions = PyFBA.parse.model_seed.reactions()
        self.assertEqual(reactions['rxn00148'].direction, '<')
        reactions = PyFBA.parse.model_seed.reactions('microbial')
        self.assertEqual(reactions['rxn00148'].direction, '=')

    def test_template_reactions(self):
        """
        Test parsing the template reactions in the model seed
        """

        PyFBA.parse.model_seed.reset_cache()
        enz = PyFBA.parse.model_seed.template_reactions('microbial')
        errstr = f"The microbial template enzymes are {len(enz)}"
        self.assertEqual(len(enz), 20351, errstr)
        allkeys = enz.keys()
        errstr = f"The allkeys is {len(allkeys)}"
        self.assertEqual(len(allkeys), 20351, errstr)

        self.assertIn('direction', enz[list(allkeys)[0]], "The model seed template data should contain direction")
        self.assertIn('enzymes', enz[list(allkeys)[0]], "The model seed template data should contain enzymes")

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
        PyFBA.parse.model_seed.reset_cache()
        reactions = PyFBA.parse.model_seed.reactions()

        self.assertEqual(len(reactions), 43774)
        is_transport = 0
        direction = {}
        for r in reactions:
            if reactions[r].is_transport:
                is_transport += 1
            direction[reactions[r].direction] = direction.get(reactions[r].direction, 0) + 1

        self.assertGreaterEqual(is_transport, 5000)

        self.assertEqual(len(direction), 3)
        self.assertGreaterEqual(direction['='], 27000)
        self.assertGreaterEqual(direction['>'], 10000)
        self.assertGreaterEqual(direction['<'], 2000)


    def test_enzymes(self):
        """ Test just the enzymes methods"""
        PyFBA.parse.model_seed.reset_cache()
        enzs = PyFBA.parse.model_seed.enzymes()
        self.assertEqual(len(enzs), 9423, f"THere are {len(enzs)} enzymes")

    def test_complexes(self):
        """ Test getting the complexes back"""
        cmplx = PyFBA.parse.model_seed.complexes()
        self.assertGreaterEqual(len(cmplx), 7000, f"THere are {len(cmplx)} enzymes")

    def test_roles(self):
        """ Test getting the complexes back"""
        rles = PyFBA.parse.model_seed.roles()
        self.assertGreaterEqual(len(rles), 7000, f"THere are {len(rles)} enzymes")

    def test_compounds_enzymes_and_reactions(self):
        """Test retrieving the compounds, enzymes, and reactions"""
        PyFBA.parse.model_seed.reset_cache()
        cpds, rcts, enzs = PyFBA.parse.model_seed.compounds_reactions_enzymes()
        # Getting all three gave 33992 compoounds, 43774 reactions, and 9423 enzymes
        self.assertEqual(len(enzs), 9423, f"THere are {len(enzs)} enzymes")
        self.assertEqual(len(rcts), 43774, f"There are {len(rcts)} reactions")
        self.assertEqual(len(cpds), 33992, f"There are {len(cpds)} compounds")

