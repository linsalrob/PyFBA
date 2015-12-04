import os
import sys
import unittest

import PyFBA

MODELSEED_DIR = ""
if 'ModelSEEDDatabase' in os.environ:
    MODELSEED_DIR = os.environ['ModelSEEDDatabase']

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
        # without template parsing rxn00001 has direction =
        # after GramNegative template parsing rxn00001 has direction >
        compounds, reactions = PyFBA.parse.model_seed.reactions()
        self.assertEqual(reactions['rxn00001'].direction, '=')
        compounds, reactions = PyFBA.parse.model_seed.reactions('gram_negative')
        self.assertEqual(reactions['rxn00001'].direction, '>')

    def test_template_reactions(self):
        """
        Test parsing the template reactions in the model seed
        """

        enz = PyFBA.parse.model_seed.template_reactions('microbial')
        self.assertGreaterEqual(len(enz), 19000, 'The microbial template has changed. Most likely the model seed has been ' +
                         'updated and the test code is wrong!')
        allkeys = enz.keys()
        self.assertGreaterEqual(len(allkeys), 19000)

        self.assertIn('direction', enz[allkeys[0]], "The model seed template data should contain direction")
        self.assertIn('enzymes', enz[allkeys[0]], "The model seed template data should contain enzymes")

    def test_compounds(self):
        """
        Test the compounds() method in model seed parsing
        """

        cmps = PyFBA.parse.model_seed.compounds()
        self.assertEqual(len(cmps), 27586, 'The compounds list has changed. Most likely the model seed has been ' +
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
        # in the current version of modelseeddatabase (11/16/2015)
        # we have the following data -
        #
        # Note that these numbers are occasionally updated, and so you may need to update the test values.
        # To mitigate this, we use >= in our comparison (in the hope none are deleted!)
        self.assertGreaterEqual(len(compounds), 45616)
        self.assertGreaterEqual(len(reactions), 34696)
        is_transport = 0
        direction = {}
        for r in reactions:
            if reactions[r].is_transport:
                is_transport += 1
            direction[reactions[r].direction] = direction.get(reactions[r].direction, 0) + 1

        self.assertGreaterEqual(is_transport, 5272)

        self.assertEqual(len(direction), 3)
        self.assertGreaterEqual(direction['<'], 3328)
        self.assertGreaterEqual(direction['>'], 12760)
        self.assertGreaterEqual(direction['='], 18608)

    def test_complexes(self):
        """Test parsing the complexes by parse.model_seed"""
        cmplxs = PyFBA.parse.model_seed.complexes()
        self.assertGreaterEqual(len(cmplxs), 4183)

    def test_roles(self):
        """Test the roles() method in parse.model_seed"""
        roles = PyFBA.parse.model_seed.roles()
        # this should have the same number of lines as
        #   wc -l Biochemistry/ModelSEEDDatabase/SOLRDump/ComplexRoles.tsv
        #   4747
        self.assertGreaterEqual(len(roles), 2350)

    def test_enzymes(self):
        """Test the enzymes() method in parse.model_seed"""
        enzs = PyFBA.parse.model_seed.enzymes()
        self.assertEqual(len(enzs), 4067)

    def test_compounds_enzymes_and_reactions(self):
        """Test retrieving the compounds, enzymes, and reactions"""
        cpds, rcts, enzs = PyFBA.parse.model_seed.compounds_reactions_enzymes()
        self.assertEqual(len(enzs), 4067)
        self.assertGreaterEqual(len(rcts), 34696)
        self.assertGreaterEqual(len(cpds), 45616)

