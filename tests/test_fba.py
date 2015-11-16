import copy
import unittest

import fba
import metabolism
import parse.model_seed


class testFBA(unittest.TestCase):
    def setUp(self):
        """
        This is run before everything else
        """

    def test_external_reactions(self):
        """Testing the fba external reactions"""
        compounds, reactions = parse.model_seed.reactions()
        # create a few external compounds
        cpds = compounds.keys()[0:10]
        model_cpds = set()
        for c in cpds:
            new_comp = copy.copy(compounds[c])
            new_comp.location = 'e'
            compounds[str(new_comp)] = new_comp
            model_cpds.add(str(new_comp))

        upsec = fba.uptake_and_secretion_reactions(model_cpds, compounds)
        self.assertEquals(len(upsec), 10)

        # now remove them and see there is nothing left
        emptyset = fba.remove_uptake_and_secretion_reactions(upsec)
        self.assertEqual(len(emptyset), 0)

    def test_create_sm(self):
        """Test the stoichiometric matrix"""

        # we need to get everything for this!
        compounds, reactions, enzymes = parse.model_seed.enzymes_and_reactions()
        self.assertIn('CoA (location: c)', compounds)
        reactions2run = reactions.keys()[0:20]
        bme = metabolism.biomass_equation('gram_negative')
        cp, rc, reactions = fba.create_stoichiometric_matrix(reactions2run, reactions, compounds, [], bme)
        self.assertEqual(len(cp), 1)
        self.assertEqual(len(rc), 2)
        self.assertEqual(len(reactions), 3)