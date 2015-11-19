import copy
import os
import unittest

import fba
import metabolism
import parse.model_seed
from parse import read_media_file


class TestFBA(unittest.TestCase):
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
        self.assertEqual(len(upsec), 10)

        # now remove them and see there is nothing left
        emptyset = fba.remove_uptake_and_secretion_reactions(upsec)
        self.assertEqual(len(emptyset), 0)

    def test_create_sm(self):
        """Test the stoichiometric matrix"""

        # we need to get everything for this!
        compounds, reactions, enzymes = parse.model_seed.compounds_reactions_enzymes()
        reactions2run = reactions.keys()[0:20]
        biomass_equation = metabolism.biomass_equation('gram_negative')
        cp, rc, reactions = fba.create_stoichiometric_matrix(reactions2run, reactions, compounds, set(),
                                                             biomass_equation)
        self.assertEqual(len(cp), 102)
        self.assertEqual(len(rc), 23)
        self.assertEqual(len(reactions), 34698)

    def test_run_fba(self):
        """Test running the fba. We build a run a complete FBA based on reaction_list.txt"""
        self.assertTrue(os.path.exists('tests/reaction_list.txt'))
        self.assertTrue(os.path.exists('media/ArgonneLB.txt'))
        compounds, reactions, enzymes = parse.model_seed.compounds_reactions_enzymes('gram_negative')
        reactions2run = set()
        with open('tests/reaction_list.txt', 'r') as f:
            for l in f:
                if l.startswith('#'):
                    continue
                if "biomass" in l.lower():
                    continue
                r = l.strip()
                if r in reactions:
                    reactions2run.add(r)
        media = read_media_file('media/ArgonneLB.txt')
        biomass = metabolism.biomass_equation('gram_negative')
        status, value, growth = fba.run_fba(compounds, reactions, reactions2run, media, biomass, verbose=False)
        self.assertTrue(growth)
        value = float('%0.3f' % value)
        self.assertEqual(value, 282.674)
