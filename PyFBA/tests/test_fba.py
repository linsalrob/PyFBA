import copy
import os
import unittest

import sys

import PyFBA

test_file_loc = ''
if os.path.exists('tests/roles.txt'):
    test_file_loc = 'tests'
elif os.path.exists('PyFBA/tests/roles.txt'):
    test_file_loc = 'PyFBA/tests'

media_file_loc = ''
if os.path.exists('../media/ArgonneLB.txt'):
    media_file_loc = '../media'
elif os.path.exists('media/ArgonneLB.txt'):
    media_file_loc = 'media'
elif 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], 'ArgonneLB.txt')):
    media_file_loc = os.environ['PYFBA_MEDIA_DIR']
else:
    sys.stderr.write("No media found. Can't proceed with testing the FBA.\n")
    sys.stderr.write("You can specify the media location by setting the PYFBA_MEDIA_DIR environment variable\n")

class TestFBA(unittest.TestCase):

    compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes('gram_negative')

    def setUp(self):
        """
        This is run before everything else
        """

    def test_external_reactions(self):
        """Testing the fba external reactions"""
        compounds, reactions = PyFBA.parse.model_seed.reactions()
        # create a few external compounds
        cpds = compounds.keys()[0:10]
        model_cpds = set()
        for c in cpds:
            new_comp = copy.copy(compounds[c])
            new_comp.location = 'e'
            compounds[str(new_comp)] = new_comp
            model_cpds.add(str(new_comp))

        upsec = PyFBA.fba.uptake_and_secretion_reactions(model_cpds, compounds)
        self.assertEqual(len(upsec), 10)

        # now remove them and see there is nothing left
        emptyset = PyFBA.fba.remove_uptake_and_secretion_reactions(upsec)
        self.assertEqual(len(emptyset), 0)

    def test_create_sm(self):
        """Test the stoichiometric matrix"""

        reactions2run = self.__class__.reactions.keys()[0:20]
        biomass_equation = PyFBA.metabolism.biomass_equation('gram_negative')
        cp, rc, reactions = PyFBA.fba.create_stoichiometric_matrix(reactions2run, self.__class__.reactions,
                                                             self.__class__.compounds, set(), biomass_equation)
        self.assertEqual(len(cp), 102)
        self.assertEqual(len(rc), 23)
        self.assertEqual(len(reactions), 34698)

    def test_run_fba(self):
        """Test running the fba. We build a run a complete FBA based on reaction_list.txt"""
        if media_file_loc == '':
            return
        self.assertTrue(os.path.exists(os.path.join(test_file_loc, 'reaction_list.txt')))
        self.assertTrue(os.path.exists(os.path.join(media_file_loc, 'ArgonneLB.txt')))
        compounds, reactions, enzymes = self.__class__.compounds, self.__class__.reactions, self.__class__.enzymes
        reactions2run = set()
        with open(os.path.join(test_file_loc, 'reaction_list.txt'), 'r') as f:
            for l in f:
                if l.startswith('#'):
                    continue
                if "biomass" in l.lower():
                    continue
                r = l.strip()
                if r in reactions:
                    reactions2run.add(r)
        media = PyFBA.parse.read_media_file(os.path.join(media_file_loc, 'ArgonneLB.txt'))
        biomass = PyFBA.metabolism.biomass_equation('gram_negative')
        status, value, growth = PyFBA.fba.run_fba(compounds, reactions, reactions2run, media, biomass, verbose=False)
        self.assertTrue(growth)
        value = float('%0.3f' % value)
        self.assertEqual(value, 282.674)
