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
    modeldata = PyFBA.parse.model_seed.parse_model_seed_data('gramnegative', verbose=True)
    media = PyFBA.parse.pyfba_media("ArgonneLB")

    def setUp(self):
        """
        This is run before everything else
        """

    def test_external_reactions(self):
        """Testing the fba external reactions"""
        compounds = self.__class__.modeldata.compounds
        # create a few external compounds
        cpds = list(compounds)[0:10]
        model_cpds = set()
        for c in cpds:
            new_comp = PyFBA.metabolism.CompoundWithLocation.from_compound(c, 'e')
            model_cpds.add(new_comp)

        upsec = PyFBA.fba.uptake_and_secretion_reactions(model_cpds, self.__class__.media)
        self.assertEqual(len(upsec), 10)

        # now remove them and see there is nothing left
        emptyset = PyFBA.fba.remove_uptake_and_secretion_reactions(upsec)
        self.assertEqual(len(emptyset), 0)

    def test_create_sm(self):
        """Test the stoichiometric matrix"""

        reactions2run = set(list(self.__class__.modeldata.reactions.keys())[0:20])
        biomass_equation = PyFBA.metabolism.biomass_equation('gram_negative')
        cp, rc, upsr = PyFBA.fba.create_stoichiometric_matrix(reactions_to_run=reactions2run,
                                                                   modeldata=self.__class__.modeldata,
                                                                   media=set(),
                                                                   biomass_equation=biomass_equation)
        # this allows some wiggle room as the data changes
        self.assertGreaterEqual(len(cp), 100)
        self.assertLessEqual(len(cp), 150)
        self.assertGreaterEqual(len(rc), 20)
        self.assertLessEqual(len(rc), 30)
        self.assertGreaterEqual(len(upsr), 34000)
        self.assertLessEqual(len(upsr), 35000)

    def test_run_fba(self):
        """Test running the fba. We build a run a complete FBA based on reaction_list.txt"""
        if media_file_loc == '':
            return
        self.assertTrue(os.path.exists(os.path.join(test_file_loc, 'reaction_list.txt')))
        self.assertTrue(os.path.exists(os.path.join(media_file_loc, 'ArgonneLB.txt')))
        reactions2run = set()
        with open(os.path.join(test_file_loc, 'reaction_list.txt'), 'r') as f:
            for lf in f:
                if lf.startswith('#'):
                    continue
                if "biomass" in lf.lower():
                    continue
                r = lf.strip()
                if r in self.__class__.modeldata.reactions:
                    reactions2run.add(r)
        media = PyFBA.parse.read_media_file(os.path.join(media_file_loc, 'ArgonneLB.txt'))
        biomass = PyFBA.metabolism.biomass_equation('gram_negative')
        status, value, growth = PyFBA.fba.run_fba(self.__class__.modeldata, reactions2run, media, biomass,
                                                  verbose=False)
        self.assertTrue(growth)
        value = float('%0.3f' % value)
        self.assertEqual(value, 340.873)
