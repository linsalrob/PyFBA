import unittest
import PyFBA
from PyFBA import lp

class TestBounds(unittest.TestCase):

    def test_row_bounds(self):
        '''Testing the assertion of row bounds and the logic'''

        # define a media compound
        media_compound = PyFBA.metabolism.CompoundWithLocation("m0011", 'Media', 'e')
        media = {media_compound}

        # define an FBA matrix of the right size
        newfba = lp.load([[1,2,3,4,5,6]])

        # make up some reactions
        ra = PyFBA.metabolism.Reaction('rA', 'reaction A')
        ra.direction = '>'
        rb = PyFBA.metabolism.Reaction('rB', 'reaction B')
        rb.direction = '<'
        rc = PyFBA.metabolism.Reaction('rC', 'reaction C')
        rc.direction = '='
        rd = PyFBA.metabolism.Reaction('rUS', 'reaction US')
        rd.direction = '='
        rd.is_uptake_secretion = True
        re = PyFBA.metabolism.Reaction('biomass_equation', 'biomass_equation')
        re.direction = '>'
        rf = PyFBA.metabolism.Reaction('rUSM', 'reaction US Media')
        rf.direction = '='
        rf.is_uptake_secretion = True
        rf.add_left_compounds({media_compound})

        # test our reactions
        all_reactions = {'ra': ra, 'rb': rb, 'rc': rc, 'rd': rd, 're': re, 'rf': rf}
        reactions2run = set(all_reactions.keys())
        rbvals = PyFBA.fba.reaction_bounds(all_reactions, reactions2run, media)

        self.assertEqual(rbvals['ra'], (0, 1000))
        self.assertEqual(rbvals['rb'], (-1000, 1000))
        self.assertEqual(rbvals['rc'], (-1000, 1000))
        self.assertEqual(rbvals['rd'], (-1000, 1000))
        self.assertEqual(rbvals['re'], (0, 1000))
        self.assertEqual(rbvals['rf'], (-1000, 1000))

    def test_col_bounds(self):
        '''Testing the assertion of column bounds'''
        # define an FBA matrix of the right size
        newfba = lp.load([[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]])
        ca = PyFBA.metabolism.CompoundWithLocation("c1", "A", "e")
        cb = PyFBA.metabolism.CompoundWithLocation("c2", "B", "c")
        cc = PyFBA.metabolism.CompoundWithLocation("c3", "C", "h")
        compounds = {ca, cb, cc}
        col_bounds = PyFBA.fba.compound_bounds(compounds)
        self.assertEqual(col_bounds[ca], (0, 0))
        self.assertEqual(col_bounds[cb], (0, 0))
        self.assertEqual(col_bounds[cc], (0, 0))


if __name__ == '__main__':
    unittest.main()
