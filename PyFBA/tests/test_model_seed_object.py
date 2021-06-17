import os
import sys
import unittest
import PyFBA

class TestModelSeedObject(unittest.TestCase):
    def test_object(self):
        msp = PyFBA.model_seed.ModelData()
        self.assertIsNone(msp.compounds)
        self.assertIsNone(msp.enzymes)
        self.assertFalse(msp.reactions)

    def test_add_to_object(self):
        # we just set these up as some dicts
        msp = PyFBA.model_seed.ModelData()
        self.assertIsNone(msp.compounds)
        self.assertIsNone(msp.enzymes)
        self.assertFalse(msp.reactions)
        msp.compounds = {'a' : 1, 'b' : 2}
        self.assertDictEqual({'a': 1, 'b': 2}, msp.compounds)
        msp.enzymes = {'c' : 3, 'd' : 4}
        self.assertDictEqual({'c': 3, 'd': 4}, msp.enzymes)
        msp.reactions['test'] = {'e' : 5, 'f' : 6}
        self.assertDictEqual({'test': {'e' : 5, 'f' : 6}}, msp.reactions)
        msp.reset()
        self.assertFalse(msp.compounds)
        self.assertIsNone(msp.enzymes)
        self.assertFalse(msp.reactions)

if __name__ == '__main__':
    unittest.main()
