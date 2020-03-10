import unittest
import PyFBA

"""
A class to test the compound class. 
"""


class TestCompound(unittest.TestCase):

    def setUp(self):
        """This method is called before every test_ method"""
        self.compound = PyFBA.metabolism.Compound("test compound", 'c')

    def test_equals(self):
        """Test that our equals function works"""
        othercompound = PyFBA.metabolism.Compound("test compound", 'c')
        self.assertEqual(self.compound, othercompound)
        othercompound.name = "Another compound"
        self.assertNotEqual(self.compound, othercompound)

    def test_hash(self):
        """Test that our hash function works"""
        self.assertEqual(hash(self.compound), hash(("test compound", 'c')))

    def test_location(self):
        """Test the location of a compound"""
        self.assertEqual(self.compound.location, 'c')

    def test_in_reactions(self):
        """Test which reactions the compound is in"""
        self.compound.add_reactions({"a", "b", "c"})
        self.compound.add_reactions({"d"})
        self.assertTrue(self.compound.has_reaction("a"))
        self.assertEqual(self.compound.number_of_reactions(), 4)
        # check that it is a set
        self.compound.add_reactions({"a", "b", "c"})
        self.assertEqual(self.compound.number_of_reactions(), 4)
        # check that we are provided with a set
        self.assertRaises(
            TypeError,
            self.compound.add_reactions,
            "A reaction"
        )
