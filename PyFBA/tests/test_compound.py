import unittest
import PyFBA

"""
A class to test the compound class. 
"""


class TestCompound(unittest.TestCase):

    def setUp(self):
        """This method is called before every test_ method"""
        self.compound = PyFBA.metabolism.Compound("t1", "test compound")
        self.compound.abbreviation = "Cool"
        self.compound.add_attribute('What', "Everything")
        self.compound_with_loc = PyFBA.metabolism.CompoundWithLocation.from_compound(self.compound, "extracellular")

    def test_equals(self):
        """Test that our equals function works"""
        othercompound = PyFBA.metabolism.Compound("t2", "test compound")
        self.assertEqual(self.compound, othercompound)
        othercompound.name = "Another compound"
        self.assertNotEqual(self.compound, othercompound)

    def test_hash(self):
        """Test that our hash function works"""
        self.assertEqual(hash(self.compound), hash(("t1", "test compound")))

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

    def test_abbreviation(self):
        """ Did we get the new abbreviation?"""
        self.assertEqual(self.compound.abbreviation, "Cool")

    def test_adding_attributes(self):
        """ Did we add the new attributes"""
        self.assertEqual(self.compound.get_attribute("What"), "Everything")

    def test_compound_with_location(self):
        """Test the location of a compound"""
        self.assertEqual(self.compound_with_loc.location, 'extracellular')

    def test_comp_with_loc_copied(self):
        """Test we copied all attributes properly"""
        self.assertEqual(self.compound_with_loc.get_attribute("What"), "Everything")
