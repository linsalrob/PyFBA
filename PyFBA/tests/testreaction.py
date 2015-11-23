import unittest
import PyFBA

"""
A class to test the reaction class. 
"""


class TestReaction(unittest.TestCase):

    def setUp(self):
        """This method is called before every test_ method"""
        self.reaction = PyFBA.metabolism.Reaction("test reaction")

    def test_reaction_name(self):
        """The name should be a string and not the null string"""
        name = self.reaction.name
        self.assertIsInstance(name, str)
        self.assertEqual(name, "test reaction")

    def test_reaction_left_compounds(self):
        """Left compounds are a set of compounds that are on the left of the reaction."""
        # check that we are provided with a set
        self.assertRaises(
            TypeError,
            self.reaction.add_left_compounds,
            "my node"
        )
        self.reaction.add_left_compounds({"a"})
        self.reaction.add_left_compounds({"b", "c"})
        self.assertEqual(self.reaction.number_of_left_compounds(), 3)
        # check to see that we have a set and not an array or other data
        # structure
        self.reaction.add_left_compounds({"a"})
        self.reaction.add_left_compounds({"b", "c"})
        self.assertEqual(self.reaction.number_of_left_compounds(), 3)

    def test_reaction_left_compounds_abundance(self):
        """Test adding the abundance of compounds on the left side of the equation"""
        self.reaction.add_left_compounds({"a"})
        # check that we raise an Exception if the compound is not there
        self.assertRaises(
            KeyError,
            self.reaction.set_left_compound_abundance,
            "c",
            1
        )
        # check that we raise an Exception if the value is a string
        self.assertRaises(
            TypeError,
            self.reaction.set_left_compound_abundance,
            "a",
            "1"
        )
        # check that it works if we add an int
        self.reaction.set_left_compound_abundance('a', 1)
        # check that it works if we add a float
        self.reaction.set_left_compound_abundance('a', 1.5)

    def test_reaction_left_compounds_abundance_retrieval(self):
        """Test getting the abundance of compounds on the left side of the equation"""
        self.reaction.add_left_compounds({"a"})
        self.reaction.set_left_compound_abundance('a', 1.5)
        # check we throw an exception 
        self.assertRaises(
            KeyError,
            self.reaction.get_left_compound_abundance,
            "c",
        )
        # check the value
        self.assertEqual(
            self.reaction.get_left_compound_abundance('a'), 1.5)
        self.assertEqual(
            self.reaction.left_abundance['a'], 1.5)

    def test_reaction_right_compounds(self):
        """Right compounds are a set of compounds that are on the right of the reaction."""
        # check that we are provided with a set
        self.assertRaises(
            TypeError,
            self.reaction.add_right_compounds,
            "my node"
        )
        self.reaction.add_right_compounds({"a"})
        self.reaction.add_right_compounds({"b", "c"})
        self.assertEqual(self.reaction.number_of_right_compounds(), 3)
        # check to see that we have a set and not an array or other data
        # structure
        self.reaction.add_right_compounds({"a"})
        self.reaction.add_right_compounds({"b", "c"})
        self.assertEqual(self.reaction.number_of_right_compounds(), 3)
    
    def test_reaction_right_compounds_abundance(self):
        """Test adding the abundance of compounds on the right side of the equation"""
        self.reaction.add_right_compounds({"a"})
        # check that we raise an Exception if the compound is not there
        self.assertRaises(
            KeyError,
            self.reaction.set_right_compound_abundance,
            "c",
            1
        )
        # check that we raise an Exception if the value is a string
        self.assertRaises(
            TypeError,
            self.reaction.set_right_compound_abundance,
            "a",
            "1"
        )
        # check that it works if we add an int
        self.reaction.set_right_compound_abundance('a', 1)
        # check that it works if we add a float
        self.reaction.set_right_compound_abundance('a', 1.5)

    def test_reaction_right_compounds_abundance_retrieval(self):
        """Test getting the abundance of compounds on the right side of the equation"""
        self.reaction.add_right_compounds({"a"})
        self.reaction.set_right_compound_abundance('a', 1.5)
        # check we throw an exception 
        self.assertRaises(
            KeyError,
            self.reaction.get_right_compound_abundance,
            "c",
        )
        # check the value
        self.assertEqual(
            self.reaction.get_right_compound_abundance('a'), 1.5)
        self.assertEqual(
            self.reaction.right_abundance['a'], 1.5)

    def test_all_compounds(self):
        """Test that we can get all the compounds back"""
        self.reaction.add_left_compounds({"x", "y", "z"})
        self.reaction.add_right_compounds({"a", "b", "c"})
        ans = {"x", "y", "z", "a", "b", "c"}
        self.assertEqual(self.reaction.all_compounds(), ans)

    def test_num_compounds(self):
        """Test that we can get the number of compounds back"""
        self.reaction.add_left_compounds({"x", "y", "z"})
        self.reaction.add_right_compounds({"a", "b", "c"})
        self.assertEqual(self.reaction.number_of_compounds(), 6)

    def test_equals(self):
        """Test the equals method defined for two reactions"""
        other_reaction = PyFBA.metabolism.Reaction("similar reaction")
        self.reaction.add_right_compounds({"a", "b", "c"})
        self.reaction.add_left_compounds({"1", "2", "3"})
        other_reaction.add_right_compounds({"a", "b", "c"})
        other_reaction.add_left_compounds({"1", "2", "3"})
        self.assertEqual(self.reaction, other_reaction)
        other_reaction.add_left_compounds({"4"})
        self.assertNotEqual(self.reaction, other_reaction)

    def test_equals_reverse(self):
        """Test the equals method defined for two reactions but when l=r and r=l"""
        other_reaction = PyFBA.metabolism.Reaction("similar reaction")
        self.reaction.add_right_compounds({"a", "b", "c"})
        self.reaction.add_left_compounds({"1", "2", "3"})
        other_reaction.add_left_compounds({"a", "b", "c"})
        other_reaction.add_right_compounds({"1", "2", "3"})
        self.assertEqual(self.reaction, other_reaction)

    def test_has(self):
        """Test for the presence of a compound in a reaction"""
        self.reaction.add_left_compounds({"1", "2", "3"})
        self.reaction.add_right_compounds({"a", "b", "c"})
        self.assertTrue(self.reaction.has("2"))
        self.assertTrue(self.reaction.has("b"))
        self.assertFalse(self.reaction.has("2b"))

    def test_opposite_sides(self):
        """Test whether two things are on opposite sides"""
        self.reaction.add_left_compounds({"1", "2", "3"})
        self.reaction.add_right_compounds({"a", "b", "c"})
        self.assertRaises(
            ValueError,
            self.reaction.opposite_sides,
            "not present",
            "a"
        )
        self.assertRaises(
            ValueError,
            self.reaction.opposite_sides,
            "1",
            "not present"
        )
        self.assertFalse(self.reaction.opposite_sides("1", "2"))
        self.assertFalse(self.reaction.opposite_sides("a", "c"))
        self.assertTrue(self.reaction.opposite_sides("1", "a"))
        self.assertTrue(self.reaction.opposite_sides("1", "b"))
        self.assertTrue(self.reaction.opposite_sides("1", "c"))

    def test_p_LR(self):
        """Test the probability of running left to right"""
        self.assertEqual(self.reaction.pLR, 0)
        self.reaction.set_probability_left_to_right(10)
        self.assertEqual(
                self.reaction.get_probability_left_to_right(), 10)
        self.reaction.pLR = 5
        self.assertEqual(self.reaction.pLR, 5)
        self.assertRaises(
            TypeError,
            self.reaction.set_probability_left_to_right,
            "not an int"
        )

    def test_p_RL(self):
        """Test the probability of running right to left"""
        self.assertEqual(self.reaction.pRL, 0)
        self.reaction.set_probability_right_to_left(20)
        self.assertEqual(
                self.reaction.get_probability_right_to_left(), 20)
        self.reaction.pRL = 15
        self.assertEqual(self.reaction.pRL, 15)
        self.assertRaises(
            TypeError,
            self.reaction.set_probability_right_to_left,
            "not an int"
        )

    def test_deltaG(self):
        """Test the delta G of the reaction"""
        self.assertEqual(self.reaction.deltaG, 0)
        self.reaction.set_deltaG(20)
        self.assertEqual(
                self.reaction.get_deltaG(), 20)
        self.reaction.deltaG = 15
        self.assertEqual(self.reaction.deltaG, 15)
        self.assertRaises(
            TypeError,
            self.reaction.set_deltaG,
            "not an int"
        )

    def test_enzymes(self):
        """Test adding enzymes to reactions and complexes to enzymes"""
        self.reaction.add_enzymes({'a', 'b', 'c'})
        self.assertTrue(self.reaction.has_enzyme("a"))
        self.assertFalse(self.reaction.has_enzyme("A"))
        self.assertEqual(self.reaction.number_of_enzymes(), 3)
        # check that it is a set
        self.reaction.add_enzymes({'a', 'b', 'c'})
        self.assertEqual(self.reaction.number_of_enzymes(), 3)
        # check that is only accepts sets
        self.assertRaises(
            TypeError,
            self.reaction.add_enzymes,
            "Not a set"
        )
        self.assertEqual({'a', 'b', 'c'}, self.reaction.all_enzymes())

    def test_pegs(self):
        """Test identifying pegs in reactions"""
        self.reaction.add_pegs({"1", "w", "2", "a"})
        self.assertTrue(self.reaction.has_peg("w"))
        self.assertFalse(self.reaction.has_peg("W"))
        self.assertRaises(
            TypeError,
            self.reaction.add_pegs,
            "Not a set"
        )

    def test_input_output(self):
        """Input and output should start out false, but we should be able to make them true"""
        self.assertFalse(self.reaction.is_input_reaction())
        self.assertFalse(self.reaction.is_output_reaction())
        self.reaction.toggle_input_reaction()
        self.reaction.toggle_output_reaction()
        self.assertTrue(self.reaction.is_input_reaction())
        self.assertTrue(self.reaction.is_output_reaction())

    def test_reverse_reaction(self):
        """Convert everything from the left side to the right side and vice versa"""
        self.reaction.add_left_compounds({'a', 'b', 'c'})
        self.reaction.add_right_compounds({'x', 'y', 'z'})
        self.reaction.set_left_compound_abundance('a', 1)
        self.reaction.set_left_compound_abundance('b', 2)
        self.reaction.set_left_compound_abundance('c', 3)
        self.reaction.set_right_compound_abundance('x', 7)
        self.reaction.set_right_compound_abundance('y', 8)
        self.reaction.set_right_compound_abundance('z', 9)
        self.reaction.set_probability_right_to_left(10)
        self.reaction.set_probability_left_to_right(20)

        self.reaction.reverse_reaction()

        self.assertEqual(self.reaction.left_compounds, {'x', 'y', 'z'})
        self.assertEqual(self.reaction.right_compounds, {'a', 'b', 'c'})
        self.assertEqual(self.reaction.get_left_compound_abundance('y'), 8)
        self.assertEqual(self.reaction.get_right_compound_abundance('b'), 2)
        self.assertEqual(self.reaction.get_probability_right_to_left(), 20)
        self.assertEqual(self.reaction.get_probability_left_to_right(), 10)


if __name__ == '__main__':
    unittest.main()

