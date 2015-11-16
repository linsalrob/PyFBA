from unittest import TestCase
from metabolism import biomass_equation

__author__ = 'Rob Edwards'


class TestBiomass_equation(TestCase):
    """
    The methods in the class really just test that the biomass equations are Reaction objects and that they
    are the right length. There is not a lot to test for biomass equations!
    """

    def test_kbase_length(self):
        """
        Test the length of the kbase biomass reactions.
        """
        bme = biomass_equation('kbase')
        self.assertEqual(bme.number_of_left_compounds(), 76)
        self.assertEqual(bme.number_of_right_compounds(), 9)

    def test_kbase_simple_length(self):
        """
        Test the length of the kbase biomass reactions.
        """
        bme = biomass_equation('kbase_simple')
        self.assertEqual(bme.number_of_left_compounds(), 73)
        self.assertEqual(bme.number_of_right_compounds(), 9)

    def test_gramnegative_length(self):
        """
        Test the length of the kbase biomass reactions.
        """
        bme = biomass_equation('gram_negative')
        self.assertEqual(bme.number_of_left_compounds(), 53)
        self.assertEqual(bme.number_of_right_compounds(), 7)
