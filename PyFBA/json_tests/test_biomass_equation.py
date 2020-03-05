from unittest import TestCase
import PyFBA


class TestBiomass_equation(TestCase):
    """
    The methods in the class really just test that the biomass_equation equations are Reaction objects and that they
    are the right length. There is not a lot to test for biomass_equation equations!
    """

    def test_kbase_length(self):
        """
        Test the length of the kbase biomass_equation reactions.
        """
        biomass_eqn = PyFBA.metabolism.biomass_equation('kbase')
        self.assertEqual(biomass_eqn.number_of_left_compounds(), 76)
        self.assertEqual(biomass_eqn.number_of_right_compounds(), 9)

    def test_kbase_simple_length(self):
        """
        Test the length of the kbase biomass_equation reactions.
        """
        biomass_eqn = PyFBA.metabolism.biomass_equation('kbase_simple')
        self.assertEqual(biomass_eqn.number_of_left_compounds(), 73)
        self.assertEqual(biomass_eqn.number_of_right_compounds(), 9)

    def test_gramnegative_length(self):
        """
        Test the length of the kbase biomass_equation reactions.
        """
        biomass_eqn = PyFBA.metabolism.biomass_equation('gram_negative')
        self.assertEqual(biomass_eqn.number_of_left_compounds(), 53)
        self.assertEqual(biomass_eqn.number_of_right_compounds(), 7)
