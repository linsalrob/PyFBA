import unittest
import PyFBA

"""
A class to test the compound class. 
"""


class TestEnzyme(unittest.TestCase):

    def setUp(self):
        """This method is called before every test_ method"""
        self.enz = PyFBA.metabolism.Enzyme('test enz')
        self.enz.roles = {'a', 'b', 'c'}

    def test_equals(self):
        """Does the equals function work"""
        otherenz = PyFBA.metabolism.Enzyme('test enz')
        otherenz.roles = {'a', 'b', 'c'}
        self.assertEqual(self.enz, otherenz)

    def test_ne(self):
        """Does the not equals work for both roles and name"""
        otherenz = PyFBA.metabolism.Enzyme('test enz')
        otherenz.roles = {'a', 'b'}
        self.assertNotEquals(self.enz, otherenz)
        otherenz = PyFBA.metabolism.Enzyme('test enz2')
        otherenz.roles = {'a', 'b', 'c'}
        self.assertNotEquals(self.enz, otherenz)

    def test_hash(self):
        """How does the hash function work"""
        h = hash('test enz')
        self.assertEqual(h, hash(self.enz))

    def test_add_roles(self):
        """Test adding the roles and that we are using sets"""
        # we should get an error unless we add a set
        self.assertRaises(
            TypeError,
            self.enz.add_roles,
            '123'
        )
        self.enz.add_roles({'m', 'n', 'o'})
        self.enz.add_roles({'m', 'n', 'o'})
        self.assertEqual(len(self.enz.roles), 6)

    def has_role(self):
        """Test the has role function"""
        self.assertTrue(self.enz.has_role('b'))
        self.assertFalse(self.enz.has_role('B'))

    def test_number_of_roles(self):
        """ Test how many roles we have"""
        self.assertEqual(self.enz.number_of_roles(), 3)

    def test_add_pegs(self):
        """Check that we can add a dictionary of pegs"""
        # check that we have to submit a dict
        self.assertRaises(
            TypeError,
            self.enz.add_pegs,
            'peg1'
        )
        # check that we have to submit a dict
        self.assertRaises(
            TypeError,
            self.enz.add_pegs,
            {'peg1', 'peg2'}
        )
        # this should work
        self.enz.add_pegs({'p1': 'a', 'p2': 'b'})
        # check that the role must be present
        self.assertRaises(
            KeyError,
            self.enz.add_pegs,
            {'p3': 'role not present!'}
        )

    def test_add_a_peg(self):
        """Test adding a single peg"""
        # we have to supply a string
        self.assertRaises(
            TypeError,
            self.enz.add_a_peg,
            {'p1', 'p2'},
            {'r1', 'r2'}
        )
        # this should work
        self.enz.add_a_peg('p1', 'a')

        # check the role is there
        self.assertRaises(
            KeyError,
            self.enz.add_a_peg,
            'p3',
            'role not present!'
        )

    def test_number_of_pegs(self):
        """ How many pegs do we have? """
        self.assertEqual(self.enz.number_of_pegs(), 0)
        self.enz.add_pegs({'p1': 'a', 'p2': 'b', 'p3': 'c'})
        self.assertEqual(self.enz.number_of_pegs(), 3)

    def test_add_reaction(self):
        """ Test adding a reaction"""
        self.assertRaises(
            TypeError,
            self.enz.add_reaction,
            {'r1', 'r2'}
        )
        self.enz.add_reaction('a reaction')
        self.assertEqual(len(self.enz.reactions), 1)
        # check this is a set
        self.enz.add_reaction('a reaction')
        self.assertEqual(len(self.enz.reactions), 1)

    def test_number_of_reactions(self):
        """Test how many reactions we have"""
        self.enz.add_reaction('a reaction')
        self.assertEqual(self.enz.number_of_reactions(), 1)
        self.enz.add_reaction("Another reaction")
        self.assertEqual(self.enz.number_of_reactions(), 2)
    
    def test_probability(self):
        """Test that we get a probability back"""
        self.enz.add_a_peg('123', 'a')
        self.assertEqual(self.enz.probability(), (1.0 * 1/3))

    def test_has_peg_for_role(self):
        """Test whether we have a peg for a role"""
        self.enz.add_a_peg('123', 'a')
        self.assertTrue(self.enz.has_peg_for_role('a'))
        self.assertFalse(self.enz.has_peg_for_role('b'))
        self.assertFalse(self.enz.has_peg_for_role('c'))


