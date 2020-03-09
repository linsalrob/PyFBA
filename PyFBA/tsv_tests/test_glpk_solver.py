import unittest
from PyFBA.tests.assertDeepAlmostEqual import assertDeepAlmostEqual
from PyFBA import lp

"""
A class to test the running of the FBA.

"""

class TestGLPKSolver(unittest.TestCase):

    def setUp(self):
        """This method is called before every test_ method"""

    def test_matrix(self):
        """Just initialize with a simple matrix"""
        mat = [
                [1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]
        ]

        lp.load(mat)

    def test_headers(self):
        """Test adding headers to the rows and columns"""
        mat = [
                [1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]
        ]

        rh = ['a', 'b', 'c']
        ch = ['w', 'x', 'y', 'z']

        lp.load(mat, rh, ch)

        rh.append('should fail')
        self.assertRaises(ValueError, lp.load, mat, rh, ch)

        rh = ['a', 'b', 'c']
        lp.load(mat, rh, ch)

        ch.append("also should fail")
        self.assertRaises(ValueError, lp.load, mat, rh, ch)

    def test_bound_rows(self):
        """Test adding tuples of boundary conditions for rows"""
        mat = [
                [1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]
        ]
        lp.load(mat)

        boundsr = []
        # test that we don't have enough values
        self.assertRaises(ValueError, lp.row_bounds, boundsr)

        # a boundary of None should  be infiniity!
        boundsr = [(None, 1000), (1000, None), (1000, 1000)]
        lp.row_bounds(boundsr)

    def test_bound_cols(self):
        """Test adding tuples of boundary conditions for columns"""
        mat = [
                [1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]
        ]
        lp.load(mat)

        boundsc = []
        # test that we don't have enough values
        self.assertRaises(ValueError, lp.col_bounds, boundsc)

        # a boundary of None should  be infiniity!
        boundsc = [(None, 1000), (1000, None), (1000, 1000), (None, None)]
        lp.col_bounds(boundsc)


    def test_objective_coeff(self):
        """Test adding the objective coefficients. Note that there
        should be as many coefficients as columns in the matrix."""

        mat = [
                [ 1.0, 1.0, 1.0],
                [10.0, 4.0, 5.0],
                [ 2.0, 2.0, 6.0],
                [ 2.0, 2.0, 6.0],
        ]
        lp.load(mat)
        lp.objective_coefficients([ 10.0, 6.0, 4.0 ])

    def test_solve(self):
        """Test the complete linear programing solution, using the 
        example from the documentation"""
        mat = [
                [ 1.0, 1.0, 1.0],
                [10.0, 4.0, 5.0],
                [ 2.0, 2.0, 6.0],
        ]
        lp.load(mat)
        lp.objective_coefficients([ 10.0, 6.0, 4.0 ])
        lp.row_bounds([(None, 100.0), (None, 600.0), (None, 300.0)])
        lp.col_bounds([(0, None), (0, None), (0, None)])
        status, result = lp.solve()
        r = "%0.3f" % result
        self.assertEqual(r, "733.333")
        self.assertEqual(status, 'opt')


    def test_primal_hash(self):
        """Test getting the primals back as a hash"""
        mat = [
                [ 1.0, 1.0, 1.0],
                [10.0, 4.0, 5.0],
                [ 2.0, 2.0, 6.0],
        ]
        rh = ['a', 'b', 'c']
        ch = ['x', 'y', 'z']

        lp.load(mat, rh, ch)
        lp.objective_coefficients([ 10.0, 6.0, 4.0 ])
        lp.row_bounds([(None, 100.0), (None, 600.0), (None, 300.0)])
        lp.col_bounds([(0, None), (0, None), (0, None)])
        status, result = lp.solve()
        col_pri = {'y': 66.66666666666666, 'x': 33.333333333333336, 'z': 0.0}
        col_res = lp.col_primal_hash()
        assertDeepAlmostEqual(self, col_pri, col_res, places=5)
        #self.assertEqual(col_pri, col_res)

    def test_primals(self):
        """Test getting the primals back as a list"""
        mat = [
                [ 1.0, 1.0, 1.0],
                [10.0, 4.0, 5.0],
                [ 2.0, 2.0, 6.0],
        ]
        rh = ['a', 'b', 'c']
        ch = ['x', 'y', 'z']

        lp.load(mat, rh, ch)
        lp.objective_coefficients([ 10.0, 6.0, 4.0 ])
        lp.row_bounds([(None, 100.0), (None, 600.0), (None, 300.0)])
        lp.col_bounds([(0, None), (0, None), (0, None)])
        status, result = lp.solve()
        col_pri = [33.333333333333336, 66.66666666666666, 0.0]
        col_res = lp.col_primals()
        assertDeepAlmostEqual(self, col_pri, col_res, places=5)
        #self.assertEqual(col_pri, col_res)
        







        


