import unittest
import glpk


class TestLinearProgramming(unittest.TestCase):

    def setUp(self):
        """This method is called before every test_ method"""
        self.lp = glpk.LPX()
        # initiate the LP and set the objective to maximize
        self.lp.erase()
        self.lp.obj.maximize = True

    def test_solve(self):
        """Test the complete linear programing solution, using the example from the documentation"""
        # define the number of rows and columns
        self.lp.rows.add(3)
        self.lp.cols.add(3)

        # set up the matrix
        self.lp.matrix = [1.0, 1.0, 1.0, 10.0, 4.0, 5.0, 2.0, 2.0, 6.0]
        # self.lp.matrix = [1.0, 10.0, 2.0, 1.0, 4.0, 2.0, 1.0, 5.0, 6.0]

        # the objective coefficient
        self.lp.obj[:] = [10.0, 6.0, 4.0]

        # add the row and column bounds
        rbounds = [(None, 100.0), (None, 600.0), (None, 300.0)]
        cbounds = [(0, None), (0, None), (0, None)]
        for i in range(3):
            self.lp.cols[i].bounds = cbounds[i]
            self.lp.rows[i].bounds = rbounds[i]

        # solve the LP
        self.lp.simplex()

        # convert the output to 3 sd to test
        r = "%0.3f" % self.lp.obj.value
        self.assertEqual(r, "733.333")
        self.assertEqual(self.lp.status, 'opt')
