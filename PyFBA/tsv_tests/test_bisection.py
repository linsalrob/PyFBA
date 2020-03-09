import unittest

import PyFBA


class BisectionTest(unittest.TestCase):
    def test_bisect(self):
        alist = [1,2,3,4,5,6]
        lista, listb = PyFBA.gapfill.bisections.bisect(alist)
        self.assertEqual(lista, [1, 3, 5])
        self.assertEqual(listb, [2, 4, 6])

    def test_percent_split(self):
        alist = [1,2,3,4,5,6]
        lista, listb = PyFBA.gapfill.bisections.percent_split(alist, 50)
        self.assertEqual(lista, [1, 2, 3])
        self.assertEqual(listb, [4, 5, 6])
        lista, listb = PyFBA.gapfill.bisections.percent_split(alist, 100.0 / 3)
        self.assertEqual(lista, [1, 2])
        self.assertEqual(listb, [3, 4, 5, 6])

    def test_split_by_clusters(self):
        alist = [1,2,3,4,5,6]
        clusters = {1:1, 3:1, 4:1, 2:2, 5:2, 6:2}
        lista, listb = PyFBA.gapfill.bisections.optimize_split_by_rclust(alist, clusters, 50)
        self.assertEqual(lista, [1, 3, 4])
        self.assertEqual(listb, [2, 5, 6])



if __name__ == '__main__':
    unittest.main()
