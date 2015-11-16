import os
import unittest
from parse import read_media_file
from metabolism import Compound

class TestReadMedia(unittest.TestCase):
    def setUp(self):
        self.assertTrue(os.path.exists('media/ArgonneLB.txt'))

    def test_read_media_file(self):
        media = read_media_file('media/ArgonneLB.txt')
        self.assertEqual(len(media), 65)
        gluc = Compound('D-Glucose', 'e')
        self.assertIn(gluc, media)


if __name__ == '__main__':
    unittest.main()
