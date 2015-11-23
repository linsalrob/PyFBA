import os
import unittest

import PyFBA
from parse import read_media_file

media_file_loc = ''
if os.path.exists('../media/ArgonneLB.txt'):
    media_file_loc = '../media'
elif os.path.exists('media/ArgonneLB.txt'):
    media_file_loc = 'media'


class TestReadMedia(unittest.TestCase):
    def setUp(self):
        self.assertTrue(os.path.exists(os.path.join(media_file_loc, 'ArgonneLB.txt')))

    def test_read_media_file(self):
        """Test reading a media file"""
        media = read_media_file(os.path.join(media_file_loc, 'ArgonneLB.txt'))
        self.assertEqual(len(media), 65)
        gluc = PyFBA.metabolism.Compound('D-Glucose', 'e')
        self.assertIn(gluc, media)


if __name__ == '__main__':
    unittest.main()
