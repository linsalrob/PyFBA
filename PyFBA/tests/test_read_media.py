import os
import unittest

import sys

import PyFBA

class TestReadMedia(unittest.TestCase):

    def test_read_media_file(self):
        """Test reading a media file"""
        media = PyFBA.parse.raw_media('ArgonneLB')
        self.assertEqual(len(media), 65)
        testm = None
        for m in  media:
            if m.name == 'D-Glucose':
                testm = m
                break
        self.assertIn(testm, media)


if __name__ == '__main__':
    unittest.main()
