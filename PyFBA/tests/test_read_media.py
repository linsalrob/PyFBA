import os
import unittest

import sys

import PyFBA

media_file_loc = ''
if os.path.exists('../media/ArgonneLB.txt'):
    media_file_loc = '../media'
elif os.path.exists('media/ArgonneLB.txt'):
    media_file_loc = 'media'
elif 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], 'ArgonneLB.txt')):
    media_file_loc = os.environ['PYFBA_MEDIA_DIR']
else:
    sys.stderr.write("No media found. Can't proceed with the tests.\n")
    sys.stderr.write("You can specify the media location by setting the PYFBA_MEDIA_DIR environment variable\n")


class TestReadMedia(unittest.TestCase):
    def setUp(self):
        self.assertTrue(os.path.exists(os.path.join(media_file_loc, 'ArgonneLB.txt')))

    def test_read_media_file(self):
        """Test reading a media file"""
        if media_file_loc == "":
            return
        media = PyFBA.parse.read_media_file(os.path.join(media_file_loc, 'ArgonneLB.txt'))
        self.assertEqual(len(media), 65)
        gluc = PyFBA.metabolism.Compound('D-Glucose', 'e')
        self.assertIn(gluc, media)


if __name__ == '__main__':
    unittest.main()
