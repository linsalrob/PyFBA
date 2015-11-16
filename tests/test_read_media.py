import os
import unittest


class TestReadMedia(unittest.TestCase):
    def setUp(self):
        self.assertTrue(os.path.exists('media/ArgonneLB.txtxxx'))

    def test_read_media_file(self):
        media = read_media_file('media/ArgonneLB.txt')
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
