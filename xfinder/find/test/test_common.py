import unittest
import os
import xfinder.find.common


class TestCommon(unittest.TestCase):


    def test_get_git_root(self):
        self.git_root = xfinder.find.common.get_git_root()
        self.assertEqual(os.path.basename(self.git_root), "Xfinder")

    def test_get_database_dir(self):
        self.database_dir= xfinder.find.common.get_database_dir()
        self.assertTrue(os.path.exists(self.database_dir))

    def test_timestamp(self):
        # Dont know how to test this
        pass
        