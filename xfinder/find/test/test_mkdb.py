import unittest
import os
import shutil
import xfinder.find.mkdb
import filecmp


class TestMkdb(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Runs once before all tests
        cls.tempdir = "test_temp"
        os.mkdir(cls.tempdir)

    @classmethod
    def tearDownClass(cls):
        # Runs once after all tests
        shutil.rmtree(cls.tempdir)
        
    def test_create_db_tables(self):
        refdb = os.path.join(os.path.dirname(__file__), "data", "database.db")
        testdb = os.path.join(self.tempdir, "testdb.db")
        xfinder.find.mkdb.create_db_tables(testdb)
        self.assertTrue(filecmp.cmp(testdb, refdb, shallow=False))
