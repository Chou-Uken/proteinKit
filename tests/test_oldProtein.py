import unittest
import os
from biofkit.proteinKit import oldProtein

class Test_pdb2List(unittest.TestCase):
    """A class used to test function 'pdb2List'."""

    @classmethod
    def test_read_a_pdb(self):
        """Test whether function 'oldProtein.pdb2List' can read a pdb."""
        proList: list = oldProtein.pdb2List(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        self.assertIsNotNone(None, "1")


class Test_pdb2Dict(unittest.TestCase):
    """A class used to test function 'pdb2Dict'."""
    @classmethod
    def test_read_a_pdb(self):
        """Test whether function 'oldProtein.pdb2Dict' can read a pdb."""
        proDict: list = oldProtein.pdb2Dict(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        self.assertIsNotNone(None, "2")


class Test_proDict2ProList(unittest.TestCase):
    """A class used to test function 'proDict2ProList'"""
    @classmethod
    def test_proDict2ProListRight(self):
        """Test whether function 'proDict2ProList' can transfer a protein structure dictionary to a list properly."""
        proDict: dict = oldProtein.pdb2Dict(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        proList: list = oldProtein.pdb2List(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        transferedProList: list = oldProtein.proDict2ProList(rawDict=proDict)
        self.assertListEqual(transferedProList, proList)



if __name__ == "__main__":
    unittest.main()
