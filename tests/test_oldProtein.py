import pytest
import os
from biofkit.proteinKit import oldProtein

class Test_pdb2List():
    """A class used to test function 'pdb2List'."""

    @classmethod
    def test_read_a_pdb(self):
        """Test whether function 'oldProtein.pdb2List' can read a pdb."""
        proList: list = oldProtein.pdb2List(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        assert (proList), "oldProtein.pdb2List did not read a PDB properly!"


class Test_pdb2Dict():
    """A class used to test function 'pdb2Dict'."""
    @classmethod
    def test_read_a_pdb(self):
        """Test whether function 'oldProtein.pdb2Dict' can read a pdb."""
        proDict: list = oldProtein.pdb2Dict(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        assert (proDict), "oldProtein.pdb2Dict did not read a PDB properly!"
    

class Test_proDict2ProList():
    """A class used to test function 'proDict2ProList'"""
    @classmethod
    def test_proDict2ProListRight(self):
        """Test whether function 'proDict2ProList' can transfer a protein structure dictionary to a list properly."""
        proDict: dict = oldProtein.pdb2Dict(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        proList: list = oldProtein.pdb2List(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        transferedProList: list = oldProtein.proDict2ProList(rawDict=proDict)
        assert (transferedProList == proList), "Convert failed!"


class Test_proListIsValid():
    """A class used to test function 'proListIsValid'"""
    @classmethod
    def test_distinguishAnInvalidList(self):
        """Test whether function 'proListIsValid' can distinguish an invalid proList"""
        proList: list = [["yuanshen"], ["qidong"]]
        with pytest.raises(IndexError):
            assert (not oldProtein.proListIsValid(proList))

    @classmethod
    def test_distinguishAValidList(self):
        """Test whether function 'proListIsValid' can distinguish a valid proList"""
        proList: list = oldProtein.pdb2List(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        assert (oldProtein.proListIsValid(proList)), "oldProtein.proListIsValid did not trust a valid proList!"
        
    @classmethod
    def test_distinguishAConvertedProList(self):
        """Test whether function 'proListIsValid' trust a converted proList"""
        proDict: dict = oldProtein.pdb2Dict(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        proList: list = oldProtein.proDict2ProList(proDict)
        assert (oldProtein.proListIsValid(proList)), "oldProtein.proListIsValid did not trust a converted proList!"


class Test_proDictIsValid():
    """A class used to test function 'proDictIsValid'"""
    @classmethod
    def test_distinguishAnInvalidDict(self):
        """Test whether function 'proDictIsValid' can distinguish an invalid proDict"""
        proDict: dict = {"yuanshen": ["qidong"]}
        with pytest.raises(KeyError):
            assert (not oldProtein.proDictIsValid(proDict))

    @classmethod
    def test_distinguishAValidDuct(self):
        """Test whether function 'proDictIsValid' can distinguish a valid proDict"""
        proDict: dict = oldProtein.pdb2Dict(os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        assert (oldProtein.proDictIsValid(proDict)), "oldProtein.proDictIsValid did not trust a valid proDict!"

    @classmethod
    def test_distinguishAConvertedProDict(self):
        """Test whether function 'proDictIsValid' trust a converted proDict"""
        proDict: dict = oldProtein.pdb2Dict(pdbFilePath=os.path.join(os.path.dirname(__file__), "pdb", "1A.pdb"))
        assert (oldProtein.proDictIsValid(proDict)), "oldProtein.proDictIsValid did not trust a converted proDict!"

if __name__ == "__main__":
    pytest.main()
