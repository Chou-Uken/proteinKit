"""
This script was created on Jan 3rd, 2024 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.
Sorry for my poor English.
"""

import os
from biofkit.proteinKit import pdbKit

# load the information of all amino-acid-residue atoms into a list. output[idx] shows the information of an atom.
def pdb2List(pdbFilePath: str, csvPath: str = '', colName: bool = False) -> list[list[int, str, str, int, str, float, float, float]]:
    """Read a PDB file into a list.

    Args:
        pdbFilePath (str): PDB file path.
        csvPath (str): If you want to output a csv file, this argument is needed.
        colName (str): True for pdbInfoColumns exist in this list, false otherwise.

    Returns:
        list: Protein structure in form of list.
    """

    output: list[list] = []
    # pdbInfoColumns: list = ['Serial', 'Atom', 'ResName', 'ResSeq', 'ChainId', 'X', 'Y', 'Z']
    if (colName):
        proteinKit: pdbKit.ProteinKit = pdbKit.ProteinKit()
        output.append(proteinKit.pdbInfoColumns)
    with open(file=pdbFilePath, mode='r') as pdbFile:
        line: str = pdbFile.readline()
        if (line.startswith('ATOM')):
            output.append([int(line[6:11].strip()), str(line[12:16].strip()), \
                            str(line[17:20].strip()), int(line[22:26].strip()), \
                            str(line[21]), float(line[30:38].strip()), \
                            float(line[38:46].strip()), float(line[46:54])])
        while (line):
            line = pdbFile.readline()
            if (line.startswith('ATOM')):
                output.append([int(line[6:11].strip()), str(line[12:16].strip()), \
                            str(line[17:20].strip()), int(line[22:26].strip()), \
                            str(line[21]), float(line[30:38].strip()), \
                            float(line[38:46].strip()), float(line[46:54])])
    # output the csv file.
    if (csvPath):
        for atomListIdx in range(len(output)):
            output[atomListIdx] = list(map(str, output[atomListIdx]))
        with open(file=csvPath, mode='w') as csvFile:
            csvFile.writelines([','.join(atomInfoList)+'\n' for atomInfoList in output])
    return (output)



# load the information of all amimo-acid-residue atoms into a dictionary, which can be converted into a dataframe with famous `pandas`.
def pdb2Dict(pdbFilePath: str) -> dict[str, list]:
    """Read a PDB file and transfer it into a dictionary.

    Args:
        pdbFilePath (str): PDB file path.

    Returns:
        dict: Protein structure information in form of dictionary.
    """

    # pdbInfoColumns: [str] = ['Serial', 'Atom', 'ResName', 'ResSeq', 'ChainId', 'X', 'Y', 'Z']
    proteinKit: pdbKit.ProteinKit = pdbKit.ProteinKit()
    keyNameList: list[str] = proteinKit.pdbInfoColumns
    output: dict[str, list] = {keyName: [] for keyName in keyNameList}
    with open(file=pdbFilePath, mode='r') as pdbFile:
        line: str = pdbFile.readline()
        if (line.startswith('ATOM')):
            output[keyNameList[0]].append(int(line[6:11].strip()))
            output[keyNameList[1]].append(str(line[12:16].strip()))
            output[keyNameList[2]].append(str(line[17:20].strip()))
            output[keyNameList[3]].append(int(line[22:26].strip()))
            output[keyNameList[4]].append(str(line[21]))
            output[keyNameList[5]].append(float(line[30:38].strip()))
            output[keyNameList[6]].append(float(line[38:46].strip()))
            output[keyNameList[7]].append(float(line[46:54].strip()))
        while (line):
            line = pdbFile.readline()
            if (line.startswith('ATOM')):
                output[keyNameList[0]].append(int(line[6:11].strip()))
                output[keyNameList[1]].append(str(line[12:16].strip()))
                output[keyNameList[2]].append(str(line[17:20].strip()))
                output[keyNameList[3]].append(int(line[22:26].strip()))
                output[keyNameList[4]].append(str(line[21]))
                output[keyNameList[5]].append(float(line[30:38].strip()))
                output[keyNameList[6]].append(float(line[38:46].strip()))
                output[keyNameList[7]].append(float(line[46:54].strip()))
    return (output)


# This funtion 'proDict2ProList' is used to tranfer pdbDict into pdbList.
def proDict2ProList(rawDict: dict) -> list:
    """Transfer a protein structure information with dictionary form into a list. It has been decrepit.

    Args:
        rawDict (dict): Protein structure information with form of dictionary.

    Returns:
        list: Protein structure information with form of list.
    """

    if (proDictIsValid(rawDict)):
        output: list = [[] for count in range(len(rawDict['Serial']))]
        for atomIdx in range(len(rawDict['Serial'])):
            output[atomIdx].append(rawDict['Serial'][atomIdx])
            output[atomIdx].append(rawDict['Atom'][atomIdx])
            output[atomIdx].append(rawDict['ResName'][atomIdx])
            output[atomIdx].append(rawDict['ResSeq'][atomIdx])
            output[atomIdx].append(rawDict['ChainId'][atomIdx])
            output[atomIdx].append(rawDict['X'][atomIdx])
            output[atomIdx].append(rawDict['Y'][atomIdx])
            output[atomIdx].append(rawDict['Z'][atomIdx])
        return (output)
    else:
        return ([])

# This function 'proList2ProDIct' is used to tranfer pdbList into pdbDict.
def proList2ProDict(rawList: list) -> dict:
    """Transfer protein structure in form of list into a dictionary.

    Args:
        rawList (list): Protein structure information in form of list.

    Returns:
        dict: Protein structure information in form of dictionary.
    """

    if (proListIsValid(rawList)):
        output: dict = {'Serial': [], 'Atom': [], 'ResName': [], 'ResSeq': [], \
                    'ChainId': [], 'X': [], 'Y': [], 'Z': []}
        for atomIdx in range(len(rawList)):
            output['Serial'].append(rawList[atomIdx][0])
            output['Atom'].append(rawList[atomIdx][1])
            output['ResName'].append(rawList[atomIdx][2])
            output['ResSeq'].append(rawList[atomIdx][3])
            output['ChainId'].append(rawList[atomIdx][4])
            output['X'].append(rawList[atomIdx][5])
            output['Y'].append(rawList[atomIdx][6])
            output['Z'].append(rawList[atomIdx][7])
        return (output)
    else:
        return ({})


# Define proList is valid or not?
def proListIsValid(proList: list) -> bool:
    """Define protein structure in form of list is valid or not.

    Args:
        proList (list): Protein structure information in form of list.

    Returns:
        bool: True for valid, false otherwise.

    Raises:
        IndexError: Not every sublist has a length of 8.
        TypeError: Atom properties (elements) have wrong type.
    """

    # form and type test
    try:
        for atomInfo in proList:
            # form test:
            if (len(atomInfo) != 8):
                raise (IndexError)
            # type test:
            if ((type(atomInfo[0]) != int) or (type(atomInfo[1]) != str) or \
                (type(atomInfo[2]) != str) or (type(atomInfo[3]) != int) or \
                (type(atomInfo[4]) != str) or (type(atomInfo[5]) != float) or \
                (type(atomInfo[6]) != float) or (type(atomInfo[7]) != float)):
                raise (TypeError)
    except (IndexError):
        print('Uncomplete information for some atom!')
        raise
    except (TypeError):
        print('Atom property type error!')
        raise
    except (Exception):
        raise
    else:
        return (True)


# Define proDict is valid or not?
def proDictIsValid(proDict: dict) -> bool:
    """Define protein structure information in form of dictionary is valid or not.

    Args:
        proDict (dict): Protein structure in form of dictionary.

    Returns:
        bool: True for proDict is valid, false otherwise.

    Raises:
        IndexError: Not each sublist has a length of 8.
        TypeError: Not each atom property (elements) has proper type.
    """

    # form and type test
    try:
        # form test
        if (len(set([len(proDict['Serial']), len(proDict['Atom']), \
                    len(proDict['ResName']), len(proDict['ResSeq']), \
                    len(proDict['ChainId']), len(proDict['X']), \
                    len(proDict['Y']), len(proDict['Z'])])) != 1):
            raise (IndexError)
        # type test
        if (any([type(i) != int for i in proDict['Serial']]) or \
            any([type(i) != str for i in proDict['Atom']]) or \
            any([type(i) != str for i in proDict['ResName']]) or \
            any([type(i) != int for i in proDict['ResSeq']]) or \
            any([type(i) != str for i in proDict['ChainId']]) or \
            any([type(i) != float for i in proDict['X']]) or \
            any([type(i) != float for i in proDict['Y']]) or \
            any([type(i) != float for i in proDict['Z']])):
            raise (TypeError)
    except (IndexError):
        print('Uncomplete information for some atom!')
        raise
    except (TypeError):
        print('Atom property type error!')
        raise
    else:
        return (True)


# This new type 'protein' is for operation pdb file.
class OldProtein:
    """A protein class which has been decrepit.

    This old class organize protein structures in form of list or dictionary.
    This class can contain several chains.

    Attributes:
        name (str): Protein name.
        infoList (list): Protein structure information in form of list.
        infoDict (dict): Protein structure information in form of dictionary.
    """

    def __init__(self, pdbFile: str = '', infoDict: dict = {}, infoList: list = [], chainId: list[str] = ['all']):
        """Inits oldProtein with pdbFile, infoDict, infoList, chainId. Not all of them are necessary.
        Args:
            pdbFile (str): The path of your PDB file.
            infoDict (dict): Organized dictionary of protein structure information.
            infoList (list): Organized list of protein structure information.
            chainId (list): Read listed chains and ignore the others.
        """

        # PDB file as INPUT
        if (pdbFile):
            self.name: str = pdbFile.split(os.sep)[-1].rstrip('.pdb')
            self.infoList: list = pdb2List(pdbFile)
            _ = proListIsValid(self.infoList)
            self.infoDict: dict = pdb2Dict(pdbFilePath=pdbFile)
        # proDict as INPUT
        elif (infoDict):
            self.name: str = 'unnamed'
            _ = proDictIsValid(infoDict)
            self.infoDict: dict = infoDict
            self.infoList: list = proDict2ProList(infoDict)
        # proList as INPUT
        elif (infoList):
            self.name: str = 'unnamed'
            _ = proListIsValid(infoList)
            self.infoList: list = infoList
            self.infoDict: dict = proList2ProDict(infoList)
        # Remove unused atom
        if (chainId != ['all']):
            for atom in self.infoList:
                if (atom[4] not in chainId):
                    self.infoList.remove(atom)

    def __str__(self) -> str:
        return ('{name}: containing Chain {chainList}'.format(name = self.name, chainList = list(set(self.infoDict['ChainId']))))

    def __repr__(self) -> str:
        return ('{name}: containing Chain {chainList}'.format(name = self.name, chainList = list(set(self.infoDict['ChainId']))))

    def __sizeof__(self) -> int:
        return (len(set(self.infoDict['ChainId'])))

    def __len__(self) -> int:
        return (len(set(self.infoDict['ChainId'])))



