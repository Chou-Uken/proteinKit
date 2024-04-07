'''
This script was created on Dec 17th, 2023 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.
Sorry for my poor English.
'''

import os
from biofkit.proteinKit.proteinClass import Atom, Residue, Peptide, Protein

class ProteinKit:

    # Dictionary used to transfer abbreviation with the first letter capitalized to shorter abbreviation.
    aaDictTHREE2One: dict[str, str] = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',    \
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',    \
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',    \
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',    \
        'SEC': 'U', 'PYL': 'O'                                          # Rare amino acid
    }

    aaDictThree2One: dict[str, str] = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',    \
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',    \
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',    \
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',    \
        'Sec': 'U', 'Pyl': 'O'                                          # Rare amino acid
    }

    # Dictionary used to transfer abbreviation to shorter abbreviation.
    aaDictthree2One: dict[str, str] = {
        'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'cys': 'C',    \
        'gln': 'Q', 'glu': 'E', 'gly': 'G', 'his': 'H', 'ile': 'I',    \
        'leu': 'L', 'lys': 'K', 'met': 'M', 'phe': 'F', 'pro': 'P',    \
        'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V',    \
        'sec': 'U', 'pyl': 'O'                                          # Rare amino acid
    }

    # Dictionary used to transfer shorter abbreviation to abbreviation with the first letter capitalized.
    aaDictOne2Three: dict[str, str] = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',    \
        'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',    \
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',    \
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',    \
        'U': 'Sec', 'O': 'Pyl'                                          # Rare amino acid
    }

    # Dictionary used to transfer shorter abbreviation to abbreviation.
    aaDictOne2three: dict[str, str] = {
        'A': 'ala', 'R': 'arg', 'N': 'asn', 'D': 'asp', 'C': 'cys',    \
        'Q': 'gln', 'E': 'glu', 'G': 'gly', 'H': 'his', 'I': 'ile',    \
        'L': 'leu', 'K': 'lys', 'M': 'met', 'F': 'phe', 'P': 'pro',    \
        'S': 'ser', 'T': 'thr', 'W': 'trp', 'Y': 'tyr', 'V': 'val',    \
        'U': 'sec', 'O': 'pyl'                                          # Rare amino acid
    }

    # PDB file infos import.
    # Column meanings:
    # Serial: Atom serial number.
    # Atom: Atom name.
    # ResName: Residue name.
    # ResSeq: Residue sequence number.
    # ChainId: Chain identifier.
    # X: Orthogonal coordinates for X in angstroms (A).
    # Y: Orthogonal coordinates for Y in angstroms (A).
    # Z: Orthogonal coordinates for Z in angstroms (A).
    pdbInfoColumns: list[str] = ['Serial', 'Atom', 'ResName', 'ResSeq', 'ChainId', 'X', 'Y', 'Z']


    def __init__(self):
        pass



# transfer protein structure file (pdb) into sequence string.
def pdb2Seq(pdbFilePath: str, fastaFilePath: str = '', fastaLineLen: int = 80) -> dict[str, str]:
    """Extract sequence information from a PDB file.

    Args:
        pdbFilePath (str): Path of PDB file.
        fastaFilePath (str): FASTA file file output. None for nothing to create.
        fastaLineLen (int): How many residue letter each line.

    Returns:
        dict[str, str]: Key is chainId. Value is sequence.
    """

    proteinKit: ProteinKit = ProteinKit()
    with open(file=pdbFilePath, mode='r') as pdbFile:
        thisChainId: str = 'defined'
        line: str = pdbFile.readline()
        chainSeq: str = ''
        resSeq: int = 0
        resName: str = ''
        output: dict[str, str] = {}

        if (line.startswith('ATOM')):
            thisChainId = line[21]
            resSeq = int(line[22:26].strip())
            resName = line[17:20]
            if (resName == 'UNK'):
                chainSeq += '-'
            else:
                chainSeq += proteinKit.aaDictTHREE2One[resName]

        while (line):
            line = pdbFile.readline()
            if (line.startswith('ATOM')):
                if (resSeq == 0):                                       # chain begin
                    thisChainId = line[21]
                    resSeq = int(line[22:26].strip())
                    resName = line[17:20]
                    if (resName == 'UNK'):
                        chainSeq += '-'
                    else:
                        chainSeq += proteinKit.aaDictTHREE2One[resName]
                elif (int(line[22:26].strip()) == resSeq):
                    continue
                else:
                    if (int(line[22:26].strip()) == resSeq + 1):
                        resSeq = int(line[22:26].strip())
                        resName = line[17:20]
                        if (resName == 'UNK'):
                            chainSeq += '-'
                        else:
                            chainSeq += proteinKit.aaDictTHREE2One[resName]
                    else:
                        gap: int = int(line[22:26].strip()) - resSeq - 1
                        gapSeq: str = '-' * gap
                        chainSeq += gapSeq
                        resSeq = int(line[22:26].strip())
                        resName = line[17:20]
                        if (resName == 'UNK'):
                            chainSeq += '-'
                        else:
                            chainSeq += proteinKit.aaDictTHREE2One[resName]
            elif (line.startswith('TER')):                              # recognizing termination with the line 'TER'
                output[thisChainId] = chainSeq
                resSeq = 0
                chainSeq = ''
    # output the fasta files.
    if (fastaFilePath):
        fileName: str = pdbFilePath.split(os.sep)[-1].rstrip('.pdb')
        with open(file=os.path.join(fastaFilePath, fileName+'.fasta'), mode='w') as fastaFile:
            for key in output.keys():
                fastaFile.write('>'+fileName+'_chain_'+key+'\n')
                thisLine: list[str] = [output[key][i:i+fastaLineLen] for i in range(0, len(output[key]), fastaLineLen)]
                for i in thisLine:
                    fastaFile.write(i + '\n')
            print(fileName+' converted!\n')
    return (output)




# load PDB file
def readPDB(pdbFile: str) -> Protein:
    """Load PDB file to get a Protein.
    Args:
        pdbFile (str): PDB file path.

    Returns:
        Protein: Class Protein in proteinClass.py
    """

    proteinKit: ProteinKit = ProteinKit()
    with open(file=pdbFile, mode='r') as pdb:
        atomBuffer: list[Atom] = []
        residueBuffer: list[Residue] = []
        peptideBuffer: list[Peptide] = []
        resSeq: int = 0
        resName: str = ''
        chainId: str = ''
        line: str = pdb.readline()
        if (line.startswith('ATOM')):
            atomBuffer.append(Atom(serial=int(line[6:11].strip()), atom=str(line[12:16].strip()), x=float(line[30:38].strip()), y=float(line[38:46].strip()), z=float(line[46:54].strip())))
            resName = str(line[17:20].strip())
            resSeq = int(line[22:26].strip())
            chainId = str(line[21])
        while (line):
            line = pdb.readline()
            if (line.startswith('ATOM')):
                if ((not chainId) or (str(line[21])==chainId)):
                    # Chain continue
                    if ((not resSeq) or (int(line[22:26].strip())==resSeq)):
                        # Residue continue
                        resName = str(line[17:20].strip())
                        resSeq = int(line[22:26].strip())
                        chainId = str(line[21])
                        atomBuffer.append(Atom(serial=int(line[6:11].strip()), atom=str(line[12:16].strip()), x=float(line[30:38].strip()), y=float(line[38:46].strip()), z=float(line[46:54])))
                    else:
                        # Residue TER
                        try:
                            residueBuffer.append(Residue(atomList=atomBuffer, resSeq=resSeq, resName=proteinKit.aaDictTHREE2One[resName]))
                        except (KeyError):
                            residueBuffer.append(Residue(atomList=atomBuffer, resSeq=resSeq, resName=resName))
                        except (Exception) as e:
                            print(e)
                            raise
                        # Initiate atomBuffer
                        atomBuffer = []
                        # Read new atom properties
                        atomBuffer.append(Atom(serial=int(line[6:11].strip()), atom=str(line[12:16].strip()), x=float(line[30:38].strip()), y=float(line[38:46].strip()), z=float(line[46:54].strip())))
                        # Read new residue/chain properties
                        resSeq = int(line[22:26].strip())
                        resName = str(line[17:20].strip())
                        chainId = str(line[21])
                else:
                    # Chain TER/Residue TER
                    try:
                        residueBuffer.append(Residue(atomList=atomBuffer, resSeq=resSeq, resName=proteinKit.aaDictTHREE2One[resName]))
                    except (KeyError):
                        residueBuffer.append(Residue(atomList=atomBuffer, resSeq=resSeq, resName=resName))
                    except (Exception) as e:
                        print(e)
                        raise
                    peptideBuffer.append(Peptide(resList=residueBuffer, chainId=chainId))
                    # Initiate atomBuffer/residueBuffer
                    atomBuffer = []
                    residueBuffer = []
                    # Read new atom properties
                    atomBuffer.append(Atom(serial=int(line[6:11].strip()), atom=str(line[12:16].strip()), x=float(line[30:38].strip()), y=float(line[38:46].strip()), z=float(line[46:54].strip())))
                    # Read new residue/chain properties
                    resSeq = int(line[22:26].strip())
                    resName = str(line[17:20].strip())
                    chainId = str(line[21])
        if (not atomBuffer):
            residueBuffer.append(Residue(atomList=atomBuffer, resSeq=resSeq, resName=proteinKit.aaDictTHREE2One[resName]))
            peptideBuffer.append(Peptide(resList=residueBuffer, chainId=chainId))
            protein: Protein = Protein(pepList=peptideBuffer, proteinName=pdbFile.split(os.sep)[-1].rstrip('.pdb'))
        else:
            try:
                peptideBuffer.append(Peptide(resList=residueBuffer, chainId=chainId))
            except Exception:
                pass
            finally:
                protein: Protein = Protein(pepList=peptideBuffer, proteinName=pdbFile.split(os.sep)[-1].rstrip('.pdb'))
        return (protein)


# Merge Peptides to generate a Protein.
def assemble(pepList: list[Peptide], name: str = 'Unnamed') -> Protein:
    """merge peptides.

    Args:
        newProtein (list[Peptide]): List of merged peptides.
        name (str): Name of new protein.

    Returns:
        Protein.
    """

    if (type(pepList) != list):
        pepList = [pepList]
    newProtein: Protein = Protein(pepList=pepList, proteinName=name)
    return (newProtein)


# Write Protein into a PDB file.
def writePDB(protein: list[Protein], filePath: str, printInfo: bool = True) -> None:
    """Write PDB file
    
    Args:
        protein (Protein): The Protein you want to save as PDB file.
        filePath (str): Output path.
        printInfo (bool): Whether printing writing output or not.
    """

    if (type(protein) != list):
        protein = [protein]
    atomCol: str = 'ATOM'
    col5_6: str = ' ' * 2
    col12: str = ' '
    col17: str = ' '
    col21: str = ' '
    col27_30: str = ' ' * 4    
    col67_76: str = ' ' * 10
    col79_80: str = ' ' * 2
    for eachProtein in protein:
        with open(file=os.path.join(filePath, eachProtein.name+'.pdb'), mode='w') as thisPDB:
            for peptide in eachProtein.pepSet:
                for residue in peptide.resSet:
                    for atom in residue.atomSet:
                        serialCol: str = ' ' * (5 - len(str(atom.serial))) + str(atom.serial).strip()
                        nameCol: str = (' ' + atom.atom) if (len(atom.atom) < 4) else atom.atom
                        resNameCol: str = residue.getName()
                        chainIDCol: str = peptide.getChainId()[0]
                        resSeqCol: str = ' ' * (4 - len(str(residue.resSeq))) + str(residue.resSeq)
                        xCol: str = ' ' * (8 - len('%.3f'%atom.getCoord()['x'])) + '%.3f'%atom.getCoord()['x']
                        yCol: str = ' ' * (8 - len('%.3f'%atom.getCoord()['y'])) + '%.3f'%atom.getCoord()['y']
                        zCol: str = ' ' * (8 - len('%.3f'%atom.getCoord()['z'])) + '%.3f'%atom.getCoord()['z']
                        occupancyCol: str = ' ' * (6 - len('%.2f'%atom.occupancy)) + '%.2f'%atom.occupancy
                        tempFactorCol: str = ' ' * (6 - len('%.2f'%atom.tempFactor)) + '%.2f'%atom.tempFactor
                        elementCol: str = ' ' * (2 - len(atom.element)) + str(atom.element)
                        line: str = ''.join([atomCol, col5_6, serialCol, col12, nameCol, col17, resNameCol, col21, chainIDCol, resSeqCol, col27_30, xCol, yCol, zCol, occupancyCol, tempFactorCol, col67_76, elementCol, col79_80])
                        thisPDB.write(line + '\n')
                thisPDB.write('TER\n')
            thisPDB.write('END\n')
