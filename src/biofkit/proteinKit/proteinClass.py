# This script was created on Feb 22nd, 2024 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
# Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.

class Atom(object):
    """A class containing number, type and coordinate of an atom.

    Attributes:
        serial (int): Id of atom.
        atom (str): Type of atom.
        x (float): Coordinate of x.
        y (float): Coordinate of y.
        z (float): Coordinate of z.
    """

    def __init__(self, serial: int, atom: str, x: float, y: float, z: float):
        """Inits Atom with its id, type and coordinate.

        Args:
            serial (int): Id of atom.
            atom (str): Type of atom.
            x (float): Coordinate of x.
            y (float): Coordinate of y.
            z (float): Coordinate of z.
        """

        self.serial: int = int(serial)
        self.atom: str = str(atom)
        self.x: float = float(x)
        self.y: float = float(y)
        self.z: float = float(z)

    def __str__(self) -> str:
        return ('Atom: {name} at location of ({x}, {y}, {z})'.format(name=self.atom, x=self.x, y=self.y, z=self.z))

    def __repr__(self) -> str:
        return ('Atom: {name} at location of ({x}, {y}, {z}))'.format(name=self.atom, x=self.x, y=self.y, z=self.z))

    def getCoord(self) -> dict[str, float]:
        """Get coordinate of some atom.

        Returns:
            dict[str, float]: x/y/z axis and its coordinate.
        """

        return ({'x': self.x, 'y': self.y, 'z': self.z})


class Residue(object):
    """A class containing number, type and its Atoms of a residue in a peptide.

    Attributes:
        resSeq (int): Id of residue.
        resName (str): Type of residue.
        atomSet (list[Atom]): Set of Atom belonging to this residue.
    """

    def __init__(self, atomList: list[Atom], resSeq: int = 0, resName: str = '*'):
        """Inits residue with its id, type and atoms.

        Args:
            resSeq (int): Id of residue.
            resName (str): Type of residue.
            atomSet (list[Atom]): Set of Atoms belonging to this residue.
        """

        self.resSeq: int = int(resSeq)
        self.resName: str = str(resName)
        # type test
        try:
            for atom in atomList:
                if (type(atom) != Atom):
                    raise (TypeError)
        except (TypeError):
            print('TypeError! atomList must only contain elements of Atom type.')
            raise
        except (Exception):
            raise
        else:
            self.atomSet: list[Atom] = atomList

    def __str__(self) -> str:
        return (self.resName + str(self.resSeq))

    def __repr__(self) -> str:
        return (self.resName + str(self.resSeq))

    def __getitem__(self, idx: int) -> Atom:
        return (self.atomSet[idx])

    def getName(self) -> str:
        """To get type of this residue.

        Returns:
            str: Type of this residue.
        """

        return (self.resName)


class Peptide(object):
    """a class of peptide with its id and set of residues belonging to it.

    Attributes:
        chainId (str): Id of this chain.
        resSet (list[Residue]): Set of Residues belonging to this chain.
    """

    def __init__(self, resList: list[Residue], chainId: str = 'A'):
        """Inits Peptide with its id and residues belonging to it.

        Args:
            chainId (str): Id of this chain.
            resSet (list[Residue]): Set of Residues belonging to this chain.
        """

        self.chainId: str = chainId
        # type test
        try:
            for res in resList:
                if (type(res) != Residue):
                    raise (TypeError)
        except (TypeError):
            print('TypeError! resList must only contain elements of Residue type.')
            raise
        except (Exception):
            raise
        else:
            self.resSet: list[Residue] = resList

    def __str__(self) -> str:
        return ('Chain {name}: {seq}'.format(name=self.chainId, seq=''.join([i.getName() for i in self.resSet])))

    def __repr__(self) -> str:
        return ('Chain {name}: {seq}'.format(name=self.chainId, seq=''.join([i.getName() for i in self.resSet])))

    def __len__(self) -> int:
        return (len(self.resSet))

    def __getitem__(self, idx: int) -> Residue:
        return (self.resSet[idx])

    def getChainId(self) -> str:
        """To get id of this chain.

        Returns:
            str: Id of this chain.
        """

        return (self.chainId)

    def cut(self, end: int, beginning: int = 0, newChainId: str = 'X'):
        """To cut a peptide at specific site.

        Args:
            end (int): Id of the residue at the end of output.
            beginning (int): Id of the residue at the beginning of output.
            newChainId (str): Id of the Peptide you created.

        Returns:
            Peptide: A new peptide.
        """

        if ((beginning == 0) and (end == len(self.resSet))):
            print('An identical peptide is created!')
        return (Peptide(resList=self.resSet[beginning: end], chainId=newChainId))



class Protein(object):
    """a class containing protein with its peptides.

    Attributes:
        name (str): Name of this protein.
        pepSet (list[Peptide]): The list of Peptides belonging to the Protein.
    """

    def __init__(self, pepList: list[Peptide], proteinName: str = 'Unnamed'):
        """Inits Protein with its Peptides and its name.

        Args:
            pepList (list[Peptide]): a list comprising of Peptides.
            name (str): Name of this protein.
        """

        self.name = proteinName
        # type test
        try:
            for pep in pepList:
                if (type(pep) != Peptide):
                    raise (TypeError)
        except (TypeError):
            print('TypeError! PepList must only contain elements of Peptide type.')
            raise
        except (Exception):
            raise
        else:
            self.pepSet: list[Peptide] = pepList

    def __str__(self) -> str:
        return (self.name + ': ' + str(self.pepSet))

    def __repr__(self) -> str:
        return (self.name + ': ' + str(self.pepSet))

    def __getitem__(self, idx: int) -> Peptide:
        return (self.pepSet[idx])

    def pick(self, chainId: list[str]) -> list[Peptide]:
        """pick specific peptide from a Protein.

        Args:
            chainId (list[str]): Id of selected peptides.

        Returns:
            list[Peptide]: a list of Peptides with given Id.
        """

        output: list[Peptide] = [peptide for peptide in self.pepSet if (peptide.getChainId() in chainId)]
        return (output)

