import re
from collections import namedtuple


class Nuclide:
    """A class that represents nuclides."""

    name_regex = r"([A-Za-z]+)[-_]*(\d+)[_]*(\w*)"
    elementnames = [
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "I",
        "Xe",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Rn",
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U",
        "Np",
        "Pu",
        "Am",
        "Cm",
        "Bk",
        "Cf",
        "Es",
        "Fm",
        "Md",
        "No",
        "Lr",
        "Rf",
        "Db",
        "Sg",
        "Bh",
        "Hs",
        "Mt",
        "Ds",
        "Rg",
        "Cn",
        "Uut",
        "Fl",
        "Uup",
        "Lv",
    ]

    elements_dict = {i + 1: name for i, name in enumerate(elementnames)}
    elements_dict = elements_dict | {name: i + 1 for i, name in enumerate(elementnames)}
    elements_dict = elements_dict | {
        name.upper(): i + 1 for i, name in enumerate(elementnames)
    }
    isomernames = ["", "M", "N", "O"]
    isomer = {0: "", 1: "M", 2: "N", 3: "O"}
    isomer |= {"": 0, "G": 0, "M": 1, "N": 2, "O": 3}
    isomer |= {"": 0, "g": 0, "m": 1, "n": 2, "o": 3}
    isomer |= {f"m{i}": i for i in range(1, 10)}

    @classmethod
    def slice_nucid(cls, nucid: int) -> tuple[int, int, int]:
        """Turns a nucid into a named tuple (Z, A, M).

        Args:
            nucid (int): A nuclide id.

        Returns:
            tuple[int, int, int]: A named zam tuple (z, a, m)

        """
        e = nucid % 10
        a = nucid % 10_000 // 10
        z = nucid // 10_000
        return namedtuple("Zam", ["z", "a", "m"])(z, a, e)

    @classmethod
    def aggregate_zam(cls, zam: tuple[int, int, int]) -> int:
        """Turns a (Z, A, M) tuple into a nuclide id.

        Args:
            zam (tuple[int, int, int]): A tuple (z, a, m)

        Returns:
            int: A nuclide id.

        """
        return 10_000 * zam[0] + 10 * zam[1] + zam[2]

    @classmethod
    def name_to_sliced_zam(cls, name: str) -> int:
        """Compute a zam tuple from a nuclide name.

        Args:
            name (str): The nuclide's name.

        Returns:
            tuple[int, int, int]: A named zam tuple (z, a, m)

        """
        groups = re.search(cls.name_regex, name).groups()
        z = cls.elements_dict[groups[0].capitalize()]
        a = int(groups[1])
        e = cls.isomer[groups[2]]
        return namedtuple("Zam", ["z", "a", "e"])(z, a, e)

    @classmethod
    def zam_to_name(cls, zam: tuple[int, int, int]) -> str:
        """Compute a nuclide's name from its (z, a, m) tuple.
        
        Args:
            zam (tuple[int, int, int]): A nuclide's (z, a ,m) tuple.

        Returns:
            str: The nuclide's name.
        """
        elem = cls.elements_dict[zam[0]]
        A = zam[1]
        isom = cls.isomer[zam[2]]
        return f"{elem}{A:d}{isom}"

    def __init__(self, initiator: str | int) -> None:
        if isinstance(initiator, str):
            self.zam = self.name_to_sliced_zam(initiator)
        elif isinstance(initiator, int):
            self.zam = self.slice_nucid(initiator)

    @property
    def zam(self):
        return self._zam

    @zam.setter
    def zam(self, rhs):
        self._zam = namedtuple("Zam", ["z", "a", "m"])(*rhs)
        self._name = self.zam_to_name(self._zam)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, rhs):
        self.zam = self.name_to_sliced_zam(rhs)

    @property
    def nucid(self):
        return self.aggregate_zam(self.zam)

    @nucid.setter
    def nucid(self, rhs):
        self.zam = self.slice_nucid(rhs)

    @property
    def z(self):
        return self.zam.z

    @property
    def a(self):
        return self.zam.a

    @property
    def m(self):
        return self.zam.m

    @property
    def symb(self):
        return self.elements_dict[self.z]

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.zam.__repr__()
