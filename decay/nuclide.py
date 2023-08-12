from multiprocessing.sharedctypes import Value
import re

class Nuclide:

    name_regex = r"([A-Za-z]+)[-]*(\d+)[_]*(\w*)"
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

    element = {i + 1: name for i, name in enumerate(elementnames)}
    element = element | {name: i + 1 for i, name in enumerate(elementnames)}
    element = element | {name.upper(): i + 1 for i, name in enumerate(elementnames)}
    isomernames = ["", "M", "N", "O"]
    isomer = {0: "", 1: "M", 2: "N", 3:"O", "M": 1, "N": 2, "G": 0, "O": 3, "": 0}
    isomer = {k: i for i, k in enumerate(isomernames)}

    @classmethod
    def getA(cls, isot: str | int) -> int:
        """Return the mass number of the nuclide.

        Args:
            isot (str | int): A isotope's name or an isotope's ZAM

        Returns:
            int: The isotope's mass number
        """
        if isinstance(isot, str):
            if "_" not in isot:
                regex = re.search(r"\D+(\d+)\D*", isot.upper())
                try:
                    return int(regex.group(1))
                except AttributeError:
                    raise ValueError(f"No mass number found for {isot}")
            else:
                _name = isot.split("_")[0]
                regex = re.search(r"\D+(\d+)", _name.upper())
                try:
                    return int(regex.group(1))
                except AttributeError:
                    raise ValueError(f"No mass number found for {isot}")
        elif isinstance(isot, int):
            assert isot >= 10_000
            return isot % 10_000 // 10
        else:
            raise ValueError("Input should be an isotope name or id.", isot)

    @classmethod
    def getZ(cls, isot: str | int) -> int:
        if isinstance(isot, str):
            if "_" not in isot:
                regex = re.search(r"(\D+)\d+\D*", isot.upper())
                try:
                    return cls.element[regex.group(1)]
                except KeyError:
                    raise KeyError(f"{isot} element does not exist.")
            else:
                _name = isot.split("_")[0]
                regex = re.search(r"(\D+)\d+", _name.upper())
                try:
                    return cls.element[regex.group(1)]
                except KeyError:
                    raise KeyError(f"{isot} element does not exist.")
        elif isinstance(isot, int):
            return isot // 10_000

    @classmethod
    def getE(cls, isot: str | int) -> int:
        if isinstance(isot, str):
            if "_" not in isot:
                regex = re.search(r"\D+\d+(\D*)", isot.upper())
                try:
                    return cls.isomer[regex.group(1)]
                except KeyError:
                    raise KeyError(f"{isot} has unreadable isomeric value.")
            else:
                try:
                    _name = isot.split("_")[1]
                    return int(_name[1:])
                except IndexError:
                    return 0
        elif isinstance(isot, int):
            return isot % 10

    @classmethod
    def getZAE(cls, isot: str | int) -> int:
        Z, A, E = cls.getZ(isot), cls.getA(isot), cls.getE(isot)
        return 10_000 * Z + 10 * A + E

    @classmethod
    def getName(cls, zae: int) -> str:
        Z, A, E = cls.getZ(zae), cls.getA(zae), cls.getE(zae)
        try:
            elem = cls.element[Z].capitalize()
        except KeyError:
            raise KeyError(f"Z id in {zae} does not map to any element.")

        try:
            isom = cls.isomer[E]
        except KeyError:
            raise KeyError(f"E id in {zae} does not map to any isomeric state.")

        return f"{elem}{A:d}{isom}"

    def __init__(self, initiator: str | int) -> None:
        self.ZAE = self.getZAE(initiator)

    @property
    def ZAE(self):
        return self._ZAE

    @ZAE.setter
    def ZAE(self, rhs):
        self._ZAE = rhs
        self._name = self.getName(rhs)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, rhs):
        self.ZAE = self.getZAE(rhs)
    
    @property
    def Z(self):
        return self.getZ(self.ZAE)

    @property
    def A(self):
        return self.getA(self.ZAE)

    @property
    def E(self):
        return self.getE(self.ZAE)
        
    

