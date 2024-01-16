import pymatgen
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import *


def get_kpoints_SC(structure: pymatgen.core.structure.Structure):
    kpoints = KPathSetyawanCurtarolo(structure).kpath
    points = []
    full_path = []
    for path in kpoints["path"]:
        full_path = full_path + path
    for point in full_path:
        points.append([point, list(kpoints['kpoints'][point])])
    return points

def get_kpoints_McQueen(structure):

    sym = SpacegroupAnalyzer(structure)
    lattice_type = sym.get_lattice_type()
    spg_symbol = sym.get_space_group_symbol()
    kpath = []
    if lattice_type == "cubic":
        if "P" in spg_symbol:
            kpath = ["\\Gamma", "M", "x", "\\Gamma", "R"]
        elif "F" in spg_symbol:
            kpath = ["L", "\\Gamma", "X", "W", "L", "K", "\\Gamma", "U"]
        elif "I" in spg_symbol:
            kpath = ["\\Gamma", "N", "P", "\\Gamma", "H", "N"]
        else:
            warn(f"Unexpected value for {spg_symbol=}")
    else:
        warn("non-cubic lattice type not implemented, using SC")
        return get_kpoints_SC(structure)

    kpoints = KPathSetyawanCurtarolo(structure).kpath
    points = []
    for point in kpath:
        points.append([point, list(kpoints['kpoints'][point])])
    return points
'''
    elif lattice_type == "tetragonal":
        if "P" in spg_symbol:
            self._kpath = self.tet()
        elif "I" in spg_symbol:
            a = self._conv.lattice.abc[0]
            c = self._conv.lattice.abc[2]
            if c < a:
                self._kpath = self.bctet1(c, a)
            else:
                self._kpath = self.bctet2(c, a)
        else:
            warn(f"Unexpected value for {spg_symbol=}")

    elif lattice_type == "orthorhombic":
        a = self._conv.lattice.abc[0]
        b = self._conv.lattice.abc[1]
        c = self._conv.lattice.abc[2]

        if "P" in spg_symbol:
            self._kpath = self.orc()

        elif "F" in spg_symbol:
            if 1 / a ** 2 > 1 / b ** 2 + 1 / c ** 2:
                self._kpath = self.orcf1(a, b, c)
            elif 1 / a ** 2 < 1 / b ** 2 + 1 / c ** 2:
                self._kpath = self.orcf2(a, b, c)
            else:
                self._kpath = self.orcf3(a, b, c)

        elif "I" in spg_symbol:
            self._kpath = self.orci(a, b, c)

        elif "C" in spg_symbol or "A" in spg_symbol:
            self._kpath = self.orcc(a, b, c)
        else:
            warn(f"Unexpected value for {spg_symbol=}")

    elif lattice_type == "hexagonal":
        self._kpath = self.hex()

    elif lattice_type == "rhombohedral":
        alpha = self._prim.lattice.parameters[3]
        if alpha < 90:
            self._kpath = self.rhl1(alpha * pi / 180)
        else:
            self._kpath = self.rhl2(alpha * pi / 180)

    elif lattice_type == "monoclinic":
        a, b, c = self._conv.lattice.abc
        alpha = self._conv.lattice.parameters[3]
        # beta = self._conv.lattice.parameters[4]

        if "P" in spg_symbol:
            self._kpath = self.mcl(b, c, alpha * pi / 180)

        elif "C" in spg_symbol:
            kgamma = self._rec_lattice.parameters[5]
            if kgamma > 90:
                self._kpath = self.mclc1(a, b, c, alpha * pi / 180)
            if kgamma == 90:
                self._kpath = self.mclc2(a, b, c, alpha * pi / 180)
            if kgamma < 90:
                if b * cos(alpha * pi / 180) / c + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2 < 1:
                    self._kpath = self.mclc3(a, b, c, alpha * pi / 180)
                if b * cos(alpha * pi / 180) / c + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2 == 1:
                    self._kpath = self.mclc4(a, b, c, alpha * pi / 180)
                if b * cos(alpha * pi / 180) / c + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2 > 1:
                    self._kpath = self.mclc5(a, b, c, alpha * pi / 180)
        else:
            warn(f"Unexpected value for {spg_symbol=}")

    elif lattice_type == "triclinic":
        kalpha = self._rec_lattice.parameters[3]
        kbeta = self._rec_lattice.parameters[4]
        kgamma = self._rec_lattice.parameters[5]
        if kalpha > 90 and kbeta > 90 and kgamma > 90:
            self._kpath = self.tria()
        if kalpha < 90 and kbeta < 90 and kgamma < 90:
            self._kpath = self.trib()
        if kalpha > 90 and kbeta > 90 and kgamma == 90:
            self._kpath = self.tria()
        if kalpha < 90 and kbeta < 90 and kgamma == 90:
            self._kpath = self.trib()

    else:
        warn(f"Unknown lattice type {lattice_type}")
        '''
