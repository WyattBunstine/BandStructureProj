import mp_api
from mp_api.client import MPRester
import pymatgen
from pymatgen.core import *
import json
import pandas as pd
import itertools

HtoEv = 27.2114 # 1 Hartree is 27 Electron volts

def search_mp_element_combos(elements: list[str], num_combos=2, print_mats=False):
    data = []
    for combo in itertools.combinations(elements, num_combos):
        API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
        with MPRester(API_KEY) as mpr:
            data += mpr.materials.summary.search(elements=list(combo), is_stable=True, all_fields=False,
                                                 fields=["composition", "material_id", "structure"])
    if print_mats:
        for mat in data:
            print(str(mat.composition) + "  " + str(mat.material_id))
    return {str(mat.material_id): mat.structure for mat in data}


def get_dftu_information(mat: Structure):
    dftu_string = ""
    species = []
    for el in mat.sites:
        for element in el.species:
            if element not in species:
                species.append(element)
    i = 0
    for specie in sorted(species):
        i += 1
        if specie.is_transition_metal:
            partial_orbitals = []
            for orbital in specie.full_electronic_structure:
                if orbital[1] == "s" and orbital[2] < 2:
                    partial_orbitals.append(0)
                elif orbital[1] == "p" and orbital[2] < 6:
                    partial_orbitals.append(1)
                elif orbital[1] == "d" and orbital[2] < 10:
                    partial_orbitals.append(2)
                elif orbital[1] == "f" and orbital[2] < 14:
                    partial_orbitals.append(3)
            for partial_orbital in partial_orbitals:
                a = (specie.ionization_energy - specie.electron_affinity)/HtoEv
                dftu_string += ("  " + str(i) + " " + str(partial_orbital) + " "
                                + f'{a:.10f}' + "\n")
    return dftu_string


def gen_cfg(matid, mat):
    return {
        "MatType": "MP",
        "MatLoc": "data/TIProj/" + matid + "/",
        "MatID": matid,
        "tasks": ["GS", "BS", "DOS"],
        "kpoints": "MP",
        "overwrite": True,
        "elkparams": {
            "primcell": ".true.",
            "vhighq": ".false.",
            "highq": ".false.",
            "trimvg": ".true.",
            "maxscl": 100,
            "mixtype": 1,
            "xctype": 22,
            "isgkmax": -1,
            "ngridk": "4 4 4",
            "bfieldc": "0.0 0.0 -0.01",
            "batch": ".true.",
            "spinorb": ".true.",
            "spinpol": ".true.",
            "nempty": 10,
            "stype": 3,
            "swidth": 0.0001,
            "dft+u": "1 5\n" + get_dftu_information(mat)
        },
        "SGU": False,
        "ROUNDING": 4,
        "numBandPoints": 1000,
        "remote": True,
        "monitor": False
    }


def main():
    mats = search_mp_element_combos(["Ni", "Cr"])
    for matid in mats.keys():
        print(gen_cfg(matid, mats[matid]))


def get_configs(elements: list[str]):
    mats = search_mp_element_combos(elements)
    configs = []
    for matid in mats.keys():
        configs.append(gen_cfg(matid, mats[matid]))
    return configs


if __name__ == "__main__":
    main()
