import sys
import os
from warnings import *
import json
import pymatgen
from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


ROUNDING = 4
numBandPoints = 400
BtoA = 0.52917721090  # 1 bohr radius is 0.52 angstroms

def get_mp_mat(MPID):
    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    # inputs
    with MPRester(API_KEY) as mpr:
        # save the relevent data to a dict
        return mpr.get_structure_by_material_id((MPID))

def gen_spacegroupin(path, mat):

    file = open(path + "/spacegroup.in", 'w')
    file.write("'" + mat.get_space_group_info()[0] + "'\n")
    latparams = list(mat.lattice.parameters)
    for param in latparams:
        latparams.append(round(latparams[0], ROUNDING))
        latparams.remove(latparams[0])
    file.write(str(latparams[0] / BtoA) + " " + str(latparams[1] / BtoA) + " " + str(latparams[2] / BtoA) + "\n")
    file.write(str(latparams[3]) + " " + str(latparams[4]) + " " + str(latparams[5]) + "\n")
    file.write("1 1 1\n.false.\n")

    elements = dict()
    for el in mat.sites:
        if not el.species in elements.keys():
            elements[el.species] = 1
        else:
            elements[el.species] += 1

    file.write(str(len(elements.keys()))+ "\n")

    for el in elements.keys():
        file.write("'" + str(el.elements[0].symbol) + "' '" + str(el.elements[0].symbol) + ".in'\n")
        file.write(str(elements[el]) + "\n")
        for site in mat.sites:
            if site.species == el:
                file.write(str(round(float(site.a), ROUNDING)) + " " + str(round(float(site.b), ROUNDING)) + " "
                           + str(round(float(site.c), ROUNDING)) + "\n")
    file.close()
    os.system("cd " + path + " && spacegroup")

def main():
    #check if the config file is given
    if len(sys.argv) < 2:
        warn("no args given")
        return -1
    #check if config exists
    config = sys.argv[1]
    if not os.path.exists(config):
        warn("config file, "+str(config)+", does not exist")
        return -1
    with open(config) as f:
        config = json.load(f)
    #get the material data, from materials proj or cif
    mat = 0
    if config["MatType"] == "MP":
        mat = get_mp_mat(config["MatLoc"])
    elif config["MatType"] == "cif":
        if not os.path.exists(config["MatLoc"]):
            warn("specified cif does not exist")
            return -1
        mat = pymatgen.core.structure.Structure.from_file(config["MatLoc"])
    else:
        warn("MatType not properly specified")
        return -1
    #create the data directory
    if not os.path.isdir("data/" + config["MatID"]):
        os.mkdir("data/" + config["MatID"])

    #ELK functions
    gen_spacegroupin("data/" + config["MatID"], mat)

    return 0




if __name__ == "__main__":
    main()