import json
from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen
from pymatgen.core import *
from pymatgen.symmetry.kpath import *
import numpy as np
import subprocess
import os
import itertools
import plot
import main
import Utils.MPSearch
import yaml

def main():

    config = "/home/wyatt/PycharmProjects/BandStructureProj/data/Li2CuO2/Li2CuO2_phonon_calc_2.json"
    with open(config) as f:
        config = json.load(f)
    mat = pymatgen.core.structure.Structure.from_file(config["CIF"])

    prim = mat.to_primitive()

    band = {}
    band["calculator"] = "elk"
    #number of special kpoints on path
    band["npath"] = 3
    #number of qpoints between each special point
    band["segment_nqpoint"] = [26,26,26]
    #reciprocal lattice
    band["reciprocal_lattice"] = [str(el) for el in prim.lattice.reciprocal_lattice.matrix.tolist()]
    #number of atoms in unit cell
    band["natom"] = len(prim.sites)
    #lattice
    band["lattice"] = [str(el) for el in prim.lattice.matrix.tolist()]
    #crystal sites
    points = []
    for site in prim.sites:
        site_dict = {}
        site_dict["symbol"] = site.specie.symbol
        site_dict["coordinates"] = str(prim.lattice.get_fractional_coords(site.coords).tolist())
        site_dict["mass"] = site.specie.atomic_mass.real
        points.append(site_dict)
    band["points"] = points

    phonons = []
    phonon_num = 0
    num_atom = 0
    bands = []
    phonon = {}
    eigenvecs = {}
    eigenvectors = []
    atom_eigenvector = []
    labels = {1: "*$\\Gamma$*", 26: "*$X$*", 51: "*$S$*", 76: "*$R$*"}
    with open("/home/wyatt/PycharmProjects/BandStructureProj/data/Li2CuO2/Li2CuO2_phonon_2/PHONON.OUT") as phonon_file:
        for line in phonon_file:
            line = line.split()
            if "q-point," in line:
                phonon_num += 1
                if len(phonon.keys()) > 0:
                    eigenvecs["eigenvector"] = eigenvectors
                    bands.append(eigenvecs)
                    phonon["band"] = bands
                    phonons.append(phonon)
                    if "label" in phonon.keys() and (phonon["label"] == "*$X$*" or phonon["label"] == "*$S$*"):
                        phonons.append(phonon)
                phonon = {}
                bands = []
                eigenvecs = {}
                eigenvectors = []
                atom_eigenvector = []
                num_atom = 0
                phonon["q-position"] = "[" + str(float(line[1])) + "," + str(float(line[2])) + "," + str(
                    float(line[3])) + "]"
                phonon["distance"] = float(line[0]) / 100
                if int(line[0]) in labels.keys():
                    phonon["label"] = labels[int(line[0])]
                atom_counter = 0
            elif "frequency" in line:
                if len(eigenvecs.keys()) > 0:
                    eigenvecs["eigenvector"] = eigenvectors
                    bands.append(eigenvecs)
                    eigenvecs = {}
                #need to convert the frequencies from hartree to Thz, conversion is 1h = 6579.6897 mEv
                eigenvecs["frequency"] = float(line[1]) * 6579.6897
                eigenvectors = []
                atom_eigenvector = []
                num_atom = 0
            elif len(line) > 3:

                eig1 = round(float(line[3]), 5)
                eig2 = round(float(line[4]), 5)
                atom_eigenvector.append("[" + str(eig1) + "," + str(eig2) + "]")
                if len(atom_eigenvector) > 2:
                    eigenvectors.append(atom_eigenvector)
                    num_atom+=1
                    atom_eigenvector = []


    eigenvecs["eigenvector"] = eigenvectors
    bands.append(eigenvecs)
    phonon["band"] = bands
    phonons.append(phonon)
    band["nqpoint"] = len(phonons)

    band["phonon"] = phonons



    with open("tmp_band.yaml", "w+") as stream:
        yaml.Dumper.ignore_aliases = lambda *args: True
        yaml.dump(band, stream, sort_keys=False)
    with open("tmp_band.yaml") as tmp:
        with open("Li2CuO2_phonon_2.yaml", "w+") as stream:
            for line in tmp.readlines():
                line = line.replace('\'', "")
                line = line.replace('*', "\'")
                stream.write(line)

if __name__ == "__main__":
    main()