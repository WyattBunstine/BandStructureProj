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
import Utils.Phonopy_conv
import yaml

ROUNDING = 4
BtoA = 0.52917721090
mpid = "mp-1079800"
API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
# with MPRester(API_KEY) as mpr:
#     mat = mpr.get_structure_by_material_id(mpid)
# config = Utils.MPSearch.gen_cfg(mpid, mat)
# config = "data/Other/CoNiO3.json"
# config = "data/KLnCuSeProj/K2Dy4Cu6Se9.json"
# config = "data/TIProj_candidates/CaLi3MnCrO6.json"
config = "data/Li2CuO2/Li2CuO2_phonon_calc_3.json"

with open(config) as f:
    config = json.load(f)
main.from_config(config)
# main.download_remote(config["MatLoc"], all=True)
# with MPRester(API_KEY) as mpr: mat = mpr.get_structure_by_material_id(config["MatID"])
# main.run_elk(config)
# mat = pymatgen.core.structure.Structure.from_file(config["CIF"])
# print(Utils.MPSearch.get_dftu_information(mat))
# main.gen_slurm(mat,config)
#print(Utils.MPSearch.get_dftu_information(mat))
#plot.plot(config, mat, show=True, spins=False, options=["EDOS", "BS"], energy_range=(-2, 2), dos_range=None,
#          el_orbs={"Co": ["d_{xy}", "d_{xz}", "d_{z^2}","d_{yz}","d_{x^2-y^2}"],
#                   "Se": ["p_x", "p_y", "p_z"], "O": ["p", "p_x", "p_y", "p_z"],}, titles=False, num_points=1000)

1/0

with open("rm.sh","w+") as sh:
    for file in os.listdir(config["MatLoc"]):
        if "DYN" in file:
            with open(config["MatLoc"] + file) as f:
                if len(f.readlines())<3:
                    sh.write("rm " + file + "\n")

1/0

with open("qpath.txt","w+") as f:
    r1 = 25
    r2 = 25
    r3 = 25
    for i in range(r1):
        f.write("  "+ str(float(i)/(2*float(r1))) + " " + str(-float(i)/(2*float(r1))) + " "+str(float(i)/(2*float(r1)))+"\n")
    for i in range(r2):
        f.write("  0.5 " + str(-1/2 + float(i)/(2*float(r2))) + " " + str(1/2 - float(i)/(2*float(r2))) + "\n")
    for i in range(r3+1):
        f.write("  " + str(1/2 - float(i)/(2*float(r3))) + " " + str(float(i)/(2*float(r3))) + " 0.0\n")



1/0

Utils.Phonopy_conv.main()

1 / 0
config_files = ["data/KLnCuSeProj/K2Ho4Cu6Se9_n=1.json",
                "data/KLnCuSeProj/K2Ho4Cu6Se9.json", "data/KLnCuSeProj/K2Ho4Cu6Se9_n=3.json"]
#configs = []
#mats = []
#for file in config_files:
#    with open(file) as f:
#        configs.append(json.load(f))
#for config in configs:
#    mats.append(pymatgen.core.structure.Structure.from_file(config["CIF"]))

#plot.tri_plot(configs,mats, show=True, spins=False, options=["EDOS", "BS"], energy_range=(-2, 2), dos_range=None)

#          el_orbs={"Co": ["d_{xy}", "d_{xz}", "d_{z^2}","d_{yz}","d_{x^2-y^2}"],
#                   "Se": ["p_x", "p_y", "p_z"], "O": ["p", "p_x", "p_y", "p_z"],}, titles=False, num_points=1000)
#          el_orbs={"Sm": ["f1", "f2", "f3", "f4", "f5", "f6", "f7"], "Fe": ["d_{x^2-y^2}"], "Cu": ["d_{x^2-y^2}"],
#                   "Ba": [], "O": ["p_x", "p_y", "p_z"]}, sites={"Ni": [1, 2, 3, 4]})
#          el_orbs={"Sm":["f1","f2","f3","f4","f5","f6","f7"],"Fe":["d_{xy}", "d_{xz}", "d_{yz}", "d_{z^2}",
#          "d_{x^2-y^2}","p_x","p_y","p_z"]}, sites={"Ni": [1, 2, 3, 4]})
#           el_orbs={"Sm":["f"],"Fe":["d"],"Cu":["d"],"Ba":[],"O":["p","s"]}, sites={"Ni": [1, 2, 3, 4]})
# main.download_remote("data/TIProj/mp-22472/", all=True)
# "d_{xy}", "d_{xz}", "d_{yz}", "d_{z^2}", "d_{x^2-y^2}"

#configs = ["K2Sm4Cu2Se7"]  # , "K2Dy4Cu2Se7", "K2Dy4Cu6Se9", "K2Er4Cu6Se9", "K2Gd4Cu6Se9", "K2Ho4Cu6Se9", "K2Lu4Cu6Se9",
# "K2Tb4Cu6Se9", "K2Yb4Cu6Se9"]
#for config in configs:
#    config = "data/KLnCuSeProj/" + config + ".json"
#    main.from_config(config)


plots = True
if plots:
    i = 0
    completed = open("completed_runs.txt", "r")
    with MPRester(API_KEY) as mpr:
        lines = completed.readlines()
        for line in lines:
            i += 1
            if i > 0:
                with open("data/TIProj/" + line[:-1] + "/" + line[:-1] + "_config.json") as f:
                    config = json.load(f)
                    mat = mpr.get_structure_by_material_id(line[:-1])
                    plot.plot(config, mat, options=["EDOS", "BS"], spins=True)
            """try:
                with open("data/TIProj/" + line[:-1] + "/" + line[:-1] + "_config.json") as f:
                    config = json.load(f)
                    mat = mpr.get_structure_by_material_id(line[:-1])
                    plot.plot(config, mat, options=["EDOS", "BS"], spins=True)
            except Exception:
                print(line[:-1] + " could not be plotted.")"""
    completed.close()
