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

ROUNDING = 4
BtoA = 0.52917721090
mpid = "mp-1079800"
API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
# with MPRester(API_KEY) as mpr:
#     mat = mpr.get_structure_by_material_id(mpid)
# config = Utils.MPSearch.gen_cfg(mpid, mat)
config = "data/TI_EX_Proj/Bi2Se3_orig.json"
config = "data/KLnCuSeProj/K2Dy4Cu6Se9.json"
with open(config) as f:
    config = json.load(f)
main.from_config(config)
# main.download_remote(config["MatLoc"], all=True)
# with MPRester(API_KEY) as mpr: mat = mpr.get_structure_by_material_id(config["MatID"])
# main.run_elk(config)
mat = pymatgen.core.structure.Structure.from_file(config["CIF"])
# print(Utils.MPSearch.get_dftu_information(mat))
# main.gen_slurm(mat,config)
# print(Utils.MPSearch.get_dftu_information(mat))
plot.plot(config, mat, show=True, spins=False, options=["OADOS", "OBS"], energy_range=(-2, 2), dos_range=None,
          el_orbs={"Cu": ["d_{xy}", "d_{xz}", "d_{z^2}","d_{yz}","d_{x^2-y^2}"], "Bi": ["s", "p_x", "p_y", "p_z"],
                   "Se": ["s", "p_x", "p_y", "p_z"], "Ta": ["d", "p_x", "p_y", "p_z"], "O": ["s","p", "p_x", "p_y", "p_z"],
                   "Cl": ["s","p", "p_x", "p_y", "p_z"]}, titles=False, num_points=1000)
#          el_orbs={"Sm": ["f1", "f2", "f3", "f4", "f5", "f6", "f7"], "Fe": ["d_{x^2-y^2}"], "Cu": ["d_{x^2-y^2}"],
#                   "Ba": [], "O": ["p_x", "p_y", "p_z"]}, sites={"Ni": [1, 2, 3, 4]})
#          el_orbs={"Sm":["f1","f2","f3","f4","f5","f6","f7"],"Fe":["d_{xy}", "d_{xz}", "d_{yz}", "d_{z^2}", "d_{x^2-y^2}","p_x","p_y","p_z"]}, sites={"Ni": [1, 2, 3, 4]})
#           el_orbs={"Sm":["f"],"Fe":["d"],"Cu":["d"],"Ba":[],"O":["p","s"]}, sites={"Ni": [1, 2, 3, 4]})
# main.download_remote("data/TIProj/mp-22472/", all=True)
# "d_{xy}", "d_{xz}", "d_{yz}", "d_{z^2}", "d_{x^2-y^2}"
1 / 0

configs = ["K2Sm4Cu2Se7"]  # , "K2Dy4Cu2Se7", "K2Dy4Cu6Se9", "K2Er4Cu6Se9", "K2Gd4Cu6Se9", "K2Ho4Cu6Se9", "K2Lu4Cu6Se9",
# "K2Tb4Cu6Se9", "K2Yb4Cu6Se9"]
for config in configs:
    config = "data/KLnCuSeProj/" + config + ".json"
    main.from_config(config)
1 / 0

plots = True
if plots:
    i = 0
    completed = open("completed_runs.txt", "r")
    with MPRester(API_KEY) as mpr:
        lines = completed.readlines()
        for line in lines:
            i += 1
            if i > 47:
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
