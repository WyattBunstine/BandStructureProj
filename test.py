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
mpid = "mp-4826"
API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
# with MPRester(API_KEY) as mpr:
#     mat = mpr.get_structure_by_material_id(mpid)
# config = Utils.MPSearch.gen_cfg(mpid, mat)
# config = "data/Other/Ta2Cl4O2_config.json"
# config = "data/Other/Ta2Cl4O2_centered_config.json"
# config = "data/mp-4826_config.json"
# config = "data/TIProj/mp-1214324/mp-1214324_config.json"
# config = "data/mp-20311_config.json"
config = "data/Ag2+Proj/Ag_Orig.json"
with open(config) as f:
    config = json.load(f)
main.from_config(config)
#with MPRester(API_KEY) as mpr: mat = mpr.get_structure_by_material_id(config["MatID"])
# main.run_elk(config)
# mat = pymatgen.core.structure.Structure.from_file(config["CIF"])
#plot.plot(config, mat, show=True, spins=True, options=["OADOS", "IDOS", "BS"], energy_range=(-3, 8),
#          el_orbs={"Na": [],"Bi": ["p_x","p_y","p_z"]}, sites={"Ni": [1, 2, 3, 4]})
#          el_orbs={"Sm": ["f1", "f2", "f3", "f4", "f5", "f6", "f7"], "Fe": ["d_{x^2-y^2}"], "Cu": ["d_{x^2-y^2}"],
#                   "Ba": [], "O": ["p_x", "p_y", "p_z"]}, sites={"Ni": [1, 2, 3, 4]})
#          el_orbs={"Sm":["f1","f2","f3","f4","f5","f6","f7"],"Fe":["d_{xy}", "d_{xz}", "d_{yz}", "d_{z^2}", "d_{x^2-y^2}","p_x","p_y","p_z"]}, sites={"Ni": [1, 2, 3, 4]})
#           el_orbs={"Sm":["f"],"Fe":["d"],"Cu":["d"],"Ba":[],"O":["p","s"]}, sites={"Ni": [1, 2, 3, 4]})
#main.download_remote(config["MatLoc"], all=True)
# "d_{xy}", "d_{xz}", "d_{yz}", "d_{z^2}", "d_{x^2-y^2}"

1/0

plots = True
if plots:
    i = 0
    completed = open("completed_runs.txt", "r")
    with MPRester(API_KEY) as mpr:
        lines = completed.readlines()
        for line in lines:
            i += 1
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
