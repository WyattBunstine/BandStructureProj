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
import matplotlib.pyplot as plt

ROUNDING = 4
BtoA = 0.52917721090
mpid = "mp-1079800"
API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
#with MPRester(API_KEY) as mpr:
#     mat = mpr.get_structure_by_material_id(mpid)
# config = Utils.MPSearch.gen_cfg(mpid, mat)
# config = "data/NKCTO/NKCTO.json"
config = ("data/TIProj_candidates/Si.json")
#config = "data/TIProj_candidates/YCr4Cu3O12.json"
# config = "data/Li2CuO2/Li2CuO2_phonon_calc_4.json"

config = "data/OrbitalPathway/CaCuO2.json"
# config = "data/NiW2B2/NiW2B2_bilayer.json"

with open(config) as f:
    config = json.load(f)
# main.from_config(config)
# main.download_remote(config["MatLoc"], all=True)
# with MPRester(API_KEY) as mpr: mat = mpr.get_structure_by_material_id(config["MatID"])
# main.run_elk(config)
mat = pymatgen.core.structure.Structure.from_file(config["CIF"])
# print(mat.lattice.reciprocal_lattice)
# main.run_phonon(mat, config, restart=True)
# print(Utils.MPSearch.get_dftu_information(mat))
# main.gen_slurm(mat,config)
# print(Utils.MPSearch.get_dftu_information(mat))

plot.plot(config, mat, show=True, spins=True, options=["OADOS", "BS"], energy_range=(-4, 6), dos_range=None,
          el_orbs={
              "Cu": ["d_{xy}", "d_{xz}", "d_{yz}", "d_{z^2}", "d_{x^2-y^2}"],  #, "p_x", "p_y", "p_z"],
              "O": [ "p_x", "p_y", "p_z"],
              "Na": ["s", "p_x", "p_y", "p_z"],
              "Bi": ["s", "p_x", "p_y", "p_z"]},
          titles=True, num_points=500)

1/0

f, ax = plt.subplots(1, 3, sharey=True, figsize=(20, 10), gridspec_kw={'width_ratios': [1, 1, 1]})


hb1_yvals = np.arange(-4, -1.9, 0.1)
xvals = -((-3 - hb1_yvals) ** 2 - 1) * 2

hb1_yvals = np.arange(-4, -2.95, 0.05)
hb2_yvals = np.arange(-0.8, 0.25, 0.05) - 0.3
hb3_yvals = np.arange(-0.2, 0.85, 0.05) + 0.3
hb4_yvals = np.arange(3, 4.05, 0.05)

c1 = "blue"
c2 = "red"

for axis in [ax[0], ax[2]]:

    axis.fill_betweenx(hb3_yvals, xvals, 0, color=c2, alpha=0.1)
    axis.plot(xvals, hb3_yvals, c=c2)

    axis.fill_betweenx(hb4_yvals, xvals, 0, color=c2, alpha=0.1)
    axis.plot(xvals, hb4_yvals, c=c2)

for axis in [ax[1], ax[2]]:
    axis.fill_betweenx(hb1_yvals, xvals, 0, color=c1, alpha=0.1)
    axis.plot(xvals, hb1_yvals, c=c1)


    axis.plot(xvals, hb2_yvals, c=c1)

ax[1].fill_betweenx(hb2_yvals, xvals, 0, color=c1, alpha=0.1)

ax[2].fill_betweenx(hb2_yvals[:-5], xvals[:-5], 0, color=c1, alpha=0.1)
ax[2].fill_betweenx(hb2_yvals[-5:], xvals[-5:], 0, color=c2, alpha=1)

ax[2].fill_betweenx(hb3_yvals[:5], xvals[:5], 0, color=c1, alpha=1)
ax[2].fill_betweenx(hb3_yvals[5:], xvals[5:], 0, color=c2, alpha=0.1)


for axis in ax:
    axis.plot([-10, 10], [0, 0], linestyle="dashed", c="black")
    axis.plot([0, 0], [-10, 10], linestyle="solid", c="black")
    axis.set_xlim((0, 3))
    axis.set_ylim((-4, 4))
    axis.set_xticks([])
    axis.set_yticks([])

ax[0].set_title("Mott-Hubbard System 1")
ax[1].set_title("Mott-Hubbard System 2")
ax[2].set_title("Combined Mott-Hubbard System")



plt.show()

1 / 0
Utils.Phonopy_conv.main()

1 / 0

for file1 in os.listdir(config["MatLoc"]):
    if "DYN" in file1:
        for file2 in os.listdir(config["MatLoc"]):
            if "DYN" in file2 and file1 != file2:
                with open(config["MatLoc"] + file1) as f1:
                    with open(config["MatLoc"] + file2) as f2:
                        same = True
                        for line in f1.readlines():
                            if f2.readline() != line:
                                same = False
                        if same:
                            print(file1 + "   " + file2)

1 / 0
with open("rm.sh", "w+") as sh:
    for file in os.listdir(config["MatLoc"]):
        if "DYN" in file:
            write = False
            with open(config["MatLoc"] + file) as f:
                if len(f.readlines()) < 3:
                    write = True
                    sh.write("rm " + file + "\n")
            if write and False:
                with open(config["MatLoc"] + file, "w+") as f:
                    if len(f.readlines()) < 3:
                        write = True
                        f.write('''  0 0 : is =    1, ia =    1, ip =    1
 0 0 : is =    1, ia =    1, ip =    2
 0 0 : is =    1, ia =    1, ip =    3
 0 0 : is =    2, ia =    1, ip =    1
 0 0 : is =    2, ia =    1, ip =    2
 0 0 : is =    2, ia =    1, ip =    3
 0 0 : is =    2, ia =    2, ip =    1
 0 0 : is =    2, ia =    2, ip =    2
 0 0 : is =    2, ia =    2, ip =    3
 0 0 : is =    3, ia =    1, ip =    1
 0 0 : is =    3, ia =    1, ip =    2
 0 0 : is =    3, ia =    1, ip =    3
 0 0 : is =    3, ia =    2, ip =    1
 0 0 : is =    3, ia =    2, ip =    2
 0 0 : is =    3, ia =    2, ip =    3
''')
1 / 0

with open("qpath.txt", "w+") as f:
    r1 = 25
    r2 = 25
    r3 = 25
    for i in range(r1):
        f.write("  " + str(float(i) / (2 * float(r1))) + " " + str(-float(i) / (2 * float(r1))) + " " + str(
            float(i) / (2 * float(r1))) + "\n")
    for i in range(r2):
        f.write(
            "  0.5 " + str(-1 / 2 + float(i) / (2 * float(r2))) + " " + str(1 / 2 - float(i) / (2 * float(r2))) + "\n")
    for i in range(r3 + 1):
        f.write("  " + str(1 / 2 - float(i) / (2 * float(r3))) + " " + str(float(i) / (2 * float(r3))) + " 0.0\n")

1 / 0

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


def get_phonon_info(config):
    with open("INFOS.txt", "w+") as sh:
        for file in os.listdir(config["MatLoc"]):
            if "INFO" in file:
                with open(config["MatLoc"] + file, "rb") as f:
                    try:  # catch OSError in case convergence info can't be found
                        f.seek(-2, os.SEEK_END)
                        seek = True
                        while seek:
                            text = f.read(23)
                            # either find the end of looping or find the convergence info
                            if text == b"RMS change in Kohn-Sham":
                                seek = False
                            else:
                                f.seek(-24, os.SEEK_CUR)  # the block to catch an error finding convergence info
                    except OSError:
                        print("Failed finding INFO.OUT convergence information.")
                        f.seek(0)
                    potential = (str(text + f.readline()) + "\n").split()
                    energy = str(f.readline()).split()
                    if float(potential[7]) > 0.00001:
                        print(file + "\nPotential: " + potential[7] + "\nEnergy:    " + energy[7] + "\n")
                        sh.write(file + "\n")
                        sh.write("Potential: " + potential[7] + "\n")
                        sh.write("Energy:    " + energy[7] + "\n")
