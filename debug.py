import json
from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter
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
import Utils.KPath
import Utils.orbitalAnalysis as orbitalAnalysis
import pandas as pd
import matplotlib.pyplot as plt

ROUNDING = 4
BtoA = 0.52917721090
mpid = "mp-22262"
API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"


def plot2d(WFFile, axislabels):
    frame = pd.read_csv(WFFile, header=None, delim_whitespace=True)
    xvals_orig = np.array(frame.loc[:, 0]) * BtoA
    yvals_orig = np.array(frame.loc[:, 1]) * BtoA

    xoffsets = []
    for i in range(len(yvals_orig)):
        if yvals_orig[i] == np.max(yvals_orig):
            xoffsets.append(xvals_orig[i])
    xoff = np.min(xoffsets)
    # do length
    xvals = np.concatenate((xvals_orig, xvals_orig + np.max(xvals_orig), xvals_orig + 2 * np.max(xvals_orig)))
    yvals = np.concatenate((yvals_orig, yvals_orig, yvals_orig))
    # do height
    xvals = np.concatenate((xvals, xvals - xoff))
    yvals = np.concatenate((yvals, yvals + np.max(yvals)))

    evals = np.array(frame.loc[:, 2]) + 0.00
    evals = np.concatenate((evals, evals, evals, evals, evals, evals))
    plt.xlim([np.min(xvals), np.max(xvals)])
    plt.ylim([np.min(yvals), np.max(yvals)])
    plt.scatter(xvals, yvals, c=np.log(evals), cmap='autumn')
    plt.xlabel(axislabels[0])
    plt.ylabel(axislabels[1])
    plt.colorbar(label=axislabels[2], pad=0, aspect=20, shrink=0.75)
    plt.gca().set_aspect(1.0)
    plt.show()


# with MPRester(API_KEY) as mpr:
#     mat = mpr.get_structure_by_material_id(mpid)
# config = Utils.MPSearch.gen_cfg(mpid, mat)
# config = "data/Other/Ta2Cl4O2_config.json"
# config = "data/mp-1070761_config.json"
# with open(config) as f:
#    config = json.load(f)
# main.run_elk(config)
# main.from_config(config)
# mat = pymatgen.core.structure.Structure.from_file("data/mp-541837/Bi2Se3PrimCell.cif")
# mat = mat.get_primitive_structure()
# print(mat)
# print(Utils.KPath.get_kpoints_SC(mat))
# w = CifWriter(mat, symprec=0.1)
# w.write_file("Ta2Cl2O_geo_opt_pymatgen.cif",)
# print(mat.get_space_group_info())
# plot.plot(config, mat, options=["ODOS", "BS"], energy_range=(-5, 8),          el_orbs={"Ta": ["d", ], "Cl": ["s", "p"], "O": ["s", "p"]})
# main.download_remote(config["MatLoc"], all=True)

# c1 = complex('-1-1j')
# c2 = complex('1+1j')
# print(c1/c2)

#plot2d("data/mp-22924/WFAB123.OUT", ["a axis (" + "$\\AA$)", "b axis (" + "$\\AA$)", r"$log(|\psi_{ik}(r)|^2)$"])

#1/0

config = "data/mp-22924_config.json"
with open(config) as f:
    config = json.load(f)
with MPRester(API_KEY) as mpr: mat = mpr.get_structure_by_material_id(config["MatID"])
#fermi_states = orbitalAnalysis.find_fermi_states(config)
#for state in fermi_states.keys():
#    print(str(state) + ": " + str(fermi_states[state]))
#orbitalAnalysis.find_band_crossings(config, mat)

points = []
lat = mat.lattice.matrix
for site in mat.sites:
    coords = site.coords
    color = 0
    if site.specie.symbol == "Na":
        color = 1

    points.append(list(coords) + [color])
    points.append(list(coords + lat[0]) + [color])
    points.append(list(coords + lat[1]) + [color])
    points.append(list(coords + lat[2]) + [color])
    points.append(list(coords + lat[0] + lat[1]) + [color])
    points.append(list(coords + lat[0] + lat[2]) + [color])
    points.append(list(coords + lat[1] + lat[2]) + [color])
    points.append(list(coords + lat[0] + lat[1] + lat[2]) + [color])
points = np.array(points)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(points[:, 0], points[:, 1], points[:, 2], c=points[:, 3], alpha=1.0, s=100)

frame = pd.read_csv("data/mp-22924/WF3D.OUT", header=None, delim_whitespace=True)
frame = frame.loc[
    ((frame[3] > 0.5 * np.average(np.array(frame.loc[:, 3]))) & (frame[3] < 100 * np.average(np.array(frame.loc[:, 3]))))]
xvals = np.array(frame.loc[:, 0]) * BtoA
yvals = np.array(frame.loc[:, 1]) * BtoA
zvals = np.array(frame.loc[:, 2]) * BtoA
xvals = np.concatenate((xvals, xvals+lat[2][0], xvals+lat[1][0], xvals+lat[0][0], xvals+lat[0][0]+lat[1][0], xvals+lat[0][0]+lat[2][0], xvals+lat[1][0]+lat[2][0]))
yvals = np.concatenate((yvals, yvals+lat[2][1], yvals+lat[1][1], yvals+lat[0][1], yvals+lat[0][1]+lat[1][1], yvals+lat[0][1]+lat[2][1], yvals+lat[1][1]+lat[2][1]))
zvals = np.concatenate((zvals, zvals+lat[2][2], zvals+lat[1][2], zvals+lat[0][2], zvals+lat[0][2]+lat[1][2], zvals+lat[0][2]+lat[2][2], zvals+lat[1][2]+lat[2][2]))

evals = np.array(frame.loc[:, 3])
evals = np.concatenate((evals, evals, evals, evals, evals, evals, evals))

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

print(len(evals))

ax.scatter(xvals, yvals, zvals, lw=0.0, alpha=0.05)
# ax.colorbar()
# ax.gca().set_aspect(1.0)
ax.set_xlim([0, 4])
ax.set_ylim([0, 4])
ax.set_zlim([0, 4])
plt.show()

1 / 0


def plot2d(WFFile, axislabels):
    frame = pd.read_csv(WFFile, header=None, delim_whitespace=True)
    xvals_orig = np.array(frame.loc[:, 0]) * BtoA
    yvals_orig = np.array(frame.loc[:, 1]) * BtoA
    xoff = []
    for i in range(len(yvals_orig)):
        if yvals_orig[i] == np.max(yvals_orig):
            xoff.append(xvals_orig[i])
    xoff = np.min(xoff)
    # xoff = [xvals_orig[i] if yvals_orig[i] == np.min(yvals_orig) else np.max(xvals_orig) for i in np.arange(len(yvals_orig))]
    # do length
    xvals = np.concatenate((xvals_orig, xvals_orig + np.max(xvals_orig), xvals_orig + 2 * np.max(xvals_orig)))
    yvals = np.concatenate((yvals_orig, yvals_orig, yvals_orig))
    # do height
    xvals = np.concatenate((xvals, xvals - xoff))
    yvals = np.concatenate((yvals, yvals + np.max(yvals)))

    evals = np.array(frame.loc[:, 2]) + 0.001
    evals = np.concatenate((evals, evals, evals, evals, evals, evals))

    plt.scatter(xvals, yvals, c=np.log(evals), cmap='autumn')
    plt.xlabel(axislabels[0])
    plt.ylabel(axislabels[1])
    plt.colorbar(label=axislabels[2])
    plt.gca().set_aspect(1.0)
    plt.show()


1 / 0
elements = ["Cl", "O", "Ta"]
col_names = ["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}", "d_{z^2}",
             "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5", "f6", "f7"]
more_species = True
species = 1
char_vals = pd.DataFrame(columns=["atom"] + col_names[2:])
while more_species:
    sites = 1
    more_sites = True
    band_file = "data/Other/Ta2Cl4O2/BAND_S"
    if species < 10:
        band_file += "0" + str(species)
    else:
        band_file += str(species)
    if os.path.exists(band_file + "_A0001.OUT"):
        bands = pd.DataFrame(columns=col_names)
        while more_sites:
            if sites < 10:
                current_band_file = band_file + "_A000" + str(sites) + ".OUT"
            else:
                current_band_file = band_file + "_A00" + str(sites) + ".OUT"
            sites += 1

            if os.path.exists(current_band_file):
                if bands.empty:
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                        names=col_names)
                else:
                    bands = pd.concat([bands, pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                                          names=col_names)], ignore_index=True)
            else:
                more_sites = False
        print(str(elements[species - 1]))
        locations = bands.loc[((bands["vec"] == 0.0) & (bands["energy"] <= 0.009) & (bands["energy"] >= -0.009))]
        new_vals = pd.DataFrame(columns=col_names)
        for column in np.arange(2, 18, 1):
            new_vals.iloc[0, column] = np.sum(locations.iloc[:, column])
        print(new_vals)
        species += 1
    else:
        more_species = False

1 / 0

vhkl = pd.read_csv("data/mp-1070761/mp-1070761_lowq/VHKL.txt", header=None, delim_whitespace=True)
pairs = []
paired_index = []
missing = []
for index, row in vhkl.iterrows():
    vec = list(row)
    if 0 not in vec:
        par = vhkl.index[(vhkl[0] == -vec[0]) & (vhkl[1] == -vec[1]) & (vhkl[2] == -vec[2])].tolist()
        if len(par) > 1:
            print("more than one equal vector: " + str(par))
        elif len(par) == 0:
            missing.append(vec)
        elif index not in paired_index:
            pairs.append([index, par[0]])
            paired_index.append(index)
            paired_index.append(par[0])

df = pd.read_csv("data/mp-1070761/mp-1070761_lowq/GWFS1.txt", header=None)
df = df.iloc[:, :-1]

parity = []
bands = [77, 78, 79, 80]
zero = complex("0+0j")
for band in bands:
    coef_vals = list(df.iloc[band - 1])
    vals = []
    standard = 0
    for pair in pairs:
        values = [complex(coef_vals[pair[0]].replace(" ", "")), complex(coef_vals[pair[1]].replace(" ", ""))]
        if not (values[0] == zero and values[1] == zero):
            if values[0] == zero or values[1] == zero:
                print("err, one val zero")
            else:
                ratio = values[0] / values[1]
                vals.append(ratio)
    if len(vals) == 0:
        parity.append(0)
    else:
        parity.append(np.average(vals))
print(parity)

1 / 0
parity = []
bands = [77, 78, 79, 80]
for band in bands:
    coef_vals = list(df.iloc[band - 1])
    vals = []
    standard = 0
    for pair in pairs:
        if (coef_vals[pair[0]][0] == "-" and coef_vals[pair[1]][0] != "-") or (
                coef_vals[pair[0]][0] != "-" and coef_vals[pair[1]][0] == "-"):
            vals.append(-1)
            if standard != 0 and standard != -1:
                print(coef_vals[pair[0]] + "    " + coef_vals[pair[1]] + "   " + str(
                    list(vhkl.iloc[pair[0]])) + "    " + str(list(vhkl.iloc[pair[1]])))
        else:
            vals.append(1)
            if standard != 0 and standard != 1:
                print(coef_vals[pair[0]] + "    " + coef_vals[pair[1]] + "   " + str(
                    list(vhkl.iloc[pair[0]])) + "    " + str(list(vhkl.iloc[pair[1]])))
        if len(vals) == 1:
            standard = vals[0]
    if len(vals) == 0:
        parity.append(0)
    else:
        parity.append(np.average(vals))
print(parity)

1 / 0

parity = []
bands = [62]
zero = complex("0+0j")
for band in bands:
    coef_vals = list(df.iloc[band - 1])
    vals = []
    standard = 0
    for pair in pairs:
        values = [complex(coef_vals[pair[0]]), complex(coef_vals[pair[1]])]
        if not (values[0] == zero and values[1] == zero):
            print(values)
            if values[0] == zero or values[1] == zero:
                print("err, one val zero")
            else:
                ratio = values[0] / values[1]
                vals.append(ratio)
    if len(vals) == 0:
        parity.append(0)
    else:
        parity.append(np.average(vals))
print(parity)

1 / 0

df = pd.read_csv("data/mp-1070761/GWFS1.txt", header=None)
df = df.iloc[:, :-1]
parity = []
bands = [77, 78, 79, 80]
for band in bands:
    coef_vals = list(df.iloc[band - 1])
    vals = []
    standard = 0
    for pair in pairs:
        if (coef_vals[pair[0]][0] == "-" and coef_vals[pair[1]][0] != "-") or (
                coef_vals[pair[0]][0] != "-" and coef_vals[pair[1]][0] == "-"):
            vals.append(-1)
            if standard != 0 and standard != -1:
                print(coef_vals[pair[0]] + "    " + coef_vals[pair[1]] + "   " + str(
                    list(vhkl.iloc[pair[0]])) + "    " + str(list(vhkl.iloc[pair[1]])))
        else:
            vals.append(1)
            if standard != 0 and standard != 1:
                print(coef_vals[pair[0]] + "    " + coef_vals[pair[1]] + "   " + str(
                    list(vhkl.iloc[pair[0]])) + "    " + str(list(vhkl.iloc[pair[1]])))
        if len(vals) == 1:
            standard = vals[0]
    if len(vals) == 0:
        parity.append(0)
    else:
        parity.append(np.average(vals))
print(parity)

1 / 0

df = pd.read_csv("data/mp-1070761/GWFS2.txt", header=None)
df = df.iloc[:, :-1]
parity = []
bands = [77, 78, 79, 80]
for band in bands:
    coef_vals = list(df.iloc[band - 1])
    vals = []
    for i in np.arange(1, len(coef_vals) - 200, 2):
        print(coef_vals[i] + "    " + coef_vals[i + 1])
        if (coef_vals[i][0] == "-" and coef_vals[i + 1][0] != "-") or (
                coef_vals[i][0] != "-" and coef_vals[i + 1][0] == "-"):
            vals.append(-1)
        else:
            vals.append(1)
    if len(vals) == 0:
        parity.append(0)
    else:
        parity.append(np.average(vals))
print(parity)

1 / 0

parity = []
for k in range(df.shape[0]):
    tmp = list(df.iloc[k])
    coef_vals = [complex(num[:-1] + "j") for num in tmp]
    vals = []
    for i in np.arange(1, len(coef_vals) - 1, 2):
        if not (coef_vals[i] == complex("0+0j") and coef_vals[i + 1] == complex("0+0j")):
            print(k)
            if coef_vals[i] == complex("0+0j") or coef_vals[i + 1] == complex("0+0j"):
                print("warn, coefficient equal to zero, other not")
            else:
                ratio = coef_vals[i] / coef_vals[i + 1]
                if not ratio.imag == 0:
                    print(ratio)
                vals.append(ratio.real)
    if len(vals) == 0:
        parity.append(0)
    else:
        parity.append(np.average(vals))
print(parity)
