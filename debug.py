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
import math
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

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
          '#17becf']


# Function to find distance
def shortest_distance(point, plane):
    (a, b, c) = np.cross(plane[0] - plane[1], plane[0] - plane[2])
    d = abs((a * point[0] + b * point[1] + c * point[2]))
    e = (math.sqrt(a * a + b * b + c * c))
    return d / e


def plot2d(WFFile, axislabels, mat: pymatgen.core.structure.Structure, plane):
    atoms = {}
    lat = mat.lattice.matrix
    vert = [(lat[0] - lat[1]) / (np.sqrt(np.sum((lat[0] + lat[1]) ** 2))),
            (lat[0] + lat[1]) / (np.sqrt(np.sum((lat[0] + lat[1]) ** 2))),
            np.cross(lat[0], lat[1]) / np.sqrt(np.sum(np.cross(lat[0], lat[1]) ** 2))]
    for site in mat.sites:
        coords = site.coords
        color = 0
        if site.specie.symbol not in atoms.keys():
            atoms[site.specie.symbol] = []
        if shortest_distance(mat.lattice.get_fractional_coords(coords), plane) < 0.1:
            for i in np.arange(4):
                coords = site.coords
                if i == 0:
                    coords = np.dot(vert, coords)
                if i == 1:
                    coords = np.dot(vert, coords + lat[0] + lat[1])
                if i == 2:
                    coords = np.dot(vert, coords + lat[2])
                if i == 3:
                    coords = np.dot(vert, coords + lat[0] + lat[1] + lat[2])

                atoms[site.specie.symbol].append(list(coords) + [color])

    frame = pd.read_csv(WFFile, header=None, delim_whitespace=True, skiprows=[0])
    xvals_orig = np.array(frame.loc[:, 0]) * BtoA
    yvals_orig = np.array(frame.loc[:, 1]) * BtoA

    reflecty = False
    if reflecty:
        xvals_orig = np.concatenate((xvals_orig, xvals_orig))
        yvals_orig = np.concatenate((yvals_orig, -yvals_orig))
        yvals_orig = np.add(yvals_orig, np.max(yvals_orig))

    xoffsets = []
    x_shift = []
    for i in range(len(yvals_orig)):
        if yvals_orig[i] == np.max(yvals_orig):
            xoffsets.append(xvals_orig[i])
        if yvals_orig[i] == np.min(yvals_orig):
            x_shift.append(xvals_orig[i])
    xoff = np.min(xoffsets)
    x_shift = np.max(x_shift) - np.min(x_shift) + x_shift[1] - x_shift[0]

    y_shift = np.max(yvals_orig)
    for i in range(len(yvals_orig)):
        if yvals_orig[i] < y_shift and yvals_orig[i] != np.min(yvals_orig):
            y_shift = yvals_orig[i]
    y_shift = np.max(yvals_orig) + y_shift

    orig_maxes = [np.max(xvals_orig), np.max(yvals_orig)]

    # do length
    xvals = np.concatenate((xvals_orig, xvals_orig + x_shift, xvals_orig + 2 * x_shift, xvals_orig + 3 * x_shift))
    yvals = np.concatenate((yvals_orig, yvals_orig, yvals_orig, yvals_orig))
    # do height
    xvals = np.concatenate((xvals, xvals + xoff, xvals + 2 * xoff))
    yvals = np.concatenate((yvals, yvals + y_shift, yvals + 2 * y_shift))

    evals = np.round(np.array(frame.loc[:, 2]), 10) + 0.0005

    if reflecty:
        evals = np.concatenate((evals, evals))

    xvals = xvals - orig_maxes[0]
    #yvals = yvals - orig_maxes[1]

    evals = np.concatenate((evals, evals, evals, evals))
    evals = np.concatenate((evals, evals, evals))
    plt.xlim([0, orig_maxes[0] * 3])
    plt.ylim([0, y_shift * 3])
    # plt.scatter(xvals, yvals, c=np.log(evals), cmap='autumn', s=0.1)
    plt.scatter(xvals, yvals, c=np.log(evals), cmap='inferno', s=0.3, vmax=None)

    plt.xlabel(axislabels[0])
    plt.ylabel(axislabels[1])
    plt.colorbar(label=axislabels[2], pad=0, aspect=20, shrink=0.75, ticks=[])

    cell_points1 = np.array([[0.0, 0.0], [3.19217, np.sqrt(5.94219 ** 2 + 3.82398 ** 2)],
                             [6.38434, np.sqrt(5.94219 ** 2 + 3.82398 ** 2)], [3.19217, 0.00000], [0.0, 0.0]])

    cell_points = np.array(
        [[0, 0], [4.78826 - 1.59609, 0], [4.78826 - 1.59609, np.sqrt((4.20151 - 1.74067) ** 2 + 3.82398 ** 2)],
         [0, np.sqrt((4.20151 - 1.74067) ** 2 + 3.82398 ** 2)], [0, 0]])

    #plt.scatter(3 * (4.78826 - 1.59609) / 2, 3 * np.sqrt((4.20151 - 1.74067) ** 2 + 3.82398 ** 2) / 2, s=20,
    #            color="white")
    #plt.scatter(3.19217, np.sqrt(5.94219 ** 2 + 3.82398 ** 2) / 2, s=20, color="blue")
    #plt.plot(np.transpose(cell_points)[0]+(4.78826 - 1.59609),
    #         np.transpose(cell_points)[1]+np.sqrt((4.20151 - 1.74067) ** 2 + 3.82398 ** 2), color="white", linewidth=2)

    i = 2
    for atom in atoms.keys():
        points = np.array(atoms[atom])
        if len(points) < 0:
            plt.scatter(points[:, 1], points[:, 2], color=colors[i], label=atom, alpha=1.0, s=20)
        i += 1
    plt.gca().set_aspect(1.0)
    plt.show()


# config = "data/TaOCl2/Ta2Cl4O2_centered_config.json"
# config = "data/TaOCl2/Ta2Cl4O2_config.json"
config = "data/OrbitalPathway/CaCuO2.json"

with open(config) as f:
    config = json.load(f)
mat = pymatgen.core.structure.Structure.from_file(config["CIF"])
centered_plane = np.array([[0.0, 0.5, 0.0], [1.0, 0.5, 0.0], [0.0, 0.0, 1.0]])
dimer_place = plane = np.array([[0.25, 0.0, 0.0], [0.25, 1.0, 0.0], [0.25, 0.0, 1.0]])

plot2d(config["MatLoc"] + "WF2D.OUT",
       ["a axis (" + "$\\AA$)", "b axis (" + "$\\AA$)", r"$\propto|\psi_{ik}(r)|^2$"], mat, dimer_place)

1/0
mat = mat.to_primitive()
points = []
species = []
i = 0
for site in mat:
    points.append(site.coords)
#print(points)

points = np.array(points)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(points[0], points[1], points[2], s=50)

frame = pd.read_csv(config["MatLoc"] + "WF3D.OUT", header=None, delim_whitespace=True, dtype="float", skiprows=[0])
frames = []
frames.append(frame.loc[
    ((frame[3] > 0.85 * np.average(np.array(frame.loc[:, 3]))) & (
            frame[3] < 0.9 * np.average(np.array(frame.loc[:, 3]))))])

frames.append(frame.loc[
    ((frame[3] > 1.5 * np.average(np.array(frame.loc[:, 3]))) & (
            frame[3] < 1.55 * np.average(np.array(frame.loc[:, 3]))))])

frames.append(frame.loc[
    ((frame[3] > 3 * np.average(np.array(frame.loc[:, 3]))) & (
            frame[3] < 3.05 * np.average(np.array(frame.loc[:, 3]))))])

for i, frame in enumerate(frames):
    xvals = np.array(frame.loc[:, 0]) * BtoA
    yvals = np.array(frame.loc[:, 1]) * BtoA
    zvals = np.array(frame.loc[:, 2]) * BtoA
    evals = np.array(frame.loc[:, 3])

    #surf = ax.plot_trisurf(xvals, yvals, zvals, linewidth=0)

    ax.scatter(xvals, yvals, zvals, lw=0.0, alpha=0.05, c=colors[i])
# ax.colorbar()
# ax.gca().set_aspect(1.0)
ax.set_xlim([0, 6])
ax.set_ylim([0, 6])
ax.set_zlim([0, 6])
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

vhkl = pd.read_csv("data/NiW2B2/mp-1070761_lowq/VHKL.txt", header=None, delim_whitespace=True)
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

df = pd.read_csv("data/NiW2B2/mp-1070761_lowq/GWFS1.txt", header=None)
df = df.iloc[:, :-1]
print(df.shape)

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
        parity.append(np.round(np.average(vals), 3))
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

df = pd.read_csv("data/NiW2B2/GWFS1.txt", header=None)
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

df = pd.read_csv("data/NiW2B2/GWFS2.txt", header=None)
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
