import itertools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import matplotlib.axes
import pandas as pd
import Utils.KPath
import pymatgen
import os
from matplotlib.colors import LinearSegmentedColormap
from pymatgen.symmetry.kpath import *
from matplotlib.collections import LineCollection

orbitals = ["s", "p_y", "p_z", "p_x", "d_{xy}", "d_{xz}", "d_{z^2}", "d_{yz}", "d_{x^2-y^2}", "f1", "f2", "f3",
            "f4", "f5", "f6", "f7"]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
          '#17becf']
line_styles = ["solid", "dotted", "dashed", "dashdot", ":"]
markers = [".", "x", ">", "<"]


def plot_spin(config, mat: pymatgen.core.structure.Structure, axis):
    mappable = None
    # plot bandlines
    more_species = True
    species = 1
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "BAND_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        if os.path.exists(band_file + "_A0001.OUT"):
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1
                if os.path.exists(current_band_file):
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True, names=["vec", "energy",
                                                                                                      "spin_up",
                                                                                                      "spin_down"])
                    # bands = bands.loc[(bands["energy"] > -5 / 27.2138) & (bands["energy"] < 8 / 27.2138)]
                    xvals = np.array(list(bands['vec']))
                    yvals = np.array(list(bands['energy'])) * 27.2138
                    char = (np.array(list(bands['spin_down'])) * 2 /
                            (np.array(list(bands['spin_up'])) + np.array(list(bands['spin_down'])))) - 1

                    mappable = axis.scatter(xvals, yvals, c=char, s=1, label='_nolegned_', cmap="autumn")
                else:
                    more_sites = False
            species += 1
        else:
            more_species = False

    # plot significant points
    sigpoints = pd.read_csv(config["MatLoc"] + "BANDLINES.OUT", header=None, delim_whitespace=True, names=["x", "y"])
    xvals = np.array(list(sigpoints['x']))
    yvals = np.array(list(sigpoints['y'])) * 27.2138
    length = 2
    ticks = []
    for i in range(0, len(xvals), length):
        axis.plot(xvals[i:i + length], yvals[i:i + length], label='_nolegned_', color="black", linestyle='dashed',
                  linewidth=0.5)
        ticks.append(xvals[i])

    sigpointlabels = []
    if config["kpoints"] == "MP":
        kpoints = Utils.KPath.get_kpoints_SC(mat)
    elif config["kpoints"] == "McQueen":
        kpoints = Utils.KPath.get_kpoints_McQueen(mat)
    else:
        warn("kpoints not properly specified")
        return -1
    for point in kpoints:
        if point[0] == "\\Gamma":
            sigpointlabels.append("\u0393")
        else:
            sigpointlabels.append("$" + str(point[0]) + "$")

    axis.set_xticks(ticks, sigpointlabels)
    axis.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black", linewidth=0.5)  # ,label="Fermi Level")
    axis.set_xlabel("Wave Vector")
    axis.set_title(str(mat.composition.to_latex_string()) + " Band Structure")
    return mappable


def plot_2_orbital_bands(config, mat: pymatgen.core.structure.Structure, axis, fig, el_orbs=None, energy_range=(-8, 5)):
    # plot bandlines

    used_colors = 0
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    if el_orbs is None or len(el_orbs) > 2:
        print("too many orbitals for 2 orbital bands, should be 2")

    orb1 = np.zeros(1)
    orb2 = np.zeros(1)

    orb1_name = ""
    orb2_name = ""

    xvals = np.zeros(1)
    yvals = np.zeros(1)

    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "BAND_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        if os.path.exists(band_file + "_A0001.OUT") and atoms[species - 1] in el_orbs.keys():
            bands = pd.DataFrame(columns=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                          "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5", "f6", "f7"])
            band_vals = None
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1

                if os.path.exists(current_band_file):
                    if bands.empty:
                        bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                            names=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                                   "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5",
                                                   "f6", "f7"])
                        #bands = bands.loc[(bands["energy"] > energy_range[0] / 27.2138) & (
                        #        bands["energy"] < energy_range[1] / 27.2138)]
                        bands = bands.reset_index()
                        band_vals = bands.to_numpy()[:, 3:]
                    else:
                        bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                            names=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                                   "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5",
                                                   "f6", "f7"])
                        #bands = bands.loc[(bands["energy"] > energy_range[0] / 27.2138) & (
                        #        bands["energy"] < energy_range[1] / 27.2138)]
                        bands = bands.reset_index()
                        band_vals += bands.to_numpy()[:, 3:]


                else:
                    more_sites = False

            xvals = np.array(list(bands['vec']))
            yvals = np.array(list(bands['energy'])) * 27.2138

            if "s" in el_orbs[atoms[species - 1]]:
                if len(orb1) == 1:
                    orb1 = bands[orbitals[0]]
                    orb1_name = atoms[species - 1] + " $s$"
                elif len(orb2) == 1:
                    orb2 = bands[orbitals[0]]
                    orb2_name = atoms[species - 1] + " $s$"
            if "p" in el_orbs[atoms[species - 1]]:
                if len(orb1) == 1:
                    orb1 = bands[orbitals[1]] + bands[orbitals[2]] + bands[orbitals[3]]
                    orb1_name = atoms[species - 1] + " $p$"
                elif len(orb2) == 1:
                    orb2 = bands[orbitals[1]] + bands[orbitals[2]] + bands[orbitals[3]]
                    orb2_name = atoms[species - 1] + " $p$"
            if "d" in el_orbs[atoms[species - 1]]:
                if len(orb1) == 1:
                    orb1 = bands[orbitals[4]] + bands[orbitals[5]] + bands[orbitals[6]] + bands[orbitals[7]] + bands[
                        orbitals[8]]
                    orb1_name = atoms[species - 1] + " $d$"
                elif len(orb2) == 1:
                    orb2 = bands[orbitals[4]] + bands[orbitals[5]] + bands[orbitals[6]] + bands[orbitals[7]] + bands[
                        orbitals[8]]
                    orb2_name = atoms[species - 1] + " $d$"

            species += 1
        elif not os.path.exists(band_file + "_A0001.OUT"):
            more_species = False
        else:
            species += 1

    points = np.array([xvals, yvals]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    segments = np.delete(segments, np.arange(999, segments.shape[0], config["numBandPoints"]), axis=0)

    colors = np.array((orb2 - orb1) / (orb2 + orb1))
    colors = np.delete(colors, np.arange(999, segments.shape[0], config["numBandPoints"]), axis=0)

    inset_ax = axis.inset_axes(
        [0.2, 0.4, 0.005, 0.005],  # [x, y, width, height] w.r.t. axes
        xlim=[1.1, 1.20], ylim=[-0.05, 0.05],  # sets viewport & tells relation to main axes
        xticklabels=[], yticklabels=[]
    )
    inset_ax.xaxis.set_ticks_position("none")
    inset_ax.yaxis.set_ticks_position("none")

    lc = LineCollection(segments, linewidths=5, cmap="seismic")
    lc.set_array(colors)
    inset_ax.add_collection(lc)

    lc = LineCollection(segments, linewidths=2, cmap="seismic")
    lc.set_array(colors)
    line = axis.add_collection(lc)

    axis.indicate_inset_zoom(inset_ax, edgecolor="blue")

    fig.colorbar(line, ax=axis, label="(" + orb2_name + " - " + orb1_name + ")/(" + orb1_name + " + " + orb2_name + ")",
                 location="right")

    # plot significant points
    sigpoints = pd.read_csv(config["MatLoc"] + "BANDLINES.OUT", header=None, delim_whitespace=True, names=["x", "y"])
    xvals = np.array(list(sigpoints['x']))
    yvals = np.array(list(sigpoints['y'])) * 27.2138
    length = 2
    ticks = []
    for i in range(0, len(xvals), length):
        axis.plot(xvals[i:i + length], yvals[i:i + length], label='_nolegned_', color="black", linestyle='dashed',
                  linewidth=0.5)
        ticks.append(xvals[i])

    sigpointlabels = []
    if config["kpoints"] == "MP":
        kpoints = Utils.KPath.get_kpoints_SC(mat)
    elif config["kpoints"] == "McQueen":
        kpoints = Utils.KPath.get_kpoints_McQueen(mat)
    else:
        kpoints = config["kpoints"]
    for point in kpoints:
        if point[0] == "\\Gamma":
            sigpointlabels.append("\u0393")
        else:
            sigpointlabels.append("$" + str(point[0]) + "$")
    axis.set_xticks(ticks, sigpointlabels)
    axis.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black", linewidth=0.5)  # ,label="Fermi Level")
    axis.set_xlabel("Wave Vector")


def plot_orbital_bands(config, mat: pymatgen.core.structure.Structure, axis, el_orbs=None, energy_range=(-8, 5)):
    # plot bandlines

    used_colors = 0
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    band_frames = pd.DataFrame(columns=["vec", "energy", "el", "orb", "val"])
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "BAND_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        if os.path.exists(band_file + "_A0001.OUT"):# and atoms[species - 1] in el_orbs.keys():
            bands = pd.DataFrame(columns=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                          "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5", "f6", "f7"])
            band_vals = None
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1

                if os.path.exists(current_band_file):
                    if bands.empty:
                        bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                            names=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                                   "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5",
                                                   "f6", "f7"])
                        bands = bands.loc[(bands["energy"] > energy_range[0] / 27.2138) & (
                                bands["energy"] < energy_range[1] / 27.2138)]
                        bands = bands.reset_index()
                        band_vals = bands.to_numpy()[:, 3:]
                    else:
                        bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                            names=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                                   "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5",
                                                   "f6", "f7"])
                        bands = bands.loc[(bands["energy"] > energy_range[0] / 27.2138) & (
                                bands["energy"] < energy_range[1] / 27.2138)]
                        bands = bands.reset_index()
                        band_vals += bands.to_numpy()[:, 3:]

                else:
                    more_sites = False

            band_vals[:, 1] = band_vals[:, 1] + band_vals[:, 2] + band_vals[:, 3]
            band_vals[:, 2] = band_vals[:, 4] + band_vals[:, 5] + band_vals[:, 6] + band_vals[:,
                                                                                    7] + band_vals[:, 8]
            band_vals[:, 3] = (band_vals[:, 9] + band_vals[:, 10] + band_vals[:, 11] + band_vals[:, 12]
                               + band_vals[:, 13] + band_vals[:, 14] + band_vals[:, 15])
            if species == 1:
                for index, row in bands.iterrows():
                    band_frames.loc[len(band_frames.index)] = ([row["vec"], row["energy"]] +
                                                               [atoms[species - 1], np.argmax(band_vals[index]),
                                                                np.max(band_vals[index])])
                    if np.argmax(band_vals[index]) > 3:
                        print("bad max val")
                        print(band_vals[index])
            else:
                for index, row in bands.iterrows():
                    if band_frames.loc[index]["val"] < np.max(band_vals[index]):
                        band_frames.loc[index] = (band_frames.loc[index]["vec"], band_frames.loc[index]["energy"],
                                                  atoms[species - 1], np.argmax(band_vals[index]),
                                                  np.max(band_vals[index]))
                        if np.argmax(band_vals[index]) > 3:
                            print("bad max val")
                            print(band_vals[index])

            species += 1
        elif not os.path.exists(band_file + "_A0001.OUT"):
            more_species = False
        else:
            species += 1
    for atom in atoms:
        if atom in el_orbs.keys():
            for index, orb in enumerate(["s", "p", "d", "f"]):
                if orb in el_orbs[atom]:
                    orb_bands = band_frames.loc[(band_frames["orb"] == index)]
                    orb_bands = orb_bands.loc[(orb_bands["el"] == atom)]
                    xvals = np.array(list(orb_bands['vec']))
                    yvals = np.array(list(orb_bands['energy'])) * 27.2138
                    axis.scatter(xvals, yvals, s=20, label=atom + "$" + orb + "$",
                                 c=colors[used_colors % 10], marker=markers[int(used_colors / 10)])
                    used_colors += 1

    # plot significant points
    sigpoints = pd.read_csv(config["MatLoc"] + "BANDLINES.OUT", header=None, delim_whitespace=True, names=["x", "y"])
    xvals = np.array(list(sigpoints['x']))
    yvals = np.array(list(sigpoints['y'])) * 27.2138
    length = 2
    ticks = []
    for i in range(0, len(xvals), length):
        axis.plot(xvals[i:i + length], yvals[i:i + length], label='_nolegned_', color="black", linestyle='dashed',
                  linewidth=0.5)
        ticks.append(xvals[i])

    sigpointlabels = []
    if config["kpoints"] == "MP":
        kpoints = Utils.KPath.get_kpoints_SC(mat)
    elif config["kpoints"] == "McQueen":
        kpoints = Utils.KPath.get_kpoints_McQueen(mat)
    else:
        kpoints = config["kpoints"]
    for point in kpoints:
        if point[0] == "\\Gamma":
            sigpointlabels.append("\u0393")
        else:
            sigpointlabels.append("$" + str(point[0]) + "$")
    axis.set_xticks(ticks, sigpointlabels)
    axis.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black", linewidth=0.5)  # ,label="Fermi Level")
    axis.set_xlabel("Wave Vector")
    #lgnd = axis.legend(loc="lower left", scatterpoints=1, fontsize=10)
    #for i in lgnd.legendHandles:
    #    i._sizes = [30]


def plot_ang_orbital_bands(config, mat: pymatgen.core.structure.Structure, axis, el_orbs=None, energy_range=(-8, 5)):
    # plot bandlines

    used_colors = 0
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    band_frames = pd.DataFrame(columns=["vec", "energy", "el", "orb", "val"])
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "BAND_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        if os.path.exists(band_file + "_A0001.OUT") and atoms[species - 1] in el_orbs.keys():
            bands = pd.DataFrame(columns=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                          "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5", "f6", "f7"])
            band_vals = None
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1

                if os.path.exists(current_band_file):
                    if bands.empty:
                        bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                            names=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                                   "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5",
                                                   "f6", "f7"])
                        bands = bands.loc[(bands["energy"] > energy_range[0] / 27.2138) & (
                                bands["energy"] < energy_range[1] / 27.2138)]
                        bands = bands.reset_index()
                        band_vals = bands.to_numpy()[:, 3:]
                    else:
                        bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                            names=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                                   "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5",
                                                   "f6", "f7"])
                        bands = bands.loc[(bands["energy"] > energy_range[0] / 27.2138) & (
                                bands["energy"] < energy_range[1] / 27.2138)]
                        bands = bands.reset_index()
                        band_vals += bands.to_numpy()[:, 3:]
                else:
                    more_sites = False
            if species == 1:
                for index, row in bands.iterrows():
                    band_frames.loc[len(band_frames.index)] = ([row["vec"], row["energy"]] +
                                                               [atoms[species - 1], np.argmax(band_vals[index]),
                                                                np.max(band_vals[index])])
            else:
                for index, row in bands.iterrows():
                    if band_frames.loc[index]["val"] < np.max(band_vals[index]):
                        band_frames.loc[index] = (band_frames.loc[index]["vec"], band_frames.loc[index]["energy"],
                                                  atoms[species - 1], np.argmax(band_vals[index]),
                                                  np.max(band_vals[index]))
            species += 1
        elif not os.path.exists(band_file + "_A0001.OUT"):
            more_species = False
        else:
            species += 1
    for atom in atoms:
        if atom in el_orbs.keys():
            for i in np.arange(16):
                if orbitals[i] in el_orbs[atom]:
                    orb_bands = band_frames.loc[(band_frames["orb"] == i)]
                    orb_bands = orb_bands.loc[(orb_bands["el"] == atom)]
                    xvals = np.array(list(orb_bands['vec']))
                    yvals = np.array(list(orb_bands['energy'])) * 27.2138
                    axis.scatter(xvals, yvals, s=20, label=atom + "$" + orbitals[i] + "$",
                                 c=colors[used_colors % 10], marker=markers[int(used_colors / 10)])
                    used_colors += 1

    # plot significant points
    sigpoints = pd.read_csv(config["MatLoc"] + "BANDLINES.OUT", header=None, delim_whitespace=True, names=["x", "y"])
    xvals = np.array(list(sigpoints['x']))
    yvals = np.array(list(sigpoints['y'])) * 27.2138
    length = 2
    ticks = []
    for i in range(0, len(xvals), length):
        axis.plot(xvals[i:i + length], yvals[i:i + length], label='_nolegned_', color="black", linestyle='dashed',
                  linewidth=0.5)
        ticks.append(xvals[i])

    sigpointlabels = []
    if config["kpoints"] == "MP":
        kpoints = Utils.KPath.get_kpoints_SC(mat)
    elif config["kpoints"] == "McQueen":
        kpoints = Utils.KPath.get_kpoints_McQueen(mat)
    else:
        kpoints = config["kpoints"]
    for point in kpoints:
        if point[0] == "\\Gamma":
            sigpointlabels.append("\u0393")
        else:
            sigpointlabels.append("$" + str(point[0]) + "$")
    axis.set_xticks(ticks, sigpointlabels)
    axis.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black", linewidth=0.5)  # ,label="Fermi Level")
    axis.set_xlabel("Wave Vector")
    #lgnd = axis.legend(loc="lower left", scatterpoints=1, fontsize=10)
    #for i in lgnd.legendHandles:
    #    i._sizes = [30]


def plot_elemental_bands(config, mat: pymatgen.core.structure.Structure, show=True):
    plt.figure(figsize=(10, 10))
    more_species = True
    species = 1
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "BAND_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        if os.path.exists(band_file + "_A0001.OUT"):
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1
                if os.path.exists(current_band_file):
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True, names=["vec", "energy",
                                                                                                      "spin_up",
                                                                                                      "spin_down"])
                    xvals = np.array(list(bands['vec']))
                    yvals = np.array(list(bands['energy'])) * 27.2138
                    plt.scatter(xvals, yvals, label='_nolegned_', color=colors[species], alpha=0.5, s=1)
                else:
                    more_sites = False
            species += 1
        else:
            more_species = False

    # plot significant points
    sigpoints = pd.read_csv(config["MatLoc"] + "BANDLINES.OUT", header=None, delim_whitespace=True, names=["x", "y"])
    xvals = np.array(list(sigpoints['x']))
    yvals = np.array(list(sigpoints['y'])) * 27.2138
    length = 2
    ticks = []
    for i in range(0, len(xvals), length):
        plt.plot(xvals[i:i + length], yvals[i:i + length], label='_nolegned_', color="black", linestyle='dashed',
                 linewidth=0.5)
        ticks.append(xvals[i])

    sigpointlabels = []
    if config["kpoints"] == "MP":
        kpoints = Utils.KPath.get_kpoints_SC(mat)
    elif config["kpoints"] == "McQueen":
        kpoints = Utils.KPath.get_kpoints_McQueen(mat)
    else:
        warn("kpoints not properly specified")
        return -1
    for point in kpoints:
        if point[0] == "\\Gamma":
            sigpointlabels.append("\u0393")
        else:
            sigpointlabels.append("$" + str(point[0]) + "$")

    plt.xticks(ticks, sigpointlabels)
    plt.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black", linewidth=0.5)  # ,label="Fermi Level")
    # plt.ylim(-5, 8)
    plt.legend()
    plt.ylabel("E-Ef (eV)")
    plt.xlabel("Wave Vector")
    plt.title(str(mat.composition.to_latex_string()) + " Band Structure")
    if show:
        plt.show()
    else:
        plt.savefig(config["MatLoc"] + str(mat.composition.to_pretty_string()) + "EB.png", dpi=250)
        plt.savefig("data/plots/" + str(mat.composition.to_pretty_string()) + "EB.png", dpi=250)


def plot_orb_DOS(config, mat: pymatgen.core.structure.Structure, axis: matplotlib.axes.Axes, both_spins=False,
                 el_orbs=None, energy_range=(-5, 8), dos_range=None, num_points=500):
    orbitals = ["s", "p", "d", "f"]
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    x_range = [0, 0]
    used_colors = 0
    range = [-1, -1]
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "PDOS_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        species_colors = 0
        if os.path.exists(band_file + "_A0001.OUT"):
            xvals = np.zeros(num_points * 32)
            yvals = np.zeros(num_points)
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1
                if os.path.exists(current_band_file):
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True, names=["energy", "dos"])
                    xvals += np.array(list(bands["dos"]))
                    yvals = np.array(list(bands['energy']))[0:num_points] * 27.2138

                else:
                    more_sites = False

            for yval in np.arange(len(yvals)):
                if yvals[yval] > energy_range[0] and range[0] == -1:
                    range[0] = yval
                if yvals[yval] > energy_range[1] and range[1] == -1:
                    range[1] = yval

            if both_spins:
                m = 0
                dos = np.zeros(num_points)
                for i in np.arange(0, num_points * 16, num_points):
                    dos += xvals[i:i + num_points]
                    m += 1
                    if (i < num_points and m == 1) or (i < num_points * 4 and m == 3) or (
                            i < num_points * 9 and m == 5) or (m == 7):
                        orbital = orbitals[int(m / 2)]
                        element = str(atoms[species - 1])
                        if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                            style = int(used_colors / len(colors))
                            color = used_colors % len(colors)
                            axis.plot(dos, yvals, label=element + " $" + orbital + '$', color=colors[color],
                                      linestyle=line_styles[style])
                            used_colors += 1
                            species_colors += 1
                        if np.max(dos[range[0]:range[1]]) > x_range[1]:
                            x_range[1] = np.max(dos[range[0]:range[1]])
                        m = 0
                        dos = np.zeros(num_points)
                if both_spins:
                    used_colors -= species_colors
                    m = 0
                    dos = np.zeros(num_points)
                    for i in np.arange(num_points * 16, num_points * 32, num_points):
                        dos += xvals[i:i + num_points]
                        m += 1
                        if (i < num_points * 17 and m == 1) or (i < num_points * 20 and m == 3) or (
                                i < num_points * 25 and m == 5) or (
                                m == 7):
                            orbital = orbitals[int(m / 2.0)]
                            element = str(atoms[species - 1])
                            if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                                style = int(used_colors / len(colors))
                                color = used_colors % len(colors)
                                axis.plot(dos, yvals, label="_nolabel_", color=colors[color],
                                          linestyle=line_styles[style])
                                used_colors += 1
                            if np.min(dos[range[0]:range[1]]) < x_range[0]:
                                x_range[0] = np.min(dos[range[0]:range[1]])
                            m = 0
                            dos = np.zeros(num_points)
            else:
                xvals = np.abs(xvals)
                element = str(atoms[species - 1])
                s_dos = np.zeros(num_points)
                p_dos = np.zeros(num_points)
                d_dos = np.zeros(num_points)
                f_dos = np.zeros(num_points)

                for i in np.arange(0, num_points * 32, num_points):
                    if i < num_points or (num_points * 15 < i < num_points * 17):
                        s_dos += xvals[i:i + num_points]
                    elif i < 4 * num_points or (num_points * 16 < i < num_points * 20):
                        p_dos += xvals[i:i + num_points]
                    elif i < 9 * num_points or (num_points * 19 < i < num_points * 25):
                        d_dos += xvals[i:i + num_points]
                    else:
                        f_dos += xvals[i:i + num_points]
                for index, dos in enumerate([s_dos, p_dos, d_dos, f_dos]):
                    if el_orbs is None or (element in el_orbs.keys() and orbitals[index] in el_orbs[element]):
                        orbital = orbitals[index]
                        style = int(used_colors / len(colors))
                        color = used_colors % len(colors)
                        axis.plot(dos, yvals, label=element + " $" + orbital + '$', color=colors[color],
                                  linestyle=line_styles[style])
                        used_colors += 1
                        if np.min(dos[range[0]:range[1]]) < x_range[0]:
                            x_range[0] = np.min(dos[range[0]:range[1]])
                        if np.max(dos[range[0]:range[1]]) > x_range[1]:
                            x_range[1] = np.max(dos[range[0]:range[1]])

            species += 1
        else:
            more_species = False

    axis.plot(x_range, [0, 0], label="_nolabel_", linestyle='dashed', color="black",
              linewidth=0.5)  # ,label="Fermi Level")
    axis.plot([0, 0], [-100, 100], label="_nolabel_", color="black", linewidth=0.5)
    axis.legend()
    if dos_range is not None:
        axis.set_xlim(dos_range)
    else:
        axis.set_xlim(x_range)
    axis.set_xlabel("Density of States")
    axis.get_xaxis().set_ticks([])


def plot_orb_ang_DOS(config, mat: pymatgen.core.structure.Structure, axis: matplotlib.axes.Axes, both_spins=False,
                     el_orbs=None, energy_range=(-5, 8), num_points=500):
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    x_range = [0, 0]
    used_styles = 0
    used_colors = 0
    range = [-1, -1]
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "PDOS_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        species_colors = 0
        if os.path.exists(band_file + "_A0001.OUT"):
            element = str(atoms[species - 1])
            if (el_orbs is None or (element in el_orbs.keys())):
                xvals = np.zeros(16 * num_points)
                yvals = np.zeros(num_points)
                while more_sites:
                    if sites < 10:
                        current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                    else:
                        current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                    sites += 1
                    if os.path.exists(current_band_file):
                        bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True,
                                            names=["energy", "dos"])
                        if not len(list(bands["dos"])) == len(xvals):
                            xvals = np.zeros(len(list(bands["dos"])))
                        xvals += np.array(list(bands["dos"]))
                        yvals = np.array(list(bands['energy']))[0:num_points] * 27.2138

                    else:
                        more_sites = False
                m = 0
                for i in np.arange(0, 16 * num_points, num_points):
                    dos = xvals[i:i + num_points]
                    orbital = orbitals[m]
                    m += 1
                    element = str(atoms[species - 1])
                    if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                        style = int(used_colors / len(colors))
                        color = used_colors % len(colors)
                        axis.plot(dos, yvals, label=element + " $" + orbital + '$', color=colors[color],
                                  linestyle=line_styles[style])
                        used_colors += 1
                        species_colors += 1
                    for yval in np.arange(len(yvals)):
                        if yvals[yval] > energy_range[0] and range[0] == -1:
                            range[0] = yval
                        if yvals[yval] > energy_range[1] and range[1] == -1:
                            range[1] = yval
                    if np.max(dos[range[0]:range[1]]) > x_range[1]:
                        x_range[1] = np.max(dos[range[0]:range[1]])

                if both_spins:
                    used_colors -= species_colors
                    m = 0
                    for i in np.arange(16 * num_points, 32 * num_points, num_points):
                        dos = xvals[i:i + num_points]
                        orbital = orbitals[m]
                        m += 1
                        element = str(atoms[species - 1])
                        if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                            style = int(used_colors / len(colors))
                            color = used_colors % len(colors)
                            axis.plot(dos, yvals, label="_nolabel_", color=colors[color], linestyle=line_styles[style])
                            used_colors += 1
                        if np.min(dos[range[0]:range[1]]) < x_range[0]:
                            x_range[0] = np.min(dos[range[0]:range[1]])
            species += 1
        else:
            more_species = False

    axis.plot(x_range, [0, 0], label="_nolabel_", linestyle='dashed', color="black",
              linewidth=0.5)  # ,label="Fermi Level")
    axis.plot([0, 0], [-100, 100], label="_nolabel_", color="black", linewidth=0.5)
    axis.legend()
    axis.set_xlim(x_range)
    axis.set_xlabel("Density of States")
    axis.set_title(str(mat.composition.to_latex_string()) + " DOS")
    axis.get_xaxis().set_ticks([])


def plot_orb_ang_site_DOS(config, mat: pymatgen.core.structure.Structure, axis: matplotlib.axes.Axes, both_spins=False,
                          el_orbs=None, ele_sites=None):
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    x_range = [0, 0]
    used_colors = 0
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "PDOS_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)

        if os.path.exists(band_file + "_A0001.OUT"):
            xvals = np.zeros(16000)
            yvals = np.zeros(500)
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1
                species_colors = 0
                if os.path.exists(current_band_file):
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True, names=["energy", "dos"])
                    xvals = np.array(list(bands["dos"]))
                    yvals = np.array(list(bands['energy']))[0:500] * 27.2138
                    m = 0
                    for i in np.arange(0, 8000, 500):
                        dos = xvals[i:i + 500]
                        orbital = orbitals[m]
                        m += 1
                        element = str(atoms[species - 1])
                        if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                            if ele_sites is None or (element in ele_sites.keys() and sites - 1 in ele_sites[element]):
                                axis.plot(dos, yvals, label=element + " " + str(sites - 1) + " $" + orbital + '$',
                                          color=colors[used_colors])
                                used_colors += 1
                                species_colors += 1
                                if np.max(dos) > x_range[1]:
                                    x_range[1] = np.max(dos)
                    if both_spins:
                        used_colors -= species_colors
                        m = 0
                        for i in np.arange(8000, 16000, 500):
                            dos = xvals[i:i + 500]
                            orbital = orbitals[m]
                            m += 1
                            element = str(atoms[species - 1])
                            if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                                if ele_sites is None or (sites - 1 in ele_sites[element]):
                                    axis.plot(dos, yvals, label="_nolabel_", color=colors[used_colors])
                                    used_colors += 1
                                    if np.min(dos) < x_range[0]:
                                        x_range[0] = np.min(dos)
                else:
                    more_sites = False

            species += 1
        else:
            more_species = False

    axis.plot(x_range, [0, 0], label="_nolabel_", linestyle='dashed', color="black",
              linewidth=0.5)  # ,label="Fermi Level")
    axis.plot([0, 0], [-100, 100], label="_nolabel_", color="black", linewidth=0.5)
    axis.legend()
    axis.set_xlabel("Density of States")
    axis.set_title(str(mat.composition.to_latex_string()) + " DOS")
    axis.get_xaxis().set_ticks([])


def plot_orb_site_DOS(config, mat: pymatgen.core.structure.Structure, axis: matplotlib.axes.Axes, both_spins=False,
                      el_orbs=None):
    orbitals = ["s", "p", "d", "f"]
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    x_range = [0, 0]
    used_colors = 0
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "PDOS_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)

        if os.path.exists(band_file + "_A0001.OUT"):
            xvals = np.zeros(16000)
            yvals = np.zeros(500)
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1
                species_colors = 0
                if os.path.exists(current_band_file):
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True, names=["energy", "dos"])
                    xvals = np.array(list(bands["dos"]))
                    yvals = np.array(list(bands['energy']))[0:500] * 27.2138
                    m = 0
                    dos = np.zeros(500)
                    for i in np.arange(0, 8000, 500):
                        dos += xvals[i:i + 500]
                        m += 1
                        if (i < 500 and m == 1) or (i < 2000 and m == 3) or (i < 4500 and m == 5) or m == 7:
                            orbital = orbitals[int(m / 2)]
                            element = str(atoms[species - 1])
                            if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                                print(element + " " + str(sites - 1) + " $" + orbital + '$' + str(used_colors))
                                axis.plot(dos, yvals, label=element + " " + str(sites - 1) + " $" + orbital + '$',
                                          color=colors[used_colors])
                                used_colors += 1
                                species_colors += 1
                            if np.max(dos) > x_range[1]:
                                x_range[1] = np.max(dos)
                            m = 0
                            dos = np.zeros(500)
                    if both_spins:
                        used_colors -= species_colors
                        m = 0
                        dos = np.zeros(500)
                        for i in np.arange(8000, 16000, 500):
                            dos += xvals[i:i + 500]
                            m += 1
                            if (i < 8500 and m == 1) or (i < 10000 and m == 3) or (i < 12500 and m == 5) or m == 7:
                                orbital = orbitals[int(m / 2.0)]
                                element = str(atoms[species - 1])
                                if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                                    axis.plot(dos, yvals, label="_nolabel_", color=colors[used_colors])
                                    used_colors += 1
                                if np.min(dos) < x_range[0]:
                                    x_range[0] = np.min(dos)
                                m = 0
                                dos = np.zeros(500)
                else:
                    more_sites = False

            species += 1
        else:
            more_species = False

    axis.plot(x_range, [0, 0], label="_nolabel_", linestyle='dashed', color="black",
              linewidth=0.5)  # ,label="Fermi Level")
    axis.plot([0, 0], [-100, 100], label="_nolabel_", color="black", linewidth=0.5)
    axis.legend()
    axis.set_xlabel("Density of States")
    axis.set_title(str(mat.composition.to_latex_string()) + " DOS")
    axis.get_xaxis().set_ticks([])

def plot_site_DOS(config, mat: pymatgen.core.structure.Structure, axis: matplotlib.axes.Axes, both_spins=False,
                      el_orbs=None, num_points=500):
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    x_range = [0, 0]
    used_colors = 0
    line_style = 0
    while more_species:
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "PDOS_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)

        if os.path.exists(band_file + "_A0001.OUT"):
            xvals = np.zeros(16000)
            yvals = np.zeros(500)
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1
                species_colors = 0
                if os.path.exists(current_band_file):
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True, names=["energy", "dos"])
                    xvals = np.array(list(bands["dos"]))
                    yvals = np.array(list(bands['energy']))[0:num_points] * 27.2138
                    dos = np.zeros(num_points)
                    for i in np.arange(0, num_points*16, num_points):
                        dos += xvals[i:i + num_points]
                    element = str(atoms[species - 1])
                    if el_orbs is None or (element in el_orbs.keys()):
                        axis.plot(dos, yvals, label=element + " " + str(sites - 1),
                                  color=colors[used_colors], linestyle=line_styles[line_style])
                        used_colors += 1
                        if used_colors >= len(colors):
                            line_style += 1
                            used_colors = 0
                        species_colors += 1
                    if np.max(dos) > x_range[1]:
                        x_range[1] = np.max(dos)

                    if both_spins:
                        used_colors -= species_colors
                        dos = np.zeros(num_points)
                        for i in np.arange(num_points*16, num_points*32, num_points):
                            dos += xvals[i:i + num_points]

                        element = str(atoms[species - 1])
                        if el_orbs is None or (element in el_orbs.keys()):
                            axis.plot(dos, yvals, label="_nolabel_", color=colors[used_colors])
                            used_colors += 1
                        if np.min(dos) < x_range[0]:
                            x_range[0] = np.min(dos)

                else:
                    more_sites = False

            species += 1
        else:
            more_species = False

    axis.plot(x_range, [0, 0], label="_nolabel_", linestyle='dashed', color="black",
              linewidth=0.5)  # ,label="Fermi Level")
    axis.plot([0, 0], [-100, 100], label="_nolabel_", color="black", linewidth=0.5)
    axis.legend()
    axis.set_xlabel("Density of States")
    axis.set_title(str(mat.composition.to_latex_string()) + " DOS")
    axis.get_xaxis().set_ticks([])



def plot_ele_DOS(config, mat: pymatgen.core.structure.Structure, axis, both_spins=False, energy_range=(-5, 8),
                 dos_range=None, num_points=500):
    atoms = sorted([el.symbol for el in mat.composition.elements])
    more_species = True
    species = 1
    used_colors = 0
    x_range = [0, 0]
    while more_species:
        species_colors = 0
        sites = 1
        more_sites = True
        band_file = config["MatLoc"] + "PDOS_S"
        if species < 10:
            band_file += "0" + str(species)
        else:
            band_file += str(species)
        if os.path.exists(band_file + "_A0001.OUT"):
            xvals = np.zeros(16 * num_points)
            yvals = np.zeros(num_points)
            while more_sites:
                if sites < 10:
                    current_band_file = band_file + "_A000" + str(sites) + ".OUT"
                else:
                    current_band_file = band_file + "_A00" + str(sites) + ".OUT"
                sites += 1
                if os.path.exists(current_band_file):
                    bands = pd.read_csv(current_band_file, header=None, delim_whitespace=True, names=["energy", "dos"])
                    if not len(list(bands["dos"])) == len(xvals):
                        xvals = np.zeros(len(list(bands["dos"])))
                    xvals += np.array(list(bands["dos"]))
                    yvals = np.array(list(bands['energy']))[0:num_points] * 27.2138
                else:
                    more_sites = False
            dos = np.zeros(num_points)
            for i in np.arange(0, num_points * 16, num_points):
                dos += xvals[i:i + num_points]
            element = str(atoms[species - 1])
            style = int(used_colors / len(colors))
            color = used_colors % len(colors)
            axis.plot(dos, yvals, label=element, color=colors[color],
                      linestyle=line_styles[style])
            used_colors += 1
            species_colors += 1
            range = [-1, -1]
            for yval in np.arange(len(yvals)):
                if yvals[yval] > energy_range[0] and range[0] == -1:
                    range[0] = yval
                if yvals[yval] > energy_range[1] and range[1] == -1:
                    range[1] = yval
            if np.max(dos[range[0]:range[1]]) > x_range[1]:
                x_range[1] = np.max(dos[range[0]:range[1]])

            dos = np.zeros(num_points)
            if both_spins:
                used_colors -= species_colors
                dos = np.zeros(num_points)
                for i in np.arange(num_points * 16, num_points * 32, num_points):
                    dos += xvals[i:i + num_points]

                style = int(used_colors / len(colors))
                color = used_colors % len(colors)
                axis.plot(dos, yvals, label="_nolabel_", color=colors[color], linestyle=line_styles[style])
                used_colors += 1
                if np.min(dos[range[0]:range[1]]) < x_range[0]:
                    x_range[0] = np.min(dos[range[0]:range[1]])

            species += 1
        else:
            more_species = False

    axis.plot(x_range, [0, 0], linestyle='dashed', color="black")  # ,label="Fermi Level")
    axis.plot([0, 0], [-100, 100], label="_nolabel_", color="black", linewidth=0.5)


    #path = [[-115, -115, 155, 155, -115], [2.75, 3.75, 3.75, 2.75, 2.75]]
    #axis.plot(path[0], path[1], color="black",linewidth=3)

    axis.legend()
    if dos_range is not None:
        axis.set_xlim(dos_range)
    else:
        axis.set_xlim(x_range)
    axis.set_xlabel("Density of States")
    axis.get_xaxis().set_ticks([])
    axis.set_title(str(mat.composition.to_latex_string()) + " DOS")


def plot_TDOS(config, mat: pymatgen.core.structure.Structure, axis, both_spins=False):
    band_file = config["MatLoc"] + "TDOS.OUT"
    if os.path.exists(band_file):
        bands = pd.read_csv(band_file, header=None, delim_whitespace=True, names=["energy", "dos"])
        xvals = np.array(list(bands["dos"]))
        yvals = np.array(list(bands['energy']))[0:500] * 27.2138
        if both_spins:
            axis.plot(xvals[0:500], yvals, color="black")  # label=r'$\uparrow$')
            axis.plot(xvals[500:1000], yvals, color="black")  # label=r'$\downarrow$')
            axis.plot([np.min(xvals), np.max(xvals)], [0, 0], linestyle='dashed', color="black", label="_nolabel_")
            axis.legend()
        else:
            axis.plot(xvals[0:500], yvals, color="black", label="TDOS")
            axis.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black", linewidth=0.5, label="_nolabel_")
        axis.plot([0, 0], [np.min(yvals), np.max(yvals)], color="black", linewidth=0.5, label="_nolabel_")
    else:
        print("Cannot find TDOS file")
    axis.set_xlabel("DOS")
    axis.get_xaxis().set_ticks([])
    axis.set_title(str(mat.composition.to_latex_string()) + " DOS")


def plot_IDOS(config, mat, axis, both_spins=False, num_points=500):
    band_file = config["MatLoc"] + "IDOS.OUT"
    if os.path.exists(band_file):
        bands = pd.read_csv(band_file, header=None, delim_whitespace=True, names=["energy", "dos"])
        xvals = np.array(list(bands["dos"]))
        yvals = np.array(list(bands['energy']))[0:num_points] * 27.2138
        if both_spins:
            axis.plot(xvals[0:num_points], yvals, color="black", label="IDOS")
            axis.plot(xvals[num_points:num_points * 2], yvals, color="black")  # label=r'$\downarrow$')
        else:
            axis.plot(xvals[0:num_points], yvals, color="black", label="IDOS")
    else:
        print("Cannot find IDOS file")


def plot_bands(config, mat: pymatgen.core.structure.Structure, axis):
    # plot bandlines
    band_file = config["MatLoc"] + "BAND.OUT"
    if "OBS" in config["tasks"]:
        band_file = config["MatLoc"] + "BAND_S01_A0001.OUT"

    if os.path.exists(band_file):
        if "OBS" in config["tasks"]:
            print(["vec", "energy"] + list(range(16)))
            bands = pd.read_csv(band_file, header=None, delim_whitespace=True,
                                names=["vec", "energy"] + list(range(16)))
        else:
            bands = pd.read_csv(band_file, header=None, delim_whitespace=True, names=["vec", "energy"])
        xvals = np.array(list(bands['vec']))[0:config["numBandPoints"]]
        yvals = np.array(list(bands['energy'])) * 27.2138
        length = config["numBandPoints"]
        for i in range(length, len(yvals), length):
            axis.plot(xvals, yvals[i:i + length], linewidth=1.0, color="black", label='_nolegned_')
    elif os.path.exists(config["MatLoc"] + "BAND_S01_A0001.OUT"):
        band_file = config["MatLoc"] + "BAND_S01_A0001.OUT"
        bands = pd.read_csv(band_file, header=None, delim_whitespace=True, names=["vec", "energy", "s1", "s2"])
        xvals = np.array(list(bands['vec']))[0:config["numBandPoints"]]
        yvals = np.array(list(bands['energy'])) * 27.2138
        length = config["numBandPoints"]
        for i in range(length, len(yvals), length):
            axis.plot(xvals, yvals[i:i + length], linewidth=1.0, color="black", label='_nolegned_')
    # plot significant points
    sigpoints = pd.read_csv(config["MatLoc"] + "BANDLINES.OUT", header=None, delim_whitespace=True, names=["x", "y"])
    xvals = np.array(list(sigpoints['x']))
    yvals = np.array(list(sigpoints['y'])) * 27.2138
    length = 2
    ticks = []
    for i in range(0, len(xvals), length):
        axis.plot(xvals[i:i + length], yvals[i:i + length], label='_nolegned_', color="black", linestyle='dashed',
                  linewidth=0.5)
        ticks.append(xvals[i])

    xval = [0.8664062500,5.370500174]
    yval = ((0.2922058652 - 0.2545507825) * 27.2138)+0.0
    axis.scatter(xval, [yval, yval], c="red", s=50)

    sigpointlabels = []
    if config["kpoints"] == "MP":
        kpoints = Utils.KPath.get_kpoints_SC(mat)
    elif config["kpoints"] == "McQueen":
        kpoints = Utils.KPath.get_kpoints_McQueen(mat)
    else:
        kpoints = config["kpoints"]
    for point in kpoints:
        if point[0] == "\\Gamma":
            sigpointlabels.append("\u0393")
        else:
            sigpointlabels.append("$" + str(point[0]) + "$")
    axis.set_xticks(ticks, sigpointlabels)
    axis.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black", linewidth=0.5)  # ,label="Fermi Level")
    axis.set_xlabel("Wave Vector")
    axis.set_title(str(mat.composition.to_latex_string()) + " Band Structure")


def plot(config, mat: pymatgen.core.structure.Structure, options=None, spins=False, el_orbs=None, energy_range=(-5, 8),
         show=False, sites=None, dos_range=None, titles=True, num_points=500):
    f, (ax2, ax1) = plt.subplots(1, 2, sharey=True, figsize=(20, 10), gridspec_kw={'width_ratios': [2, 1]})
    f.subplots_adjust(wspace=0)

    if "IDOS" in options:
        plot_IDOS(config, mat, ax1, both_spins=spins, num_points=num_points)

    if "ODOS" in options:
        plot_orb_DOS(config, mat, ax1, both_spins=spins, el_orbs=el_orbs, energy_range=energy_range,
                     dos_range=dos_range, num_points=num_points)
    elif "OSDOS" in options:
        plot_orb_site_DOS(config, mat, ax1, both_spins=spins, el_orbs=el_orbs, num_points=num_points)
    elif "SDOS" in options:
        plot_site_DOS(config, mat, ax1, both_spins=spins, el_orbs=el_orbs, num_points=num_points)
    elif "OADOS" in options:
        plot_orb_ang_DOS(config, mat, ax1, both_spins=spins, el_orbs=el_orbs, energy_range=energy_range,
                         num_points=num_points)
    elif "OASDOS" in options:
        plot_orb_ang_site_DOS(config, mat, ax2, both_spins=spins, el_orbs=el_orbs, ele_sites=sites,
                              num_points=num_points)
    if "EDOS" in options:
        plot_ele_DOS(config, mat, ax1, both_spins=spins, energy_range=energy_range, dos_range=dos_range,
                     num_points=num_points)
    if "TDOS" in options:
        plot_TDOS(config, mat, ax2, both_spins=spins, num_points=num_points)

    if "SBS" in options:
        plot_spin(config, mat, ax2)
    elif "OBS" in options:
        plot_orbital_bands(config, mat, ax2, el_orbs=el_orbs, energy_range=energy_range)
    elif "2OBS" in options:
        plot_2_orbital_bands(config, mat, ax2, f, el_orbs=el_orbs, energy_range=energy_range)
    else:
        plot_bands(config, mat, ax2)

    if not titles:
        ax1.set_title("")
        ax2.set_title("")
        #ax3.set_title("")
    else:
        tmp = mat.copy()
        tmp.remove_oxidation_states()
        form = tmp.composition.reduced_composition.to_latex_string().replace("$_{1}$", "")
        ax2.set_title(form + " Band Structure")
        #ax3.set_title(form + " DOS")
        ax1.set_title(form + " DOS")
    ax2.set_ylabel("$E-E_f (eV)$")
    plt.ylim(energy_range)
    if show:
        plt.show()
    else:
        plt.savefig(config["MatLoc"] + "plt.png")
        plt.savefig("data/plots/" + mat.formula + "plt.png")


def tri_plot(configs, mats, options=None, spins=False, el_orbs=None, energy_range=(-5, 8),
             show=False, sites=None, dos_range=None, titles=True, num_points=500):
    f, axs = plt.subplots(1, 3, sharey=True, figsize=(30, 5), gridspec_kw={'width_ratios': [1, 1, 1]})
    f.subplots_adjust(wspace=0)

    for i in np.arange(3):
        plot_bands(configs[i], mats[i], axs[i])
        axs[i].set_title("")
    axs[0].set_ylabel("$E-E_f (eV)$")
    plt.ylim(energy_range)
    plt.show()
