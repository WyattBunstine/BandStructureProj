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


orbitals = ["s", "p_y", "p_z", "p_x", "d_{xy}", "d_{xz}", "d_{z^2}", "d_{yz}", "d_{x^2-y^2}", "f1", "f2", "f3",
            "f4", "f5", "f6", "f7"]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
          '#17becf']
line_styles = ["solid", "dotted", "dashed", "dashdot"]


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


def plot_orbital_bands(config, mat: pymatgen.core.structure.Structure, axis, el_orbs=None, energy_range=(-8, 5)):
    # plot bandlines

    used_colors = 0
    atoms = sorted([el.symbol for el in mat.composition.elements])
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

            new_vals = pd.DataFrame(columns=["vec", "energy", "s", "p_y", "p_z", "p_x", "d_{xy}", "d_{yz}",
                                             "d_{z^2}", "d_{xz}", "d_{x^2-y^2}", "f1", "f2", "f3", "f4", "f5",
                                             "f6", "f7", "orb_max"])
            for index, row in bands.iterrows():
                new_vals.loc[len(new_vals.index)] = [row["vec"], row["energy"]] + list(band_vals[index]) + [
                    np.argmax(band_vals[index])]
            for i in np.arange(16):
                if orbitals[i] in el_orbs[atoms[species - 1]]:
                    orb_bands = new_vals.loc[(new_vals["orb_max"] == i)]
                    xvals = np.array(list(orb_bands['vec']))
                    xvals += species*0.001
                    yvals = np.array(list(orb_bands['energy'])) * 27.2138
                    axis.scatter(xvals, yvals, s=1, label=atoms[species - 1] + "$" + orbitals[i] + "$",
                                 c=colors[used_colors%10], zorder=10)
                    used_colors += 1
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
    orbitals = ["s", "p", "d"]
    forbitals = False
    if not el_orbs == None:
        for key in el_orbs.keys():
            if "f" in el_orbs[key] and not forbitals:
                orbitals.append("f")
                forbitals = True
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
            xvals = np.zeros(num_points*32)
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
            m = 0
            dos = np.zeros(num_points)
            for i in np.arange(0, num_points*16, num_points):
                dos += xvals[i:i + num_points]
                m += 1
                if (i < num_points and m == 1) or (i < num_points*4 and m == 3) or (i < num_points*9 and m == 5) or (m == 7 and forbitals):
                    orbital = orbitals[int(m / 2)]
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
                    m = 0
                    dos = np.zeros(num_points)
            if both_spins:
                used_colors -= species_colors
                m = 0
                dos = np.zeros(num_points)
                for i in np.arange(num_points*16, num_points*32, num_points):
                    dos += xvals[i:i + num_points]
                    m += 1
                    if (i < num_points*17 and m == 1) or (i < num_points*20 and m == 3) or (i < num_points*25 and m == 5) or (
                            m == 7 and forbitals):
                        orbital = orbitals[int(m / 2.0)]
                        element = str(atoms[species - 1])
                        if el_orbs is None or (element in el_orbs.keys() and orbital in el_orbs[element]):
                            style = int(used_colors / len(colors))
                            color = used_colors % len(colors)
                            axis.plot(dos, yvals, label="_nolabel_", color=colors[color], linestyle=line_styles[style])
                            used_colors += 1
                        if np.min(dos[range[0]:range[1]]) < x_range[0]:
                            x_range[0] = np.min(dos[range[0]:range[1]])
                        m = 0
                        dos = np.zeros(num_points)

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
                xvals = np.zeros(16*num_points)
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
                m = 0
                for i in np.arange(0, 16*num_points, num_points):
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
                    for i in np.arange(16*num_points, 32*num_points, num_points):
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
            xvals = np.zeros(16*num_points)
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
            for i in np.arange(0, num_points*16, num_points):
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
                for i in np.arange(num_points*16, num_points*32, num_points):
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
            axis.plot(xvals[num_points:num_points*2], yvals, color="black")  # label=r'$\downarrow$')
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
        plot_orb_site_DOS(config, mat, ax3, both_spins=spins, el_orbs=el_orbs, num_points=num_points)
    elif "OADOS" in options:
        plot_orb_ang_DOS(config, mat, ax1, both_spins=spins, el_orbs=el_orbs, energy_range=energy_range, num_points=num_points)
    elif "OASDOS" in options:
        plot_orb_ang_site_DOS(config, mat, ax3, both_spins=spins, el_orbs=el_orbs, ele_sites=sites, num_points=num_points)
    if "EDOS" in options:
        plot_ele_DOS(config, mat, ax1, both_spins=spins, energy_range=energy_range, dos_range=dos_range, num_points=num_points)
    if "TDOS" in options:
        plot_TDOS(config, mat, ax3, both_spins=spins, num_points=num_points)

    if "SBS" in options:
        plot_spin(config, mat, ax2)
    elif "OBS" in options:
        plot_orbital_bands(config, mat, ax2, el_orbs=el_orbs, energy_range=energy_range)
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
        ax3.set_title(form + " DOS")
        ax1.set_title(form + " DOS")
    ax2.set_ylabel("$E-E_f (eV)$")
    plt.ylim(energy_range)
    if show:
        plt.show()
    else:
        plt.savefig(config["MatLoc"] + "plt.png")
        plt.savefig("data/plots/" + mat.formula + "plt.png")


'''
            if i == 76000:
                axis.plot(xvals, yvals[i:i + 1000], label='parity @ Z = 1', color="blue", linewidth=1.0)
            elif i == 77000:
                axis.plot(xvals, yvals[i:i + 1000], label='_nolegned_', color="blue", linewidth=1.0)
            elif i == 78000:
                axis.plot(xvals, yvals[i:i + 1000], label='parity @ Z = -1', color="red", linewidth=1.0)
            elif i == 79000:
                axis.plot(xvals, yvals[i:i + 1000], label='_nolegned_', color="red", linewidth=1.0)
            else:
                axis.plot(xvals, yvals[i:i + 1000], linewidth=1.0, color="black", label='_nolegned_')
            '''
