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
import pandas as pd
import matplotlib.pyplot as plt


def find_fermi_states(config):
    nkpt = 0
    nstate = 0
    kpoints = {}
    if os.path.exists(config["MatLoc"] + "EIGVAL.OUT"):
        with open(config["MatLoc"] + "EIGVAL.OUT", "r") as f:
            for line in f.readlines():
                if "nkpt" in line:
                    nkpt = int(line[:line.find(":")].strip(" "))
                elif "nstsv" in line:
                    nstate = int(line[:line.find(":")].strip(" "))
                elif "k-point" in line:
                    pt = line[:line.find(":")].split()
                    lastkpt = int(pt[0])
                    kpoints[int(pt[0])] = {"coords": [float(pt[1]), float(pt[2]), float(pt[3])], "eigenvalues": {}}
                elif len(line) > 10 and "state" not in line:
                    eigenval = line.split()
                    kpoints[lastkpt]["eigenvalues"][int(eigenval[0])] = [float(eigenval[1]), float(eigenval[2])]
    else:
        warn(config["MatLoc"] + "EIGVAL.OUT not found.")
        return -1
    fermiStates = {}
    for kpt in np.arange(1, nkpt):
        added = False
        states = kpoints[kpt]["eigenvalues"]
        for state in np.arange(1, nstate):
            if 0.0 < states[state][1] < 1.0:
                if not added:
                    added = True
                    fermiStates[kpt] = {"coords": kpoints[kpt]["coords"], "eigenvalues": []}
                fermiStates[kpt]["eigenvalues"].append([state, states[state][0], states[state][1]])

    return fermiStates


def find_band_crossings(config, mat):
    if config["kpoints"] == "MP":
        kpoints = Utils.KPath.get_kpoints_SC(mat)
    elif config["kpoints"] == "McQueen":
        kpoints = Utils.KPath.get_kpoints_McQueen(mat)
    else:
        kpoints = config["kpoints"]

    fermi_states = find_fermi_states(config)
    for i in range(len(kpoints) - 1):
        pt1 = kpoints[i]
        pt2 = kpoints[i + 1]
        coord = -1
        if pt1[1][0] == pt2[1][0] and pt1[1][1] == pt2[1][1]:
            coord = 2
        elif pt1[1][2] == pt2[1][2] and pt1[1][1] == pt2[1][1]:
            coord = 0
        elif pt1[1][2] == pt2[1][2] and pt1[1][0] == pt2[1][0]:
            coord = 1
        if not coord == -1:
            for state in fermi_states.keys():
                coords = fermi_states[state]["coords"]
                if (coord == 0 and coords[1] == pt1[1][1] and coords[2] == pt1[1][2]) or \
                        (coord == 1 and coords[0] == pt1[1][0] and coords[2] == pt1[1][2]) or \
                        (coord == 2 and coords[1] == pt1[1][1] and coords[0] == pt1[1][0]):
                    if pt1[0] == "\\Gamma":
                        pt1[0] = "\u0393"
                    if pt2[0] == "\\Gamma":
                        pt2[0] = "\u0393"
                    print("Band crossing found between " + pt1[0] + str(pt1[1]) + " and " + pt2[0] + str(pt2[1]))
                    print("kpt: " + str(state) + ", coords: " + str(coords))
                    print("state, occupancy")
                    for eigenval in fermi_states[state]["eigenvalues"]:
                        print(str(eigenval[0]) + ": " + str(eigenval[2]))
