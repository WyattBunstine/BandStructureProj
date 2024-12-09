import json
import sys
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
import main as BSP_main
import yaml
import matplotlib.pyplot as plt


HtoEv = 27.2114

def main():
    formula_LaNiO3 = [
        [[1 / 2, "data/LaNiO/La2O3.json", 2.5], [1, "data/LaNiO/NiO.json", 2], [1 / 4, "data/LaNiO/O2.json", 0.5]],
        [[1, "data/LaNiO/LaNiO3.json", 5]]]

    formula_La2NiO4 = [[[1, "data/LaNiO/La2O3.json", 5], [1, "data/LaNiO/NiO.json", 2]],
                       [[1, "data/LaNiO/La2NiO4.json", 7]]]

    formula_La3Ni2O7 = [
        [[1.5, "data/LaNiO/La2O3.json", 7.5], [2, "data/LaNiO/NiO.json", 4], [1 / 4, "data/LaNiO/O2.json", 0.5]],
        [[1, "data/LaNiO/La3Ni2O7.json", 12]]]

    formula_La4Ni3O10 = [
        [[2, "data/LaNiO/La2O3.json", 10], [3, "data/LaNiO/NiO.json", 6], [1 / 2, "data/LaNiO/O2.json", 1]],
        [[1, "data/LaNiO/La4Ni3O10.json", 17]]]

    formula_La5Ni4O13 = [
        [[2.5, "data/LaNiO/La2O3.json", 12.5], [4, "data/LaNiO/NiO.json", 8], [3 / 4, "data/LaNiO/O2.json", 1.5]],
        [[1, "data/LaNiO/La5Ni4O13_SO.json", 22]]]

    formula_La6Ni5O16 = [
        [[3, "data/LaNiO/La2O3.json", 15], [5, "data/LaNiO/NiO.json", 10], [1, "data/LaNiO/O2.json", 2]],
        [[1, "data/LaNiO/La6Ni5O16_SO.json", 27]]]

    formula = formula_La5Ni4O13

    run = False
    download = False
    if run:
        for input in formula[0]:
            with open(input[1]) as f:
                config = json.load(f)
                if not os.path.exists(config["MatLoc"]):
                    BSP_main.from_config(config)
                    print("Runnning input " + str(input[0]) + "*" + config["MatID"])
                else:
                    print("Not running " + config["MatID"] + ", path exists.")
        for output in formula[1]:
            with open(output[1]) as f:
                config = json.load(f)
                if not os.path.exists(config["MatLoc"]):
                    BSP_main.from_config(config)
                    print("Runnning output" + str(output[0]) + "*" + config["MatID"])
                else:
                    print("Not running " + config["MatID"] + ", path exists.")
    elif download:
        for input in formula[0]:
            with open(input[1]) as f:
                config = json.load(f)
                BSP_main.download_remote(config["MatLoc"])
        for output in formula[1]:
            with open(output[1]) as f:
                config = json.load(f)
                BSP_main.download_remote(config["MatLoc"])
    else:
        input_energy = 0.0
        for input in formula[0]:
            with open(input[1]) as f:
                config = json.load(f)
                info = open(config["MatLoc"] + "INFO.OUT", "r")
                energy = 0.0
                atoms = 0
                for line in info.readlines():
                    if "Total number of atoms per unit cell :" in line:
                        line = line.split()
                        atoms = int(line[-1])
                    else:
                        line = line.split()
                        if len(line) > 1:
                            if line[0] == "total" and line[1] == "energy":
                                energy = float(line[-1])
                info.close()
                eva = energy / atoms
                print(
                    config["MatID"] + f" total energy: {energy:.3f}, number of atoms: {atoms}, E/A: {eva:.3f}".format(
                        energy, atoms, eva))
                input_energy += energy * input[2] / atoms
        output_energy = 0.0
        for output in formula[1]:
            with open(output[1]) as f:
                config = json.load(f)
                info = open(config["MatLoc"] + "INFO.OUT", "r")
                energy = 0.0
                atoms = 0
                for line in info.readlines():
                    if "Total number of atoms per unit cell :" in line:
                        line = line.split()
                        atoms = int(line[-1])
                    else:
                        line = line.split()
                        if len(line) > 1:
                            if line[0] == "total" and line[1] == "energy":
                                energy = float(line[-1])
                info.close()
                eva = energy / atoms
                print(
                    config["MatID"] + f" total energy: {energy:.3f}, number of atoms: {atoms}, E/A: {eva:.3f}".format(
                        energy, atoms, eva))
                output_energy += energy * output[2] / atoms
        input_epa = input_energy / formula[1][0][2]
        output_epa = output_energy / formula[1][0][2]

        print(f"Input energy:  {input_energy:.6f}, E/A: {input_epa:.6f}".format(input_energy, input_epa))
        print(f"Output energy: {output_energy:.6f}, E/A: {output_epa:.6f}".format(output_energy, output_epa))
        #get all the energies


def plot_formation_energies():
    energies = {"$LaNiO_3$": [-2047.813075, -2047.811413], "$La_2NiO_4$": [-2686.957220, -2686.950539],
                "$La_3Ni_2O_7$": [-2420.647159, -2420.641827], "$La_4Ni_3O_{10}$": [-2310.990076, -2310.985315],
                "$La_5Ni_4O_{13}$": [-2251.177121, -2251.167008]}
    xvals = []
    yvals = []
    tick_locs = []
    tick_labels = []
    i = 0
    for material in energies.keys():
        xvals.append(i)
        yvals.append((energies[material][1] - energies[material][0])*HtoEv)
        tick_locs.append(i)
        tick_labels.append(material)
        i+=1
    plt.scatter(xvals, yvals,)
    plt.title("Energies above hull of Nicklates vs $La_2O_3 + NiO + O_2$")
    plt.xlabel("Material")
    plt.ylabel("Energy Above Hull (eV/A)")
    plt.xticks(ticks=tick_locs, labels=tick_labels)
    plt.show()



if __name__ == "__main__":
    plot_formation_energies()

    formula_KLaTiO4 = [[[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/TiO2.json", 6],
                        [1, "data/IridateProj/K2CO3.json", 6]],
                       [[2, "data/IridateProj/KLaTiO4.json", 14], [1, "data/IridateProj/CO2.json", 3]]]
    formula_KLaIrO4 = [
        [[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/IrO2.json", 6],
         [1, "data/IridateProj/K2CO3.json", 6]],
        [[2, "data/IridateProj/KLaIrO4.json", 14], [1, "data/IridateProj/CO2.json", 3]]]

    formula_La2KIrO6 = [
        [[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/IrO2.json", 6],
         [1, "data/IridateProj/K2CO3.json", 6]],
        [[1, "data/IridateProj/La2KIrO6.json", 10], [1, "data/IridateProj/KIrO3.json", 5],
         [1, "data/IridateProj/CO.json", 2]]]
    formula_La2KTiO6 = [
        [[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/TiO2.json", 6],
         [1, "data/IridateProj/K2CO3.json", 6]],
        [[1, "data/IridateProj/La2KTiO6.json", 10], [1, "data/IridateProj/KTiO3.json", 5],
         [1, "data/IridateProj/CO.json", 2]]]

    formula_NaLaIrO4 = [
        [[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/IrO2.json", 6],
         [1, "data/IridateProj/Na2CO3.json", 6]],
        [[2, "data/IridateProj/NaLaIrO4.json", 14], [1, "data/IridateProj/CO2.json", 3]]]
    formula_La2NaIrO6 = [
        [[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/IrO2.json", 6],
         [1, "data/IridateProj/Na2CO3.json", 6]],
        [[1, "data/IridateProj/La2NaIrO6.json", 10], [1, "data/IridateProj/NaIrO3.json", 5],
         [1, "data/IridateProj/CO.json", 2]]]

    formula_LiLaIrO4 = [
        [[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/IrO2.json", 6],
         [1, "data/IridateProj/Li2CO3.json", 6]],
        [[2, "data/IridateProj/LiLaIrO4.json", 14], [1, "data/IridateProj/CO2.json", 3]]]
    formula_La2LiIrO6 = [
        [[1, "data/IridateProj/La2O3.json", 5], [2, "data/IridateProj/IrO2.json", 6],
         [1, "data/IridateProj/Li2CO3.json", 6]],
        [[1, "data/IridateProj/La2LiIrO6.json", 10], [1, "data/IridateProj/LiIrO2.json", 4],
         [1, "data/IridateProj/CO2.json", 3]]]
