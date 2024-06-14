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
import yaml
import sys


def main():


    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    starting_mat = "mp-1178424"
    starting_mat = "mp-550306"
    starting_mat = "mp-1077840"
    kept_site = ["Bi", "Te"]

    neighbor_dist = 4

    with MPRester(API_KEY) as mpr:
        starting_mat = mpr.materials.get_structure_by_material_id(starting_mat)

    sites = []
    for site in starting_mat.sites:
        neighbors = starting_mat.get_neighbors(site, neighbor_dist)
        neighbor_count = {}
        for neighbor in neighbors:
            if neighbor.specie is not site.specie:
                if neighbor.specie in neighbor_count.keys():
                    neighbor_count[neighbor.specie] +=1
                else:
                    neighbor_count[neighbor.specie] = 1
        neighbor_max = "H"
        neighbor_max_num = 1
        for neighbor in neighbor_count.keys():
            if neighbor_count[neighbor] > neighbor_max_num:
                neighbor_max_num = neighbor_count[neighbor]
                neighbor_max = neighbor
        sites.append([site.specie, neighbor_max])

    elements = {}
    for site in sites:
        if site[0] in elements.keys():
            elements[site[0]].append(site[1])
        else:
            elements[site[0]] = [site[1]]
    kept_pairs = []
    pairs = []
    non_max_pairs = []
    for element in elements.keys():
        if element.symbol not in kept_site:
            pairs.append([element,max(set(elements[element]), key=elements[element].count)])
        else:
            kept_pairs.append([element, max(set(elements[element]), key=elements[element].count)])
        for el in elements[element]:
            if el != max(set(elements[element]), key=elements[element].count):
                non_max_pairs.append([element, el])

    print(sites)
    print(elements)
    print(kept_pairs)
    print(pairs)

    nonmax_pairs_FE = []
    for pair in non_max_pairs + kept_pairs + pairs:
        if pair[0].symbol != pair[1].symbol:
            possible_mats = mpr.materials.summary.search(elements=[pair[0].symbol, pair[1].symbol],
                                                         is_stable=True, all_fields=False,
                                                         fields=["composition", "material_id",
                                                                 "formation_energy_per_atom"],
                                                         num_elements=2)

            formation_energy = 0.0
            best_comp = ""
            for possible_mat in possible_mats:
                if formation_energy > possible_mat.formation_energy_per_atom:
                    best_comp = possible_mat.composition
                    formation_energy = possible_mat.formation_energy_per_atom
            if len(possible_mats) > 0:
                nonmax_pairs_FE.append([best_comp, formation_energy])
    print("Kept pairs formation energies:")
    print(nonmax_pairs_FE)
    print("\n\n")
    nonmax_pairs_FE = []



    percent_agree = 0.25

    for pair in pairs:
        replacement_options = []
        pair_el = pair[0]
        for z in np.arange(1, 85):
            el = pymatgen.core.Element.from_Z(z)
            if el.atomic_radius is not None and el.number != pair_el.number and el.X is not None:
                if (abs(el.atomic_radius - pair_el.atomic_radius)/(2*pair_el.atomic_radius) < percent_agree and
                        abs(el.electron_affinity - pair_el.electron_affinity) / (2 * pair_el.electron_affinity) < percent_agree and
                    abs(el.X - pair_el.X) / (2 * pair_el.X) < percent_agree and
                    abs(el.ionization_energy - pair_el.ionization_energy) / (2 * pair_el.ionization_energy) < percent_agree):
                    print("\n\npossible substitution: " + str(el))
                    with MPRester(API_KEY) as mpr:
                        possible_mats = mpr.materials.summary.search(elements=[el.symbol,pair[1].symbol],
                                                                     is_stable=True,all_fields=False,
                                                                     fields=["composition", "material_id", "formation_energy_per_atom"],
                                                                     num_elements=2)

                        formation_energy = 0.0
                        best_comp = ""
                        for possible_mat in possible_mats:
                            if formation_energy > possible_mat.formation_energy_per_atom:
                                best_comp = possible_mat.composition
                                formation_energy = possible_mat.formation_energy_per_atom
                        if best_comp != "":
                            print("lowest formation energy per atom composition: " + best_comp.reduced_formula + " with FE/A: " + str(formation_energy))
                        #kept_pairs_FE = nonmax_pairs_FE
                        kept_pairs_FE = []
                        for kept_pair in kept_pairs:
                            if kept_pair[1].symbol != pair[1].symbol:
                                possible_mats = mpr.materials.summary.search(elements=[el.symbol, kept_pair[1].symbol],
                                                                             is_stable=True, all_fields=False,
                                                                             fields=["composition", "material_id",
                                                                                     "formation_energy_per_atom"],
                                                                             num_elements=2)

                                kept_formation_energy = 0.0
                                kept_best_comp = ""
                                for possible_mat in possible_mats:
                                    if kept_formation_energy > possible_mat.formation_energy_per_atom:
                                        kept_best_comp = possible_mat.composition
                                        kept_formation_energy = possible_mat.formation_energy_per_atom
                                if len(possible_mats) > 0:
                                    kept_pairs_FE.append([kept_best_comp,kept_formation_energy])
                        allow = True
                        for formation_pair in kept_pairs_FE:
                            if formation_energy < formation_pair[1]:
                                print(best_comp.reduced_formula + " will not work: formation energy/a: " + str(formation_energy))
                                print("but " + formation_pair[0].reduced_formula + " formation energy/a: " + str(formation_pair[1]))
                                allow = False
                        if allow and best_comp != "":
                            print(best_comp.reduced_formula + " may work: formation energy/a: " + str(formation_energy))
                            replacement_options.append(el)
        print("replacement options for " + str(pair_el) + ": " + str(replacement_options))


if __name__ == "__main__":
    main()
