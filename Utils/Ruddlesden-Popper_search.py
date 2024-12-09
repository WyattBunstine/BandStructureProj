import mp_api
import numpy as np
import pandas
from mp_api.client import MPRester
import pymatgen
from pymatgen.core import *
from pymatgen.analysis import *
from pymatgen.analysis.chemenv.coordination_environments import coordination_geometry_finder
from pymatgen import *
import json
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from sklearn import neural_network
from sklearn.cluster import KMeans
import json
import time
import multiprocessing
from monty.serialization import loadfn


def customwarn(message, category, filename, lineno, file=None, line=None):
    1 + 1  # sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))


HtoEv = 27.2114  # 1 Hartree is 27 Electron volts


def coordination_filter(mat: pymatgen.core.Structure, data):
    columns = ["n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH"]
    GF = coordination_geometry_finder.LocalGeometryFinder()
    GF.setup_structure(mat)
    # Get the StructureEnvironments
    ces = GF.compute_coordination_environments(mat, only_cations=False)
    for i in range(len(mat.sites)):
        if data[0] == 0:
            site = mat.sites[i]
            ce = ces[i]
            if site.specie.symbol == data[1]:
                ce_12_coord = False
                for possible_env in ce:
                    if ("12" in possible_env["ce_symbol"] or "11" in possible_env["ce_symbol"] or "10" in
                            possible_env["ce_symbol"] or "9" in possible_env["ce_symbol"]
                            or "8" in possible_env["ce_symbol"]):
                        ce_12_coord = True
                if not ce_12_coord:
                    print(mat.composition.reduced_composition.to_pretty_string() + " is not a perovskite")
                    return False
            if site.specie.symbol == data[2]:
                ce_6_coord = False
                for possible_env in ce:
                    if "6" in possible_env["ce_symbol"]:
                        ce_6_coord = True
                if not ce_6_coord:
                    print(mat.composition.reduced_composition.to_pretty_string() + " is not a perovskite")
                    return False
        else:
            site = mat.sites[i]
            ce = ces[i]
            #print(site)
            #print(ce)
            if site.specie.symbol == data[1]:
                ce_12_coord = False
                ce_9_coord = False
                for possible_env in ce:
                    if "12" in possible_env["ce_symbol"]:
                        ce_12_coord = True
                    if "9" in possible_env["ce_symbol"]:
                        ce_9_coord = True
                if not ce_12_coord and not ce_9_coord:
                    print(mat.composition.reduced_composition.to_pretty_string() + " is not a Ruddlesden-Popper")
                    return False
            if site.specie.symbol == data[2]:
                ce_6_coord = False
                for possible_env in ce:
                    if "6" in possible_env["ce_symbol"]:
                        ce_6_coord = True
                if not ce_6_coord:
                    print(mat.composition.reduced_composition.to_pretty_string() + " is not a Ruddlesden-Popper")
                    return False

    if data[0] == 0:
        print(mat.composition.reduced_composition.to_pretty_string() + " is a perovskite")
    else:
        print(mat.composition.reduced_composition.to_pretty_string() + " is a Ruddlesden-Popper")
    return True


def get_df_column(mat):
    columns = "n", "A", "B", "X", "Aox", "Box", "Cox", "Ra", "Rb", "Rc", "exp", "EAH"
    n = 0
    A = ""
    B = ""
    X = ""
    comp = mat.composition.reduced_composition.remove_charges()
    elements = comp.as_dict()
    for el in elements.keys():
        tmp = pymatgen.core.periodic_table.get_el_sp(el)
        if tmp.is_chalcogen or tmp.is_halogen:
            n = (elements[el] - 1) / 3
            X = el
            if not n.is_integer():
                n = 0
                # print("n problem: " + str(n) + " " + str(type(n)))
    for el in elements.keys():
        tmp = pymatgen.core.periodic_table.get_el_sp(el)
        if not (tmp.is_chalcogen or tmp.is_halogen):
            if (elements[el] == n or n == 0) and B == "":
                B = el
            else:
                A = el
    if pymatgen.core.periodic_table.get_el_sp(B).atomic_radius > pymatgen.core.periodic_table.get_el_sp(
            A).atomic_radius:
        tmp = A
        A = B
        B = tmp

    if A == "" or B == "" or X == "":
        print("Err with: " + mat.composition)

    data = [n, A, B, X]
    oxi_state = [0, 0, 0, 0, 0, 0]

    oxi_fail = False

    mat_is_struct = False

    if type(mat) is not pymatgen.core.Structure:
        struct = mat.structure
    else:
        mat_is_struct = True
        struct = mat

    try:
        BVAnalyzer = pymatgen.analysis.bond_valence.BVAnalyzer()
        struct = BVAnalyzer.get_oxi_state_decorated_structure(struct)
        for site in struct:
            1  #print(site)
    except Exception:
        oxi_fail = True

    A_sites = 0
    B_sites = 0
    X_sites = 0

    if oxi_fail:
        A_sites = 1
        B_sites = 1
        X_sites = 1
        for site in struct.sites:
            if site.specie.symbol == A:
                oxi_state[3] = site.specie.atomic_radius
            elif site.specie.symbol == B:
                oxi_state[4] = site.specie.atomic_radius
            elif site.specie.symbol == X:
                oxi_state[5] = site.specie.atomic_radius
    else:
        for site in struct.sites:
            if site.specie.symbol == A:
                A_sites += 1
                if oxi_state[0] != 0:
                    oxi_state[0] = (site.specie.oxi_state + oxi_state[0])
                else:
                    oxi_state[0] = site.specie.oxi_state
                rad = 0
                if site.specie.oxi_state in site.specie.ionic_radii.keys():
                    rad = float(site.specie.ionic_radii[int(site.specie.oxi_state)])
                else:
                    rad = float(site.specie.atomic_radius) - site.specie.oxi_state * 0.1
                oxi_state[3] = oxi_state[3] + rad


            elif site.specie.symbol == B:
                B_sites += 1
                if oxi_state[1] != 0:
                    oxi_state[1] = (site.specie.oxi_state + oxi_state[1])
                else:
                    oxi_state[1] = site.specie.oxi_state
                rad = 0
                if site.specie.oxi_state in site.specie.ionic_radii.keys():
                    rad = float(site.specie.ionic_radii[int(site.specie.oxi_state)])
                else:
                    print(site)
                    rad = float(site.specie.atomic_radius) - site.specie.oxi_state * 0.1
                oxi_state[4] = oxi_state[4] + rad


            elif site.specie.symbol == X:
                X_sites += 1
                if oxi_state[2] != 0:
                    oxi_state[2] = (site.specie.oxi_state + oxi_state[2])
                else:
                    oxi_state[2] = site.specie.oxi_state
                rad = 0
                if site.specie.oxi_state in site.specie.ionic_radii.keys():
                    rad = float(site.specie.ionic_radii[int(site.specie.oxi_state)])
                else:
                    rad = float(site.specie.atomic_radius) - site.specie.oxi_state * 0.1
                oxi_state[5] = oxi_state[5] + rad

    oxi_state[0] = oxi_state[0] / A_sites
    oxi_state[1] = oxi_state[1] / B_sites
    oxi_state[2] = oxi_state[2] / X_sites
    oxi_state[3] = oxi_state[3] / A_sites
    oxi_state[4] = oxi_state[4] / B_sites
    oxi_state[5] = oxi_state[5] / X_sites

    data = data + oxi_state
    if not mat_is_struct and mat.theoretical:
        data.append(False)
    else:
        data.append(True)
    if not mat_is_struct and mat.energy_above_hull is not None:
        data.append(mat.energy_above_hull)
    else:
        data.append(0.0)
    # print(data)
    return data


def get_df_column_shannon(mat, theoretical, EAH):
    columns = "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH"
    n = 0
    A = ""
    B = ""
    X = ""
    comp = mat.composition.reduced_composition.remove_charges()
    elements = comp.as_dict()
    maxn = 0
    for el in elements.keys():
        if elements[el] > maxn:
            maxn = elements[el]
            X = el
    if maxn == 3:
        n = 0
    else:
        n = int((maxn - 1) / 3)

    for el in elements.keys():
        tmp = pymatgen.core.periodic_table.get_el_sp(el)
        if (elements[el] == n or (n == 0 and elements[el] == 1)) and B == "":
            B = el
        elif (elements[el] == n + 1) and A == "":
            A = el
        elif el != X:
            print(el)
            print("Err with: " + mat.composition.reduced_composition.to_pretty_string())
            return [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    if A == "" or B == "" or X == "":
        print("Err with: " + mat.composition.reduced_composition.to_pretty_string())
        return [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    data = [n, A, B, X]
    oxi_state = [0, 0, 0, 0, 0, 0, 0]
    oxi_fail = False

    try:
        BVAnalyzer = pymatgen.analysis.bond_valence.BVAnalyzer()
        struct = BVAnalyzer.get_oxi_state_decorated_structure(mat)
    except Exception:
        oxi_fail = True
        struct = mat

    A_sites = 0
    B_sites = 0
    X_sites = 0

    if oxi_fail:
        A_sites = 1
        B_sites = 1
        X_sites = 1
        for site in struct.sites:
            if site.specie.symbol == A:
                try:
                    if 2 in site.specie.element.common_oxidation_states:
                        site.specie.oxi_state = 2
                    elif 3 in site.specie.element.common_oxidation_states:
                        site.specie.oxi_state = 3
                    else:
                        site.specie.oxi_state = 1
                    rad9 = float(site.specie.get_shannon_radius("IX"))
                    rad12 = float(site.specie.get_shannon_radius("XII"))
                except Exception:
                    rad9 = float(site.specie.atomic_radius)
                    rad12 = float(site.specie.atomic_radius)
                oxi_state[3] = rad9
                oxi_state[4] = rad12
            elif site.specie.symbol == B:
                oxi_state[5] = site.specie.atomic_radius
            elif site.specie.symbol == X:
                oxi_state[6] = site.specie.atomic_radius
    else:
        for site in struct.sites:
            if site.specie.symbol == A:
                A_sites += 1
                if oxi_state[0] != 0:
                    oxi_state[0] = (site.specie.oxi_state + oxi_state[0])
                else:
                    oxi_state[0] = site.specie.oxi_state
                rad9 = 0
                rad12 = 0
                try:
                    rad9 = float(site.specie.get_shannon_radius("IX"))
                except Exception:
                    rad9 = 0
                try:
                    rad12 = float(site.specie.get_shannon_radius("XII"))
                except Exception:
                    rad12 = 0
                oxi_state[3] = oxi_state[3] + rad9
                oxi_state[4] = oxi_state[4] + rad12


            elif site.specie.symbol == B:
                B_sites += 1
                if oxi_state[1] != 0:
                    oxi_state[1] = (site.specie.oxi_state + oxi_state[1])
                else:
                    oxi_state[1] = site.specie.oxi_state
                rad = 0
                try:
                    if site.specie.element.Z < 27:
                        rad = float(site.specie.get_shannon_radius("VI", "Low Spin"))
                    else:
                        rad = float(site.specie.get_shannon_radius("VI", "High Spin"))
                except Exception:
                    #print(site)
                    rad = float(site.specie.atomic_radius) - site.specie.oxi_state * 0.1
                oxi_state[5] = oxi_state[5] + rad


            elif site.specie.symbol == X:
                X_sites += 1
                if oxi_state[2] != 0:
                    oxi_state[2] = (site.specie.oxi_state + oxi_state[2])
                else:
                    oxi_state[2] = site.specie.oxi_state
                rad = 0
                try:
                    rad = float(site.specie.get_shannon_radius("VI"))
                except Exception:
                    #print(site)
                    #print(mat)
                    rad = float(site.specie.atomic_radius) - site.specie.oxi_state * 0.1
                oxi_state[6] = oxi_state[6] + rad

    oxi_state[0] = oxi_state[0] / A_sites
    oxi_state[1] = oxi_state[1] / B_sites
    oxi_state[2] = oxi_state[2] / X_sites
    oxi_state[3] = oxi_state[3] / A_sites
    oxi_state[4] = oxi_state[4] / A_sites
    oxi_state[5] = oxi_state[5] / B_sites
    oxi_state[6] = oxi_state[6] / X_sites

    data = data + oxi_state
    data.append(theoretical)
    data.append(EAH)

    # print(data)
    return data


def search_all():
    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    with MPRester(API_KEY) as mpr:
        mats = mpr.materials.summary.search(num_elements=(3, 3), all_fields=False,
                                                 formula=["A2BX4", "A3B2X7", "A4B3X10", "A5B4X13", "A6B5X16",
                                                          "A7B6X19", "ABX3"],
                                                 fields=["composition", "material_id", "structure",
                                                         "energy_above_hull",
                                                         "theoretical"])

        df = pandas.DataFrame(columns=["struct", "exp", "EAH"])
        #known = pd.read_csv("Ruddlesden-Popper_expanded_data_cleaned.csv")
        #mat_comps = []
        #for index, row in known.iterrows():
        #    mat_comps.append(pymatgen.core.Structure.from_dict(json.loads(row["struct"])).composition.reduced_composition.to_pretty_string())

        num_structs = 0
        for i in range(len(mats)):
            mat = mats[i]
            #if mat.composition.reduced_composition.to_pretty_string() not in mat_comps:
            df.loc[num_structs] = [json.dumps(mat.structure.as_dict()), mat.theoretical, mat.energy_above_hull]
            num_structs += 1
        df.to_csv("All_mats_RP_Stoich.csv", index=False)

def search():
    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    with MPRester(API_KEY) as mpr:
        perovs = []
        for spacegroup in [12, 71, 63, 69, 65, 47, 74, 139, 140, 127, 123, 221, 204, 167, 137, 62, 15, 11, 2, 161]:
            perovs = perovs + mpr.materials.summary.search(spacegroup_number=spacegroup, num_elements=(3, 3),
                                                           formula=["ABX3"],
                                                           theoretical=False, all_fields=False,
                                                           fields=["composition", "material_id", "structure",
                                                                   "energy_above_hull", "theoretical"])
            #data = data + mpr.materials.summary.search(spacegroup_number=spacegroup, num_elements=(3, 3),
            #                                           theoretical=True, energy_above_hull=(0.0, 0.25), all_fields=False,
            #                                           fields=["composition", "material_id", "structure",
            #                                                   "energy_above_hull", "theoretical"])

        print("number of experimentally observed perovskites: " + str(len(perovs)))
        RPS = []
        for spacegroup in [139, 64, 138, 134, 66, 142, 56, 48, 61, 14, 13, 86, 63, 12, 71, 69, 65, 47, 74, 140, 127,
                           123, 221, 36, 68, 136, 60, 204, 167, 137, 62, 15, 11, 2, 55]:
            RPS = RPS + mpr.materials.summary.search(spacegroup_number=spacegroup, num_elements=(3, 3),
                                                     theoretical=False,
                                                     all_fields=False,
                                                     formula=["A2BX4", "A3B2X7", "A4B3X10", "A5B4X13", "A6B5X16",
                                                              "A7B6X19"],
                                                     fields=["composition", "material_id", "structure",
                                                             "energy_above_hull",
                                                             "theoretical"])
        print("number of possible Ruddleden-Poppers: " + str(len(RPS)))

        df = pandas.DataFrame(columns=["struct", "exp", "EAH"])

        num_structs = 0
        for i in range(len(perovs)):
            mat = perovs[i]
            df.loc[num_structs] = [json.dumps(mat.structure.as_dict()), mat.theoretical, mat.energy_above_hull]
            num_structs += 1
        for i in range(len(RPS)):
            mat = RPS[i]
            df.loc[num_structs] = [json.dumps(mat.structure.as_dict()), mat.theoretical, mat.energy_above_hull]
            num_structs += 1
        df.to_csv("All_mats.csv", index=False)


def old_search():
    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    with MPRester(API_KEY) as mpr:
        data = []
        for spacegroup in [12, 71, 63, 69, 65, 47, 74, 139, 140, 127, 123, 221, 204, 167, 137, 62, 15, 11, 2]:
            data = data + mpr.materials.summary.search(spacegroup_number=spacegroup, num_elements=(3, 3),
                                                       theoretical=False, all_fields=False,
                                                       fields=["composition", "material_id", "structure",
                                                               "energy_above_hull", "theoretical"])
            #data = data + mpr.materials.summary.search(spacegroup_number=spacegroup, num_elements=(3, 3),
            #                                           theoretical=True, energy_above_hull=(0.0, 0.25), all_fields=False,
            #                                           fields=["composition", "material_id", "structure",
            #                                                   "energy_above_hull", "theoretical"])
        i = 0
        perovs = []
        perov_forms = []
        for mat in data:
            perov = True
            X = False
            AB = 0
            elements = mat.composition.reduced_composition.as_dict()
            for el in elements.keys():
                tmp = pymatgen.core.periodic_table.get_el_sp(el)
                if tmp.is_chalcogen or tmp.is_halogen:
                    if elements[el] != 3:
                        perov = False
                    else:
                        X = True
                else:
                    if elements[el] != 1:
                        perov = False
                    else:
                        AB += 1
            if perov and X and AB == 2 and mat.composition.reduced_composition.to_pretty_string() not in perov_forms:
                perov_forms.append(mat.composition.reduced_composition.to_pretty_string())
                # print(mat.composition.to_pretty_string())
                perovs.append(mat)
                i += 1

        print("number of experimentally observed perovskites: " + str(i))
        fams = []
        for perov in perovs:
            family_members = []
            elements = [str(el) for el in perov.composition.elements]
            data = []
            for spacegroup in [139, 64, 138, 134, 66, 142, 56, 48, 61, 14, 13, 86, 63, 12, 71, 69, 65, 47, 74, 140, 127,
                               123, 221, 36, 68, 136, 60, 204, 167, 137, 62, 15, 11, 2, 55]:
                data = data + mpr.materials.summary.search(spacegroup_number=spacegroup, num_elements=(3, 3),
                                                           elements=elements,
                                                           all_fields=False,
                                                           fields=["composition", "material_id", "structure",
                                                                   "energy_above_hull",
                                                                   "theoretical"])
            member_forms = []
            for mat in data:
                fam = True
                n = 0

                elements = mat.composition.reduced_composition.as_dict()
                for el in elements.keys():
                    tmp = pymatgen.core.periodic_table.get_el_sp(el)
                    if tmp.is_chalcogen or tmp.is_halogen:
                        n = (elements[el] - 1) / 3
                        if not n.is_integer():
                            fam = False
                            #print("n problem: " + str(n) + " " + str(type(n)))
                foundn = False
                if fam:
                    for el in elements.keys():
                        tmp = pymatgen.core.periodic_table.get_el_sp(el)
                        if not (tmp.is_chalcogen or tmp.is_halogen):
                            if elements[el] == n:
                                if foundn:
                                    #print("foundn problem")
                                    fam = False
                                foundn = True
                            if elements[el] != n and elements[el] != n + 1:
                                fam = False
                                #print("n problem " + el + " " + str(elements[el]) + " n: " + str(n))
                if fam and mat.composition.reduced_composition.to_pretty_string() not in member_forms:
                    member_forms.append(mat.composition.reduced_composition.to_pretty_string())
                    family_members.append(mat)
                else:
                    1  #print(mat.composition.to_pretty_string())

            fams.append(family_members)

        df = pandas.DataFrame(columns=["struct", "exp", "EAH"])

        num_structs = 0
        for i in range(len(perovs)):
            mat = perovs[i]
            df.loc[num_structs] = [json.dumps(mat.structure.as_dict()), mat.theoretical, mat.energy_above_hull]
            num_structs += 1
            for member in fams[i]:
                df.loc[num_structs] = [json.dumps(member.structure.as_dict()), member.theoretical,
                                       member.energy_above_hull]
                num_structs += 1
        df.to_csv("Ruddlesden-Popper_all_sg_exp_structures.csv", index=False)


def make_raw_csv():
    rawdf = pd.read_csv("Ruddlesden-Popper_expanded.csv")

    df = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                 "isperov"])

    start = time.time()
    num_structs = 0
    for index, row in rawdf.iterrows():
        try:
            mat = pymatgen.core.Structure.from_dict(json.loads(row["struct"]))
            data = get_df_column_shannon(mat, row["exp"], row["EAH"])
            df.loc[num_structs] = [json.dumps(mat.as_dict())] + data + [coordination_filter(mat, data)]
            num_structs += 1

            if index % 25 == 0:
                current = time.time()
                elapsed = start - current
                projection = elapsed * len(rawdf) / index
                print("\n\npercent complete: " + str(round(100 * index / len(rawdf), 2)) + "%, elapsed time: "
                      + str(round(elapsed / 60, 2)) + "Minutes. Projected time remained: " + str(
                    round(projection / 60, 2)) + "Minutes.\n\n")
        except Exception:
            1 + 1
    df.to_csv("Ruddlesden-Popper_expanded_data.csv", index=False)



def make_raw_csv_all():
    rawdf = pd.read_csv("All_mats_nonRP.csv")

    df = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                 "isperov"])

    start = time.time()
    num_structs = 0
    for index, row in rawdf.iterrows():
        try:
            mat = pymatgen.core.Structure.from_dict(json.loads(row["struct"]))
            data = get_df_column_shannon(mat, row["exp"], row["EAH"])
            df.loc[num_structs] = [json.dumps(mat.as_dict())] + data + [False]
            num_structs += 1

            if index % 25 == 0:
                current = time.time()
                elapsed = start - current
                projection = elapsed * len(rawdf) / index
                print("\n\npercent complete: " + str(round(100 * index / len(rawdf), 2)) + "%, elapsed time: "
                      + str(round(elapsed / 60, 2)) + "Minutes. Projected time remained: " + str(
                    round(projection / 60, 2)) + "Minutes.\n\n")
        except Exception:
            1 + 1
    df.to_csv("Ruddlesden-Popper_expanded_data.csv", index=False)


def Proc_Basic_Batch(df, queue, proc_num):
    warnings.showwarning = customwarn
    t1 = time.time()
    outdf = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                 "isperov"])

    start = time.time()
    num_structs = 0
    for index, row in df.iterrows():
        try:
            mat = pymatgen.core.Structure.from_dict(json.loads(row["struct"]))
            data = get_df_column_shannon(mat, row["exp"], row["EAH"])
            outdf.loc[num_structs] = [json.dumps(mat.as_dict())] + data + [coordination_filter(mat, data)]
            num_structs += 1
        except Exception:
            1 + 1

    queue.put(outdf)
    outdf.to_csv("Ruddlesden-Popper_expanded_data_proc_" + str(proc_num) + ".csv", index=False)
    print("subprocess " + str(proc_num) + " time: " + str(round(time.time() - t1, 4)) + " seconds.")
    return


def make_raw_csv_parallel():
    warnings.showwarning = customwarn
    rawdf = pd.read_csv("Ruddlesden-Popper_expanded.csv")
    rawdf = rawdf.sample(frac=1).reset_index(drop=True)
    num_processes = 48
    batch_size = len(rawdf) / num_processes

    df = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                 "isperov"])

    start = time.time()
    t1 = time.time()
    sub_procs = []
    sub_frames = multiprocessing.Queue()
    i = 0
    for sub_frame in np.array_split(rawdf, int(len(rawdf) / batch_size)):
        i += 1
        p = multiprocessing.Process(target=Proc_Basic_Batch,
                                    args=(sub_frame, sub_frames, i,))
        p.start()
        sub_procs.append(p)
    print("time to start " + str(i) + " subprocesses: " + str(round(time.time() - t1, 1)) + " seconds")
    t1 = time.time()
    while sub_frames.qsize() < i:
        time.sleep(0.5)
    print("time waiting for subprocesses " + str(round(time.time() - t1, 1)) + " seconds")
    t1 = time.time()
    while not sub_frames.empty():
        if len(df) == 0:
            df = sub_frames.get()
        else:
            df = pd.concat([df, sub_frames.get()], ignore_index=True)
    print("time constructing full frame " + str(round(time.time() - t1, 1)) + " seconds")
    df.columns = ["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                  "isperov"]
    df.to_csv("Ruddlesden-Popper_expanded_data.csv", index=False)
    t1 = time.time()
    for p in sub_procs:
        p.join(timeout=100)
    print("time terminating subprocesses " + str(round(time.time() - t1, 1)) + " seconds")
    df.columns = ["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                  "isperov"]
    df.to_csv("Ruddlesden-Popper_expanded_data_1.csv", index=False)

    print("database constructions time: " + str(round(time.time() - start, 1)))
    return 1
    df = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                 "isperov"])
    for file in os.listdir():
        if "_data_proc_" in file:
            if len(df) == 0:
                df = pd.read_csv(file)
            else:
                df = pd.concat([df, pd.read_csv(file)], ignore_index=True)
    df.to_csv("Ruddlesden-Popper_expanded_data.csv")


def get_nn_input(row):
    feat = []
    #feat.append(row["n"])
    #feat.append(row["isperov"])
    feat.append(row["RaIX"])
    feat.append(row["RaXII"])
    feat.append(row["Rb"])
    feat.append(row["Rc"])

    elA = pymatgen.core.periodic_table.get_el_sp(row["A"])
    elB = pymatgen.core.periodic_table.get_el_sp(row["B"])
    elX = pymatgen.core.periodic_table.get_el_sp(row["X"])

    feat.append(elA.X)
    feat.append(elB.X)
    feat.append(elX.X)
    #feat.append(0)
    #feat.append(0)
    #feat.append(0)

    #poly = PolynomialFeatures(2, interaction_only=True)
    #feat = poly.fit_transform([feat])[0]

    return feat


def plot():

    #"",3,Ce,Zn,O,4.0,2.0,-2.0,1.36,1.0,0.74,1.4,False,0.5098736521428568,True

    show_theoretical = False
    show_non_perov = True
    labels = True

    raw_df = pd.read_csv("Ruddlesden-Popper_expanded_data_cleaned.csv")

    df = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                 "isperov"])
    maxn = {}
    neq0 = []
    for index, row in raw_df.iterrows():
        name = row["A"] + row["B"] + row["X"]
        if row["n"] == 0:
            neq0.append(name)

        if name not in maxn.keys() and (not row["isperov"] or row["exp"]):
            maxn[name] = 0
        if row["isperov"] and name not in maxn.keys() and (not row["exp"]):
            maxn[name] = row["n"]
        elif row["isperov"] and maxn[name] < row["n"] and (not row["exp"]):
            maxn[name] = row["n"]

    added = []
    for index, row in raw_df.iterrows():
        name = str(row["n"]) + row["A"] + row["B"] + row["X"]
        if name not in added:
            df.loc[len(df)] = list(row)
            added.append(name)

    list_data = []
    list_target = []

    ICSD_states = {}
    module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
    all_data = loadfn(f"{module_dir}/../venv/lib/python3.11/site-packages/pymatgen/analysis/icsd_bv.yaml")
    total_data = {}
    for sp, data in all_data["occurrence"].items():
        if Species.from_str(sp).element in total_data.keys():
            total_data[Species.from_str(sp).element] += data
        else:
            total_data[Species.from_str(sp).element] = data
        if Species.from_str(sp).element in ICSD_states.keys():
            ICSD_states[Species.from_str(sp).element].append(Species.from_str(sp).oxi_state)
        else:
            ICSD_states[Species.from_str(sp).element] = [Species.from_str(sp).oxi_state]
    oxi_probs = {Species.from_str(sp): data / total_data[Species.from_str(sp).element] for sp, data in
                            all_data["occurrence"].items()}


    for index, row in df.iterrows():
        elA = pymatgen.core.periodic_table.get_el_sp(row["A"])
        elB = pymatgen.core.periodic_table.get_el_sp(row["B"])
        elX = pymatgen.core.periodic_table.get_el_sp(row["X"])
        name = row["A"] + row["B"] + row["X"]
        Prob = 0
        if elA.is_alkali and elX.is_halogen:
            Prob = 1#oxi_probs[Species.from_str(elB.symbol + "2+")]
        elif elA.is_alkali and elX.is_chalcogen:
            Prob = 1#(oxi_probs[Species.from_str(elB.symbol + "5+")] + oxi_probs[Species.from_str(elB.symbol + "6+")]) / 2
        elif elA.is_alkaline and elX.is_chalcogen:
            Prob = 1#oxi_probs[Species.from_str(elB.symbol + "4+")]
        elif elA.is_lanthanoid and elX.is_chalcogen:
            Prob = 1#(oxi_probs[Species.from_str(elB.symbol + "2+")] + oxi_probs[Species.from_str(elB.symbol + "3+")]) / 2
        elif elA.is_alkali and elX.symbol == "H":
            Prob = 1#oxi_probs[Species.from_str(elB.symbol + "2+")]
        else:
            continue

        if row["isperov"] and row["RaIX"] > 0 and row["RaXII"] > 0 and row["n"] == maxn[name] and Prob == 1 and row["n"] > 0:
            list_data.append(get_nn_input(row))
            list_target.append(maxn[name])

    norm_factor = np.max(np.sqrt(np.sum(np.array(list_data)**2, axis=1)))
    list_data = np.array(list_data) / norm_factor

    #kmeans = KMeans(n_clusters=30, max_iter=1000, random_state=0, n_init="auto").fit(list_data)

    nn = neural_network.MLPRegressor(hidden_layer_sizes=(100, ), verbose=True, tol=0.00000001, max_iter=1)

    #nn = linear_model.LinearRegression()
    nn.fit(list_data, list_target)

    #order = np.argsort(np.abs(nn.coef_))

    poly = PolynomialFeatures(2, interaction_only=True)
    poly.fit_transform([[1, 1, 1, 1, 1, 1, 1]])
    names = poly.get_feature_names_out()
    #for coef in order:
        #print(str(coef) + " " + str(names[coef]) + " " + str(nn.coef_[coef]))


    fig, ax = plt.subplots()

    x_vals_exp = []
    y_vals_exp = []
    c_exp = []

    x_vals_t = []
    y_vals_t = []
    c_t = []

    x_vals_non_perov = []
    y_vals_non_perov = []
    c_non_perov = []

    tot = 0
    num = 0

    for index, row in df.iterrows():
        if row["n"] >= 0:

            n = row["n"]
            if row["X"] == "H":
                row["Rc"] = 1.4
            if row["B"] == "Mn" and row["Box"] == 6:
                row["Rb"] = 0.255
            elA = pymatgen.core.periodic_table.get_el_sp(row["A"])
            elB = pymatgen.core.periodic_table.get_el_sp(row["B"])
            elX = pymatgen.core.periodic_table.get_el_sp(row["X"])
            electroA = elA.X
            electroB = elB.X
            electroX = elX.X

            if not elX.is_chalcogen:
                continue

            no_IX_XII = False

            if row["RaIX"] == 0.0 and row["RaXII"] == 0.0:
                row["RaIX"] = 10
                row["RaXII"] = -10
                no_IX_XII = True
                continue
            if row["RaIX"] == 0.0:
                row["RaIX"] = row["RaXII"] * 0.911
                no_IX_XII = True
                continue
            if row["RaXII"] == 0.0:
                row["RaXII"] = row["RaIX"] / 0.911
                no_IX_XII = True
                continue

            try:
                elA = Species(elA, row["Aox"])
                elB = Species(elB, row["Box"])
                elX = Species(elX, row["Cox"])

                row["RaIX"] = elA.get_shannon_radius("IX", radius_type="crystal")
                row["RaXII"] = elA.get_shannon_radius("XII", radius_type="crystal")
                row["Rb"] = elB.get_shannon_radius("VI", radius_type="crystal")
                row["Rc"] = elX.get_shannon_radius("VI", radius_type="crystal")

            except Exception:
                continue


            if n == 0:
                A_ave = (row["RaXII"])
                X_ave = (3 * electroX + electroB + electroA) / 5
                X_A_ave = electroA / (12 * electroX)
            elif n == 1:
                A_ave = (row["RaIX"])
                X_A_ave = electroA / (9 * electroX)
                X_ave = (4 * electroX + electroB + 2 * electroA) / 7
            elif n == 2:
                A_ave = ((2 * row["RaIX"] + row["RaXII"]) / 3)
                X_A_ave = (2 * electroA / (9 * electroX) + electroA / (12 * electroX)) / 3
                X_ave = (7 * electroX + 2 * electroB + 3 * electroA) / 12
            elif n == 3:
                A_ave = ((2 * row["RaIX"] + 2 * row["RaXII"]) / 4)
                X_A_ave = (2 * electroA / (9 * electroX) + 2 * electroA / (12 * electroX)) / 4
                X_ave = (10 * electroX + 3 * electroB + 4 * electroA) / 17
            elif n == 4:
                A_ave = ((2 * row["RaIX"] + 3 * row["RaXII"]) / 5)
                X_A_ave = (2 * electroA / (9 * electroX) + 3 * electroA / (12 * electroX)) / 5
                X_ave = (13 * electroX + 4 * electroB + 5 * electroA) / 22
            elif n == 5:
                A_ave = ((2 * row["RaIX"] + 4 * row["RaXII"]) / 6)
                X_A_ave = (2 * electroA / (9 * electroX) + 4 * electroA / (12 * electroX)) / 6
                X_ave = (16 * electroX + 5 * electroB + 6 * electroA) / 27
            elif n == 6:
                A_ave = ((2 * row["RaIX"] + 5 * row["RaXII"]) / 7)
                X_A_ave = (2 * electroA / (9 * electroX) + 5 * electroA / (12 * electroX)) / 7
                X_ave = (19 * electroX + 6 * electroB + 7 * electroA) / 32
            elif n == 7:
                A_ave = ((2 * row["RaIX"] + 6 * row["RaXII"]) / 8)
                X_A_ave = (2 * electroA / (9 * electroX) + 6 * electroA / (12 * electroX)) / 8
                X_ave = (22 * electroX + 7 * electroB + 8 * electroA) / 37
            if A_ave == 0:
                A_ave = 0.0001

            A_ave_gs = (A_ave + row["Rb"]) / (row["Rb"] + row["Rc"])

            goldschmidt_tolerance_factor = (row["RaXII"] + row["Rc"]) / (np.sqrt(2) * (row["Rb"] + row["Rc"]))

            layer_volume_ratio = (2 * row["RaIX"] ** 3 + 1 * row["Rb"] ** 3 + 4 * row["Rc"] ** 3) / (
                        row["RaXII"] ** 3 + row["Rb"] ** 3 + 3 * row["Rc"] ** 3)
            layer_volume_ratio = (row["RaIX"] ** 3 + row["Rc"] ** 3) / (
                    row["RaXII"] ** 3 + row["Rb"] ** 3 + 3 * row["Rc"] ** 3)

            a_to_c_ratio = (row["Rb"] + row["Rc"])/(2*row["Rb"] + 2*row["Rc"]+ 2*row["RaIX"]+0*row["RaXII"])

            #iondif = pymatgen.core.periodic_table.get_el_sp(row["B"]).ionization_energies[2] - pymatgen.core.periodic_table.get_el_sp(row["B"]).ionization_energies[1]

            t = nn.predict([(np.array(get_nn_input(row)) / norm_factor)])[0]
            t = row["RaXII"] - row["RaIX"]

            Prob = -3

            try:

                if elA.is_alkali and elX.is_halogen:
                    Prob = oxi_probs[Species.from_str(elB.symbol + "2+")]
                elif elA.is_alkali and elX.is_chalcogen:
                    Prob = (oxi_probs[Species.from_str(elB.symbol + "5+")]+oxi_probs[Species.from_str(elB.symbol + "6+")])/2
                elif elA.is_alkaline and elX.is_chalcogen:
                    Prob = oxi_probs[Species.from_str(elB.symbol + "4+")]
                elif elA.is_lanthanoid and elX.is_chalcogen:
                    Prob = (oxi_probs[Species.from_str(elB.symbol + "2+")]+oxi_probs[Species.from_str(elB.symbol + "3+")])/2
                elif elA.is_alkali and elX.symbol == "H":
                    Prob = oxi_probs[Species.from_str(elB.symbol + "2+")]
                else:
                    t = -1
            except Exception:
                t = t

            y_val = Prob
            # y_val = (2*row["Rb"] + 2*row["Rc"]+ 2*row["RaIX"]+0*row["RaXII"])
            # t = goldschmidt_tolerance_factor
            y_val = (n + index * 0.0001)

            #y_val = (row["RaIX"]+row["RaXII"])/row["Rc"]

            name = row["A"] + row["B"] + row["X"]

            if (not row["exp"]) and row["isperov"]:
                name = row["A"] + row["B"] + row["X"]
                if name != "":
                    x_vals_exp.append(t)
                    y_vals_exp.append(y_val)
                    if no_IX_XII or name not in neq0:
                        c_exp.append(1.5)
                    else:
                        c_exp.append(maxn[name])

                    if name in maxn.keys() and t > 0:
                        tot += (t - maxn[name]) ** 2
                        num += 1
            elif (row["exp"]) and (row["isperov"]):
                name = row["A"] + row["B"] + row["X"]
                if name != "":
                    x_vals_t.append(t)
                    y_vals_t.append(y_val)
                    if no_IX_XII:
                        c_t.append(1.5)
                    else:
                        c_t.append(maxn[name])
            elif (not row["exp"]) and (not row["isperov"]):
                name = row["A"] + row["B"] + row["X"]
                if name != "":
                    x_vals_non_perov.append(t)
                    y_vals_non_perov.append(y_val)
                    if no_IX_XII or name not in neq0 or maxn[name] == 0:
                        c_non_perov.append(3.5)
                    else:
                        c_non_perov.append(maxn[name])

            if (((not row["exp"]) and (row["isperov"] or show_non_perov) and row["n"] > -1) or (
                    show_theoretical and row["n"] > 0)) and labels:
                text = r"$" + row["A"] + "_" + str(int(n + 1)) + row["B"] + "_" + str(int(n)) + row["X"] + "_{" + str(
                    int(3 * n + 1)) + "}$ "
                if False:
                    ax.text(t, y_val, text, rotation=50, fontsize=10, alpha=0.25, color="red")
                else:
                    ax.text(t, y_val, text, rotation=50, fontsize=10, alpha=0.25, color="black")
                if row["A"] + row["B"] + row["X"] == "CeZnO" or row["A"] + row["B"] + row["X"] == "CeCuO":
                    print(text + " " + str(t))

    scatter = ax.scatter(x_vals_exp, y_vals_exp, c=c_exp, marker="o", cmap="Dark2", vmin=0, s=100, vmax=5,
                         label="Experimental structures")
    if show_theoretical:
        ax.scatter(x_vals_t, y_vals_t, c=c_t, marker="x", cmap="Dark2", vmin=0, s=100, vmax=5,
                   label="Theoretical structures")

    if show_non_perov:
        ax.scatter(x_vals_non_perov, y_vals_non_perov, c=c_non_perov, marker="_", cmap="Dark2", vmin=0, s=100, vmax=5,
                   label="Non Perovskite/Ruddlesden-Popper structures")

    #print("MSE: " + str(round(tot/num, 3)))

    plt.legend()
    ax.set_xlabel("Tolerance factor (T)")
    ax.set_ylabel("Homology (n number)")
    #ax.set_ylim([-0.1, 5.5])
    #ax.set_xlim([0.0, 1])
    ax.set_title(r"Tolerance Factor")
    ax.grid()

    fig.colorbar(scatter, ax=ax, label="Max Homology")

    plt.show()

    #plt.savefig("RP.jpg", format="jpg", dpi=5000)


def plot_3d():
    show_theoretical = False

    raw_df = pd.read_csv("Ruddlesden-Popper_all_sg_exp_structures_hsls_all.csv")
    data = pd.DataFrame(columns=["n", "Ra", "Rb", "Rc", "A/B", "A/X", "B/X"])
    alldata = pd.DataFrame(columns=["n", "Ra", "Rb", "Rc", "A/B", "A/X", "B/X"])

    df = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH"])
    for index, row in raw_df.iterrows():
        if row["Aox"] == 2.0 or True:
            df.loc[len(df)] = list(row)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    x_vals_exp = []
    y_vals_exp = []
    z_vals_exp = []
    c_exp = []

    x_vals_t = []
    y_vals_t = []
    z_vals_t = []
    c_t = []

    for index, row in df.iterrows():
        if row["n"] >= 0 and row["B"] != "Tc":
            n = row["n"]
            A = row["RaIX"]
            B = row["Rb"]
            X = row["Rc"]
            electroA = pymatgen.core.periodic_table.get_el_sp(row["A"]).X
            electroB = pymatgen.core.periodic_table.get_el_sp(row["B"]).X
            electroX = pymatgen.core.periodic_table.get_el_sp(row["X"]).X
            electroA = 3 * electroX - electroB - electroA

            goldschmidt_tolerance_factor = (row["RaIX"] + row["Rc"]) / (np.sqrt(2) * (row["Rb"] + row["Rc"]))
            tf = (3 * np.sqrt(2) * row["Rc"] + 2 * np.sqrt(6) * ((row["RaIX"] + row["Rc"]))) / (
                    9 * ((row["Rb"] + row["Rc"])))

            if n == 0:
                A_ave = (row["RaXII"])
            elif n == 1:
                A_ave = (row["RaIX"])
            elif n == 2:
                A_ave = ((2 * row["RaIX"] + row["RaXII"]) / 3)
            elif n == 3:
                A_ave = ((2 * row["RaIX"] + 2 * row["RaXII"]) / 4)
            elif n == 4:
                A_ave = ((2 * row["RaIX"] + 3 * row["RaXII"]) / 5)
            elif n == 5:
                A_ave = ((2 * row["RaIX"] + 4 * row["RaXII"]) / 6)

            if row["RaIX"] == 0.0:
                row["RaIX"] = -10
            if row["RaXII"] == 0.0:
                row["RaXII"] = 10

            if A_ave == 0:
                A_ave = 0.0001
            A_ave_gs = (A_ave + row["Rc"]) / (row["Rb"] + row["Rc"])

            diff_shannon_tolerance_factor = ((row["RaIX"] - row["RaXII"]) + row["Rc"]) / (
                    row["Rb"] + row["Rc"])
            ratio_tolerance_factor = 1
            # t = (row["Rb"] + row["RaXII"])/2 - ((row["RaIX"] + row["Rb"]) /2 +  row["Rc"])

            t1 = row["RaIX"] / row["Rc"]
            t2 = goldschmidt_tolerance_factor

            if not row["exp"]:
                c_exp.append(row["EAH"])
                x_vals_exp.append(t1)
                y_vals_exp.append(n)
                z_vals_exp.append(t2)
            else:
                c_t.append(row["EAH"])
                x_vals_t.append(t1)
                y_vals_t.append(n)
                z_vals_t.append(t2)

            if not row["exp"] or (show_theoretical and row["n"] > 1) and t1 > -4:
                text = r"$" + row["A"] + "_" + str(int(n + 1)) + row["B"] + "_" + str(int(n)) + row["X"] + "_{" + str(
                    int(3 * n + 1)) + "}$"
                # print(text)
                ax.text(t1, n, t2, text, rotation=50, fontsize=10)

    scatter = ax.scatter(x_vals_exp, y_vals_exp, z_vals_exp, c=c_exp, marker="o", cmap="winter", vmin=0, s=100,
                         vmax=0.25,
                         label="Experimental structures")

    plt.legend()
    ax.set_xlabel("Tolerance factor (T)")
    ax.set_ylabel("Homology (n number)")
    ax.set_zlabel("Goldschmidt")
    ax.set_title(r"Tolerance Factor for $A_{n+1}B_nO_{3n+1}$ Ruddlesden-Popper series")
    ax.grid()
    ax.set_xlim([0.0, 4.0])
    ax.set_xlim([0.0, 4.0])
    ax.set_xlim([0.0, 4.0])

    plt.show()


if __name__ == "__main__":
    #plot()
    search_all()
    # make_raw_csv_parallel()
    1 / 0
    rawdf = pd.read_csv("Ruddlesden-Popper_expanded_data.csv")
    rawdf = rawdf.drop(columns=["Unnamed: 0"])

    df = pandas.DataFrame(
        columns=["struct", "n", "A", "B", "X", "Aox", "Box", "Cox", "RaIX", "RaXII", "Rb", "Rc", "exp", "EAH",
                 "isperov"])

    start = time.time()
    num_structs = 0
    for index, row in rawdf.iterrows():
        try:
            if row["n"] > 0:
                mat = pymatgen.core.Structure.from_dict(json.loads(row["struct"]))
                swapped = False

                for el in mat.composition.reduced_composition.as_dict().keys():
                    if mat.composition.reduced_composition.as_dict()[el] == row["n"] and row["A"] == el:
                        row_list = list(row)
                        tmp = row_list[2]
                        row_list[2] = row_list[3]
                        row_list[3] = tmp
                        data = get_df_column_shannon(mat, row["exp"], row["EAH"])
                        if data[1] != row_list[2]:
                            print("Issue with " + mat.composition.reduced_composition.to_pretty_string())
                            1 / 0
                        df.loc[num_structs] = [json.dumps(mat.as_dict())] + data + [row["isperov"]]
                        num_structs += 1
                        swapped = True
                        #print(mat.composition.to_pretty_string())
                if not swapped:
                    df.loc[num_structs] = row
                    num_structs += 1
            else:
                df.loc[num_structs] = row
                num_structs += 1

            if index % 25 == 0 and index > 0 and False:
                current = time.time()
                elapsed = start - current
                projection = elapsed * len(rawdf) / index
                print("\npercent complete: " + str(round(100 * index / len(rawdf), 2)) + "%, elapsed time: "
                      + str(round(elapsed / 60, 2)) + "Minutes. Projected time remained: " + str(
                    round(projection / 60, 2)) + "Minutes.\n")
        except Exception:
            print(row)
            1 / 0
    df.to_csv("Ruddlesden-Popper_expanded_data_cleaned.csv", index=False)
