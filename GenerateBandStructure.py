import shutil

from mp_api.client import MPRester
import json
import os
import pandas as pd
from io import StringIO
import plot
import numpy as np
import sys
import pymatgen
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import *
import Utils.KPath

ROUNDING = 4
numBandPoints = 400
BtoA = 0.52917721090  # 1 bohr radius is 0.52 angstroms


def generateSpacegroupin(path, mat):
    file = open(path + "spacegroup.in", 'w+')

    file.write("'" + mat["symmetry"] + "'\n")
    structure = pymatgen.core.structure.Structure.from_dict(mat["structure"])
    latparams = list(structure.lattice.parameters)
    for param in latparams:
        latparams.append(round(latparams[0], ROUNDING))
        latparams.remove(latparams[0])
    file.write(str(latparams[0] / BtoA) + " " + str(latparams[1] / BtoA) + " " + str(latparams[2] / BtoA) + "\n")
    file.write(str(latparams[3]) + " " + str(latparams[4]) + " " + str(latparams[5]) + "\n")
    file.write("1 1 1\n.false.\n")
    file.write(str(mat["nelements"]) + "\n")

    elements = dict()
    for el in structure.sites:
        if not el.species in elements.keys():
            elements[el.species] = 1
        else:
            elements[el.species] += 1

    for el in elements.keys():
        file.write("'" + str(el.elements[0].symbol) + "' '" + str(el.elements[0].symbol) + ".in'\n")
        file.write(str(elements[el]) + "\n")
        for site in structure.sites:
            if site.species == el:
                file.write(str(round(float(site.a), ROUNDING)) + " " + str(round(float(site.b), ROUNDING)) + " "
                           + str(round(float(site.c), ROUNDING)) + "\n")
    file.close()
    os.system("cd " + path + " && spacegroup")


def generateBandStructure(path, mat, tasks=["GSC", "genBandStructure"]):
    if "GSC" in tasks:
        if os.path.exists(path + "/BAND.OUT"):
            os.remove(path + "/BAND.OUT")
        file = open(path + "elk.in", 'w+')

        basicfile = open("elkinbasic.txt", 'r')
        for line in basicfile:
            file.write(line)
        basicfile.close()

        # adding k points to plot
        file.write("  " + str(len(mat["kpoints"])) + " " + str(numBandPoints) + "\n")
        for point in mat["kpoints"]:
            coords = "  " + str(point[1][0]) + " " + str(point[1][1]) + " " + str(point[1][2]) + "\n"
            file.write(coords)
        file.write("\n\n")

        geometryfile = open(path + "GEOMETRY.OUT", 'r')
        for line in geometryfile:
            file.write(line)
        geometryfile.close()

        file.close()
        os.system("cd " + path + " && . /opt/intel/oneapi/setvars.sh && /home/wyatt/Downloads/elk/elk-8.8.26/src/elk")

    if "genBandStructure" in tasks and 1 == 0:
        if os.path.exists(path + "/BAND.OUT"):
            os.remove(path + "/BAND.OUT")
        file = open(path + "elk.in", 'w+')

        basicfile = open("elkinbasic.txt", 'r')
        for line in basicfile:
            file.write(line)
        basicfile.close()

        '''removed for now
        #adding k points to plot
        file.write("  " + str(len(mat["kpoints"])) + " " + str(numBandPoints)+"\n")
        for point in mat["kpoints"]:
            coords = "  " + str(point[1][0]) + " " + str(point[1][1]) + " " + str(point[1][2]) + "\n"
            file.write(coords)
        file.write("\n\n")
        '''
        # add the ground state calculation task:
        file.write("tasks\n  20\n\n")

        # add the plot block
        file.write("plot1d\n")
        file.write("  " + str(len(mat["kpoints"])) + " " + str(numBandPoints) + "\n")
        for point in mat["kpoints"]:
            coords = "  " + str(point[1][0]) + " " + str(point[1][1]) + " " + str(point[1][2]) + "\n"
            file.write(coords)
        file.write("\n\n")

        geometryfile = open(path + "GEOMETRY.OUT", 'r')
        for line in geometryfile:
            file.write(line)
        geometryfile.close()

        file.close()
        os.system("cd " + path + " && . /opt/intel/oneapi/setvars.sh && /home/wyatt/Downloads/elk/elk-8.8.26/src/elk")


def GenerateMatJson(path, MPID, option=["mp"]):
    if os.path.exists(path):
        shutil.rmtree(path)
    if not os.path.isdir("data"):
        os.mkdir("data")
    if not os.path.isdir("data/" + MPID):
        os.mkdir("data/" + MPID)

    mat = dict()
    tmp = dict()
    if option[0] == "mp":
        API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
        # inputs
        with MPRester(API_KEY) as mpr:
            data = mpr.summary.search(material_ids=[MPID])
            # save the relevent data to a dict
            tmp = dict(data[0])
            requirements = ["nelements", "formula_pretty"]
            for key in tmp.keys():
                if key in requirements:
                    mat[key] = tmp[key]
                if key == "symmetry":
                    mat[key] = tmp[key].symbol
                if key == "structure":
                    sga = SpacegroupAnalyzer(tmp[key])
                    data = sga.get_refined_structure()
                    mat[key] = data.as_dict()

    else:
        structure = pymatgen.core.structure.Structure.from_file(option[0]).get_reduced_structure()
        mat["structure"] = structure.as_dict()
        mat["symmetry"] = structure.get_space_group_info()[0]
        mat["nelements"] = len(np.unique(structure.species))
        mat["formula_pretty"] = structure.composition.reduced_formula
    structure = pymatgen.core.structure.Structure.from_dict(mat["structure"])

    #mat["kpoints"] = Utils.KPath.get_kpoints(structure)

    #Another old way, correct for SC band structures
    kpoints = KPathSetyawanCurtarolo(structure).kpath
    points = []
    for point in kpoints['path'][0]:
        points.append([point, list(kpoints['kpoints'][point])])
    mat["kpoints"] = points

    '''old way of generating k points probably not correct
    kpoints = KPathSetyawanCurtarolo(structure).get_kpoints()
    points = []
    for point in range(len(kpoints[0])):
        if not kpoints[1][point] == '':
            coords = []
            for coord in kpoints[0][point]:
                coords.append(round(coord, ROUNDING))
            if len(points) == 0 or not points[-1][0] == kpoints[1][point]:
                latticeCoords = []
                for coord in structure.lattice.reciprocal_lattice.get_fractional_coords(coords):
                    latticeCoords.append(round(coord, ROUNDING))
                points.append([kpoints[1][point],latticeCoords])
    mat["kpoints"] = points
    '''
    file = open(path + MPID + ".json", 'w+')
    json.dump(mat, file, skipkeys=True)
    file.close()


def oldmain():
    ID = ''
    file = ''
    stop = False
    op = ["mp"]
    if "--MP" in sys.argv:
        if len(sys.argv) > sys.argv.index("--MP") + 1:
            ID = sys.argv[sys.argv.index("--MP") + 1]
            if ID[0:1] == "--":
                print("MPID not properly specified")
                stop = True
        else:
            print("MPID not properly specified")
            stop = True
    elif "--CIF" in sys.argv:
        if len(sys.argv) > sys.argv.index("--CIF") + 2:
            ID = sys.argv[sys.argv.index("--CIF") + 1]
            file = sys.argv[sys.argv.index("--CIF") + 2]
            op = [file]
            if ID[0:1] == "--" or not file[-4:] == ".cif":
                print("CIF not properly specified, file or id incorrect")
                stop = True
        else:
            print("CIF ID and file not specified")
            stop = True
    else:
        print("Material not properly specified, use --MP or --CIF")
        stop = True

    if len(sys.argv) < 2:
        print("no commands passed, options are:\n  --genBands\n  --plotBands\n  --genMatJson")
    elif not stop:
        path = "data/" + ID + "/"
        if "--genMatJson" in sys.argv:
            GenerateMatJson(path, ID, option=op)
        if not os.path.isfile(path + ID + ".json"):
            print("material JSON file not present, use --genMatJson to generate")
        else:
            file = open(path + ID + ".json", 'r')
            data = json.load(file)
            file.close()
            if "--genBands" in sys.argv:
                generateSpacegroupin(path, data)
                generateBandStructure(path, data)
            if "--plotBands" in sys.argv:
                if not os.path.isfile(path + "BAND.OUT"):
                    print("Bandlines not generated, use --genBands to generate before plotting")
                else:
                    plot.plot(path, data)
            if "--plotDOS" in sys.argv:
                if not os.path.isfile(path + "TDOS.OUT"):
                    print("Bandlines not generated, use --genBands to generate before plotting")
                else:
                    plot.plotDOS(path, data)

