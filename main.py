import sys
import os
from warnings import *
import json
import pymatgen
from mp_api.client import MPRester
import Utils.KPath
import subprocess
import time
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROUNDING = 4
numBandPoints = 400
BtoA = 0.52917721090  # 1 bohr radius is 0.52 angstroms


def validate_config(config):
    recs = ["MatType", "MatLoc", "MatID", "tasks", "kpoints", "elkparams", "ROUNDING", "numBandPoints", "remote",
            "monitor"]
    for rec in recs:
        if rec not in config.keys():
            return False
    return True


def get_mp_mat(MPID):
    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    # inputs
    with MPRester(API_KEY) as mpr:
        # save the relevent data to a dict
        return mpr.get_structure_by_material_id((MPID))


def gen_spacegroupin(path, mat):
    file = open(path + "/spacegroup.in", 'w')
    file.write("'" + mat.get_space_group_info()[0] + "'\n")
    latparams = list(mat.lattice.parameters)
    for param in latparams:
        latparams.append(round(latparams[0], ROUNDING))
        latparams.remove(latparams[0])
    file.write(str(latparams[0] / BtoA) + " " + str(latparams[1] / BtoA) + " " + str(latparams[2] / BtoA) + "\n")
    file.write(str(latparams[3]) + " " + str(latparams[4]) + " " + str(latparams[5]) + "\n")
    file.write("1 1 1\n.false.\n")

    elements = dict()
    for el in mat.sites:
        if not el.species in elements.keys():
            elements[el.species] = 1
        else:
            elements[el.species] += 1

    file.write(str(len(elements.keys())) + "\n")

    for el in elements.keys():
        file.write("'" + str(el.elements[0].symbol) + "' '" + str(el.elements[0].symbol) + ".in'\n")
        file.write(str(elements[el]) + "\n")
        for site in mat.sites:
            if site.species == el:
                file.write(str(round(float(site.a), ROUNDING)) + " " + str(round(float(site.b), ROUNDING)) + " "
                           + str(round(float(site.c), ROUNDING)) + "\n")
    file.close()
    os.system("cd " + path + " && spacegroup")


def gen_elk_in(path, mat, config):
    file = open(path + "elk.in", 'w+')

    for param in config["elkparams"].keys():
        file.write(param + "\n  ")
        file.write(str(config["elkparams"][param]) + "\n\n")

    if "BS" in config["tasks"]:
        file.write("plot1d\n")
        # adding k points to plot
        kpoints = Utils.KPath.get_kpoints_SC(mat)
        file.write("  " + str(len(kpoints)) + " " + str(config["numBandPoints"]) + "\n")
        for point in kpoints:
            coords = "  " + str(point[1][0]) + " " + str(point[1][1]) + " " + str(point[1][2]) + "\n"
            file.write(coords)
        file.write("\n\n")

    geometryfile = open(path + "GEOMETRY.OUT", 'r')
    for line in geometryfile:
        file.write(line)
    geometryfile.close()

    file.close()


def run_elk(path, config):
    if config["remote"]:
        subprocess.run(
            ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "elk/elk-8.8.26/BSPOperDir;", "mkdir", config["MatID"]])
        subprocess.run(["scp", path + "elk.in",
                        "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + config["MatID"]
                        + "/elk.in"])
        subprocess.run(["scp", "elk.slurm",
                        "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + config[
                            "MatID"] + "/elk.slurm"])
        subprocess.run(
            ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "elk/elk-8.8.26/BSPOperDir/" + config["MatID"] + ";",
             "sbatch elk.slurm"])
        if config["monitor"]:
            cont = True
            failures = 0
            while cont:
                time.sleep(30)
                subprocess.run(
                    ["scp", "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + config[
                        "MatID"] + "/INFO.OUT", path + "INFO.OUT"])
                if os.path.exists(path + "INFO.OUT"):
                    failures = 0
                    with open(path + "INFO.OUT", 'rb') as f:
                        try:  # catch OSError in case of a one line file
                            f.seek(-2, os.SEEK_END)
                            while f.read(1) != b'\n':
                                f.seek(-2, os.SEEK_CUR)
                        except OSError:
                            f.seek(0)
                        last_line = f.readline().decode()
                    print("found")
                    return 0
                else:
                    failures += 1
                    print("failed coping INFO.OUT " + str(failures) + " times.")
                    if failures == 5:
                        print("Cannot copy INFO.OUT, aborting")
                        return -1
    else:
        os.system("cd " + path + " && . /opt/intel/oneapi/setvars.sh && /home/wyatt/Downloads/elk/elk-8.8.26/src/elk")

    return 0


def main():
    # check if the config file is given
    if len(sys.argv) < 2:
        warn("no args given")
        return -1
    # check if config exists
    config = sys.argv[1]
    if not os.path.exists(config):
        warn("config file, " + str(config) + ", does not exist")
        return -1
    with open(config) as f:
        config = json.load(f)
    if not validate_config(config):
        warn("config file invalid")
        return -1
    # get the material data, from materials proj or cif
    mat = 0
    if config["MatType"] == "MP":
        mat = get_mp_mat(config["MatLoc"])
    elif config["MatType"] == "cif":
        if not os.path.exists(config["MatLoc"]):
            warn("specified cif does not exist")
            return -1
        mat = pymatgen.core.structure.Structure.from_file(config["MatLoc"])
    else:
        warn("MatType not properly specified")
        return -1
    # create the data directory
    if not os.path.isdir("data/" + config["MatID"]):
        os.mkdir("data/" + config["MatID"])

    # ELK functions
    if "GS" in config["tasks"]:
        gen_spacegroupin("data/" + config["MatID"] + "/", mat)
        gen_elk_in("data/" + config["MatID"] + "/", mat, config)
        match run_elk("data/" + config["MatID"] + "/", config):
            case 0:
                print("Elk run success")
            case -1:
                print("Elk monitoring timeout")
                return -1
            case _:
                print("Unknown return from elk run")
                return -1

    return 0


if __name__ == "__main__":
    main()
