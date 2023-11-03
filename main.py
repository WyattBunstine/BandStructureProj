import sys
import os
from warnings import *
import json
import pymatgen
from pymatgen.core import *
from mp_api.client import MPRester
import Utils.KPath
import Utils.MPSearch
import subprocess
import time
import shutil
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROUNDING = 4
numBandPoints = 400
BtoA = 0.52917721090  # 1 bohr radius is 0.52 angstroms
HtoEv = 27.2114 # 1 Hartree is 27 Electron volts

def validate_config(config):
    """This function validates the config file used to run elk.

    :param config: a dictionary containing configuration information
    :return: True if config contains necessary keys, false otherwise
    """

    # default required configuration options
    recs = ["MatType", "MatLoc", "MatID", "tasks", "kpoints", "elkparams", "ROUNDING", "numBandPoints", "remote",
            "monitor"]
    # check that each is present, if one is missing, return false
    for rec in recs:
        if rec not in config.keys():
            return False
    # otherwise return true
    return True


def get_mp_mat(MPID: str) -> pymatgen.core.structure.Structure:
    """This function retreives a material from materials project from a materials project ID

    :param MPID: the string of the mpid, formatted like "mp-xxxxx"
    :return: the pymatgen structure of the material
    """
    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    with MPRester(API_KEY) as mpr:
        return mpr.get_structure_by_material_id(MPID)


def gen_spacegroupin(path, mat: Structure):
    """Generate the GEOMETRY.OUT file for elk. First generates spacegroup.in, then calls the spacegroup utility

    :param path: the path where the spacegroup.in file should be generated
    :param mat: the material for which to generate the file
    :return: 0 on success, -1 otherwise
    """
    file = open(path + "spacegroup.in", 'w')
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

    return 0


def gen_elk_in(path, mat, config):
    """Generate the elk.in file to run DFT with

    :param path: the path with the GEOMETRY.OUT file where elk.in will be generated
    :param mat: the material for which to run elk
    :param config: the configuration dictionary with elk parameters
    :return: 0 on success, -1 otherwise
    """

    # open the file
    file = open(path + "elk.in", 'w+')
    # write basic params from config
    for param in config["elkparams"].keys():
        file.write(param + "\n  ")
        file.write(str(config["elkparams"][param]) + "\n\n")

    # add appropriate tasks
    file.write("tasks\n")
    if "GS" in config["tasks"]: file.write("  0\n")
    if "BS" in config["tasks"]: file.write("  20\n")
    file.write("\n")

    # if band structure is to be generated, add plotting parameters and add k points
    if "BS" in config["tasks"]:
        file.write("plot1d\n")
        # adding k points to plot
        kpoints = Utils.KPath.get_kpoints_SC(mat)
        file.write("  " + str(len(kpoints)) + " " + str(config["numBandPoints"]) + "\n")
        for point in kpoints:
            coords = "  " + str(point[1][0]) + " " + str(point[1][1]) + " " + str(point[1][2]) + "\n"
            file.write(coords)
        file.write("\n\n")
    # add path to the elk files
    if config["remote"]:
        file.write("sppath\n  '/home/wbunsti1/elk/elk-8.8.26/species/'\n\n")
    else:
        file.write("sppath\n  '/home/wyatt/elk-8.8.26/species/'\n\n")
    # add the crystal geometry from GEOMETRY.OUT file

    geometry = False
    if geometry:
        geometryfile = open(path + "GEOMETRY.OUT", 'r')
        for line in geometryfile:
            file.write(line)
        geometryfile.close()
    else:
        struct = mat.get_primitive_structure()
        file.write("avec\n")
        for vec in struct.lattice.matrix:
            file.write("  ")
            for param in vec:
                a = round(param / BtoA, 10)
                file.write(f'{a:.10f}'+"  ")
            file.write("\n")

        species = {}
        for site in struct.sites:
            for specie in site.species:
                if str(specie) in species.keys():
                    species[str(specie)].append(struct.lattice.get_fractional_coords(site.coords))
                else:
                    species[str(specie)] = [struct.lattice.get_fractional_coords(site.coords)]

        file.write("\natoms\n")
        file.write("  " + str(len(species.keys()))+"\n")
        for specie in species.keys():
            file.write("'" + specie + ".in'\n  " + str(len(species[specie]))+"\n")
            for site in species[specie]:
                file.write("  ")
                for coord in site:
                    a = round(coord, ROUNDING)
                    file.write(f'{a:.10f}'+"  ")
                file.write("\n")

    file.close()
    return 0


def run_elk(path, config):
    """Function that actually runs elk, can run remotely on rockfish or locally

    :param path: the path to elk.in
    :param config: the config file for the DFT run
    :return: 0 on success, -1 on failure or timeout
    """
    # block for remote running
    if config["remote"]:
        # create directory for this DFT run
        out = subprocess.run(
            ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "elk/elk-8.8.26/BSPOperDir;", "mkdir", config["MatID"]])
        # if the directory exists, clean it
        if not out.returncode == 0:
            # if the old file should be saved
            if not config["overwrite"]:
                subprocess.run(
                    ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "elk/elk-8.8.26/BSPOperDir;", "zip", "-r", "old/" +
                     config["MatID"] + ".zip", config["MatID"]])
            # remove the old folder
            subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "elk/elk-8.8.26/BSPOperDir;", "rm", "-r",
                 config["MatID"] + ";", "rm", "-r", config["MatID"] + ";"])
            # create a fresh one
            subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "elk/elk-8.8.26/BSPOperDir;", "mkdir",
                 config["MatID"]])
        # copy elk.in and elk.slurm
        subprocess.run(["scp", path + "elk.in",
                        "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + config["MatID"]
                        + "/elk.in"])
        subprocess.run(["scp", "elk.slurm",
                        "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + config[
                            "MatID"] + "/elk.slurm"])
        # run the process
        subprocess.run(
            ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "elk/elk-8.8.26/BSPOperDir/" + config["MatID"] + ";",
             "sbatch elk.slurm"])
        # block to monitor progress
        if config["monitor"]:
            # track time it takes
            start = time.time()
            cont = True
            failures = 0
            # while it is not converged, out of loops, or fails
            while cont:
                # wait for 30 seconds before updating
                time.sleep(30)
                # get rid of old INFO.OUT
                if os.path.exists(path + "INFO.OUT"):
                    os.remove(path + "INFO.OUT")
                # copy INFO.OUT file from remote
                subprocess.run(
                    ["scp", "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + config[
                        "MatID"] + "/INFO.OUT", path + "INFO.OUT"])
                # if it was copied successfully
                if os.path.exists(path + "INFO.OUT"):
                    failures = 0
                    # open info file
                    with open(path + "INFO.OUT", 'rb') as f:
                        try:  # catch OSError in case convergence info can't be found
                            f.seek(-2, os.SEEK_END)
                            seek = True
                            while seek:
                                text = f.read(23)
                                # either find the end of looping or find the convergence info
                                if text == b"RMS change in Kohn-Sham":
                                    seek = False
                                elif text == b"Elk version 8.8.26 stop":
                                    seek = False
                                    cont = False
                                else:
                                    f.seek(-24, os.SEEK_CUR)
                            # if the convergence info was found, print an update
                            if cont:
                                f.seek(-23, os.SEEK_CUR)
                                conv_crit = f.readline().decode() + f.readline().decode()
                                seek = True
                                while seek:
                                    text = f.read(15)
                                    if text == b"Loop number :  ":
                                        seek = False
                                    else:
                                        f.seek(-16, os.SEEK_CUR)
                                print("Elapsed time: " + str(time.time() - start)[0:5] + " seconds. Loop Number: "
                                      + f.read(4).decode())
                                print(conv_crit, end="")
                        # the block to catch an error finding convergence info
                        except OSError:
                            print("Failed finding INFO.OUT convergence information.")
                            f.seek(0)
                # if it fails coping INFO.OUT, allow for breaking
                else:
                    failures += 1
                    print("Failed coping INFO.OUT " + str(failures) + " times.")
                    if failures == 5:
                        wait = input("Continue waiting?(y/n)")
                        if wait == "y" or wait == "Y":
                            failures = 0
                        else:
                            return -1

            # copy the remote files upon completion
            print("ELK finished. Coping remote files")
            subprocess.run(
                ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + config[
                    "MatID"] + "/*", path])
    else:
        os.system("cd " + path + " && . /opt/intel/oneapi/setvars.sh && /home/wyatt/Downloads/elk/elk-8.8.26/src/elk")

    return 0


def from_config(config):
    """This functions runs DFT on a material using elk with parameters from a given config file.
    Config specifies everything about the material and the DFT params

    :param config: the config file location
    :return: 0 on success, -1 otherwise
    """
    save_cfg = False
    # check if config exists
    if not type(config) == dict:
        if not os.path.exists(config):
            warn("config file, " + str(config) + ", does not exist")
            return -1
        with open(config) as f:
            config = json.load(f)
    else:
        save_cfg = True
    if not validate_config(config):
        warn("config invalid")
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

    if os.path.isdir("data/" + config["MatID"]):
        if not config["overwrite"]:
            subprocess.run(
                ["zip", "-r", "data/old/" + config["MatID"] + ".zip", "data/" + config["MatID"]])
        shutil.rmtree("data/" + config["MatID"])
    os.mkdir("data/" + config["MatID"])
    if save_cfg:
        with open("data/" + config["MatID"] + '/' + config["MatID"] + '_config.json', 'w') as f:
            json.dump(config, f)

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


def download_remote(matid: str):
    subprocess.run(
        ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + matid + "/*",
         "data/"+matid])
    return os.path.exists("data/"+matid+"/BAND.OUT")


def main():
    run = False
    elements = ["Ni", "Cr"]
    if run:
        cfgs = Utils.MPSearch.get_configs(elements)
        if os.path.exists("running_configs"+str(elements)+".txt"):
            warn("running_configs"+str(elements)+".txt exists already, make sure to download and delete file")
            return 0
        running = open("running_configs"+str(elements)+".txt", "w+")
        completed = open("completed_runs.txt", "r")
        comp_lines = list(completed.readlines())
        completed.close()
        for config in cfgs:
            if config["MatID"] not in comp_lines:
                running.write(config["MatID"]+"\n")
                from_config(config)
        running.close()
    else:
        completed = open("completed_runs.txt", "a")
        downloads = open("running_configs" + str(elements) + ".txt", "r")
        for line in downloads.readlines():
            if download_remote(line[:-1]):
                completed.write(line+"\n")
        completed.close()
        downloads.close()


if __name__ == "__main__":
    main()
