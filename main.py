import sys
import os
from warnings import *
import json

import numpy as np
import pymatgen
from pymatgen.core import *
from mp_api.client import MPRester
import Utils.KPath
import Utils.MPSearch
import subprocess
import time
import shutil
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import plot

ROUNDING = 4
numBandPoints = 400
BtoA = 0.52917721090  # 1 bohr radius is 0.52 angstroms
HtoEv = 27.2114  # 1 Hartree is 27 Electron volts
API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"


def validate_config(config):
    """This function validates the config file used to run elk.

    :param config: a dictionary containing configuration information
    :return: True if config contains necessary keys, false otherwise
    """

    # default required configuration options
    recs = ["MatType", "MatLoc", "MatID", "tasks", "kpoints", "elkparams", "ROUNDING", "numBandPoints", "remote",
            "monitor", "SGU"]
    # check that each is present, if one is missing, return false
    for rec in recs:
        if rec not in config.keys():
            warn(rec + " not found in config")
            return False
    if "PHONON" in config["tasks"]:
        if len(config["tasks"]) > 1:
            warn("Phonon calculations should be performed on their own, additional tasks present")
            return False
        if "phonon_workers" not in config.keys():
            warn("phonon_workers missing for phonon calcuations")
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


def gen_elk_in(path, mat: pymatgen.core.structure.Structure, config):
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
    if "SO" in config["tasks"]: file.write("  2\n")
    if "DOS" in config["tasks"]: file.write("  10\n")
    if "BS" in config["tasks"]: file.write("  20\n")
    if "OBS" in config["tasks"]: file.write("  22\n")
    if "SBS" in config["tasks"]: file.write("  23\n")
    if "PHONON" in config["tasks"]: file.write("  200\n")
    if "WANNIER" in config["tasks"]: file.write("  550\n")
    for task in config["tasks"]:
        if type(task) == int:
            file.write("  " + str(task) + "\n")
    file.write("\n")

    # if band structure is to be generated, add plotting parameters and add k points
    if "BS" in config["tasks"] or "OBS" in config["tasks"] or "SBS" in config["tasks"] or "PHONON" in config["tasks"]:
        file.write("plot1d\n")
        # adding k points to plot
        if config["kpoints"] == "MP":
            kpoints = Utils.KPath.get_kpoints_SC(mat)
        elif config["kpoints"] == "McQueen":
            kpoints = Utils.KPath.get_kpoints_McQueen(mat)
        else:
            kpoints = config["kpoints"]
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

    # add the crystal geometry from GEOMETRY.OUT file, this only works for some crystals
    if config["SGU"]:
        geometryfile = open(path + "GEOMETRY.OUT", 'r')
        for line in geometryfile:
            file.write(line)
        geometryfile.close()
    # adding geometry directly, this is computed from lattice vectors and Wyckoff positions
    else:
        # add the lattice vectors
        struct = mat
        if config["elkparams"]["primcell"] == ".true.":
            struct = mat.get_primitive_structure()
        file.write("avec\n")
        for vec in struct.lattice.matrix:
            file.write("  ")
            for param in vec:
                a = round(param / BtoA, 10)
                file.write(f'{a:.10f}' + "  ")
            file.write("\n")
        # get the atoms and their positions
        species = {}
        for site in struct.sites:
            for specie in site.species:
                try:
                    if str(specie.element) in species.keys():
                        species[str(specie.element)].append(struct.lattice.get_fractional_coords(site.coords))
                    else:
                        species[str(specie.element)] = [struct.lattice.get_fractional_coords(site.coords)]
                except AttributeError:
                    if str(specie) in species.keys():
                        species[str(specie)].append(struct.lattice.get_fractional_coords(site.coords))
                    else:
                        species[str(specie)] = [struct.lattice.get_fractional_coords(site.coords)]
        # print positions to elk.in file
        file.write("\natoms\n")
        file.write("  " + str(len(species.keys())) + "\n")
        for specie in sorted(species.keys()):
            file.write("'" + specie + ".in'\n  " + str(len(species[specie])) + "\n")
            for site in species[specie]:
                file.write("  ")
                for coord in site:
                    a = round(coord, ROUNDING)
                    file.write(f'{a:.10f}' + "  ")
                file.write("\n")

    file.close()
    return 0


def gen_slurm(mat: pymatgen.core.structure.Structure, config):
    if config["slurmparams"]["job-name"] == 0:
        config["slurmparams"]["job-name"] = mat.formula
    prefix = ("#!/bin/bash -l\n#SBATCH --job-name=" + str(config["MatID"]).replace(" ", "")
              + "\n#SBATCH --time=" + config["slurmparams"]["time"] + "\n")
    suffix = '''#SBATCH --ntasks-per-node=48
#SBATCH --nodes=1
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=wbunsti1@jhu.edu
# MARCC defaults to gcc, must load intel module
module load intel/2022.2
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
export OMP_STACKSIZE=25G
ulimit -Ss unlimited
mpirun /home/wbunsti1/elk/elk-8.8.26/src/elk >& elk.log\n'''

    file = open(config["MatLoc"] + "elk.slurm", 'w+')
    file.write(prefix + suffix)
    file.close()


def run_elk(config):
    """Function that actually runs elk, can run remotely on rockfish or locally

    :param path: the path to elk.in
    :param config: the config file for the DFT run
    :return: 0 on success, -1 on failure or timeout
    """
    # block for remote running
    if config["remote"]:
        # create directory for this DFT run
        out = subprocess.run(
            ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "mkdir",
             config["MatLoc"] + ";", "exit"])
        # if the directory exists, clean it
        if not out.returncode == 0:
            # if the old file should be saved
            if not config["overwrite"]:
                subprocess.run(
                    ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "zip", "-r", ".old/" +
                     config["MatID"] + ".zip", config["MatLoc"] + ";", "exit"])
            # remove the old folder
            subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "rm", "-r",
                 config["MatLoc"] + ";", "rm", "-r", config["MatLoc"] + ";", "exit"])
            # create a fresh one
            subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "mkdir",
                 config["MatLoc"] + ";", "exit"])
        # copy elk.in and elk.slurm
        subprocess.run(["scp", config["MatLoc"] + "elk.in",
                        "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + config["MatLoc"]
                        + "elk.in"])
        subprocess.run(["scp", config["MatLoc"] + "elk.slurm",
                        "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + config[
                            "MatLoc"] + "elk.slurm"])
        # run the process
        subprocess.run(
            ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/" + config["MatLoc"] + ";",
             "sbatch elk.slurm; exit"])
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
                if os.path.exists(config["MatLoc"] + "INFO.OUT"):
                    os.remove(config["MatLoc"] + "INFO.OUT")
                # copy INFO.OUT file from remote
                subprocess.run(
                    ["scp", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + config[
                        "MatLoc"] + "INFO.OUT", config["MatLoc"] + "INFO.OUT"])
                # if it was copied successfully
                if os.path.exists(config["MatLoc"] + "INFO.OUT"):
                    failures = 0
                    # open info file
                    with open(config["MatLoc"] + "INFO.OUT", 'rb') as f:
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
            time.sleep(30)
            print("ELK finished. Coping remote files")
            download_remote(config["MatLoc"])
    else:
        os.system(
            "cd " + config["MatLoc"] + " && . /opt/intel/oneapi/setvars.sh && mpirun /home/wyatt/elk-8.8.26/elk.sh")

    return 0


def run_phonon(mat: pymatgen.core.structure.Structure, config, restart=False):
    if config["remote"]:
        # create directory for this DFT run
        out = 0
        if not restart:
            out = subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "mkdir",
                 config["MatLoc"] + ";", "exit"])
        # if the directory exists, clean it
        if not restart and not out.returncode == 0:
            # if the old file should be saved
            if not config["overwrite"]:
                subprocess.run(
                    ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "zip", "-r", "old/" +
                     config["MatLoc"] + ".zip", config["MatLoc"] + ";", "exit"])
            # remove the old folder
            subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "rm", "-r",
                 config["MatLoc"] + ";", "rm", "-r", config["MatLoc"] + ";", "exit"])
            # create a fresh one
            subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/;", "mkdir",
                 config["MatLoc"] + ";", "exit"])
        # copy elk.in
        subprocess.run(["scp", config["MatLoc"] + "elk.in",
                        "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + config["MatLoc"]
                        + "elk.in"])

        for i in range(1, config["phonon_workers"] + 1):
            config["slurmparams"]["job-name"] = (mat.formula.replace(" ", "") + "_phonon_" + str(i) +
                                                 "/" + str(config["phonon_workers"]))

            prefix = ("#!/bin/bash -l\n#SBATCH --job-name=" + str(config["slurmparams"]["job-name"])
                      + "\n#SBATCH --time=" + config["slurmparams"]["time"] + "\n")

            middle = ('''#SBATCH --ntasks-per-node=48
#SBATCH --nodes=1
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=wbunsti1@jhu.edu
# MARCC defaults to gcc, must load intel module
module load intel/2022.2
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
export OMP_STACKSIZE=50G
ulimit -Ss unlimited\n''')
            suffix = "mpirun /home/wbunsti1/elk/elk-8.8.26/src/elk >& elk_phonon_" + str(i) + "_of_" + str(
                config["phonon_workers"]) + ".log\n"
            file = open(config["MatLoc"] + "elk_phonon_" + str(i) + ".slurm", 'w+')
            file.write(prefix + middle + suffix)
            file.close()
            subprocess.run(["scp", config["MatLoc"] + "elk_phonon_" + str(i) + ".slurm",
                            "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + config[
                                "MatLoc"] + "elk_phonon_" + str(i) + ".slurm"])
            subprocess.run(
                ["ssh", "wbunsti1@login.rockfish.jhu.edu", "cd", "/data/tmcquee2/wbunsti1/" + config["MatLoc"] + ";",
                 "sbatch elk_phonon_" + str(i) + ".slurm; exit"])
    else:
        print("You want to run a phonon calculation remotely")
        return -1
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
        mat = get_mp_mat(config["MatID"])
    elif config["MatType"] == "cif":
        if not "CIF" in config.keys():
            warn("CIF not in config file")
            return -1
        if not os.path.exists(config["CIF"]):
            warn("specified cif, " + config["CIF"] + ", does not exist")
            return -1
        mat = pymatgen.core.structure.Structure.from_file(config["CIF"])
    else:
        warn("MatType not properly specified")
        return -1
    # create the data directory

    if os.path.isdir(config["MatLoc"]):
        if not config["overwrite"]:
            subprocess.run(
                ["zip", "-r", "data/.old/" + config["MatID"] + ".zip", config["MatLoc"]])
        shutil.rmtree(config["MatLoc"])
    os.mkdir(config["MatLoc"])
    if save_cfg:
        with open(config["MatLoc"] + config["MatID"] + '_config.json', 'w') as f:
            json.dump(config, f)

    # ELK functions
    if config["SGU"]: gen_spacegroupin(config["MatLoc"], mat)
    gen_elk_in(config["MatLoc"], mat, config)
    if config["remote"]:
        if "PHONON" in config["tasks"]:
            run_phonon(mat, config)
            return 0
        else:
            gen_slurm(mat, config)
    if not "PHONON" in config["tasks"]:
        match run_elk(config):
            case 0:
                1 + 1
            case -1:
                print("Elk monitoring timeout")
                return -1
            case _:
                print("Unknown return from elk run")
                return -1
    return 0


def download_remote(loc: str, all=False, prefix=None):
    """
    This simply downloads a materials documents from remote.
    :param loc: The material location for which to download
    :return: if download was successful
    """

    #["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/home/wbunsti1/elk/elk-8.8.26/BSPOperDir/" + loc + "*", loc])
    if prefix is not None:
        if isinstance(prefix, str):
            subprocess.run(
                ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + prefix + "*", loc])
        elif isinstance(prefix, list):
            for pref in prefix:
                subprocess.run(
                    ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + pref + "*", loc])
    elif all:
        subprocess.run(
            ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + "*", loc])
    else:
        subprocess.run(
            ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + "INFO.OUT",
             loc])
        subprocess.run(
            ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + "BAND*",
             loc])
        subprocess.run(
            ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + "PDOS*",
             loc])
        subprocess.run(
            ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + "TDOS*",
             loc])
        subprocess.run(
            ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + "elk.log",
             loc])
        subprocess.run(
            ["scp", "-r", "wbunsti1@login.rockfish.jhu.edu:/data/tmcquee2/wbunsti1/" + loc + "GEO*",
             loc])

    return os.path.exists(loc + "BANDLINES.OUT")


def main():
    pause = False
    run = False
    download = True
    plots = True

    elements = ["Cu", "Au", "Sc", "Y"]
    if run:
        cfgs = Utils.MPSearch.get_configs(elements)
        if os.path.exists("running_configs" + str(elements) + ".txt"):
            warn("running_configs" + str(elements) + ".txt exists already, make sure to download and delete file")
            return 0
        running = open("running_configs" + str(elements) + ".txt", "w+")
        if not os.path.exists("completed_runs.txt"):
            tmp = open("completed_runs.txt", "w")
            tmp.close()
        completed = open("completed_runs.txt", "r")
        comp_lines = list(completed.readlines())
        completed.close()
        for config in cfgs:
            if config["MatID"] not in comp_lines:
                if pause:
                    time.sleep(30)
                running.write(config["MatID"] + "\n")
                from_config(config)
        running.close()
    elif download:
        completed = open("completed_runs.txt", "a")
        downloads = open("running_configs" + str(elements) + ".txt", "r")
        for line in downloads.readlines():
            if pause:
                time.sleep(30)
            if download_remote("data/TIProj/" + line[:-1] + "/"):
                completed.write(line)
            else:
                print(line[:-1] + " failed to download")
        completed.close()
        downloads.close()

    if plots:
        i = 0
        completed = open("completed_runs.txt", "r")
        with MPRester(API_KEY) as mpr:
            lines = completed.readlines()
            for line in lines:
                i += 1
                if i > 0:
                    with open("data/TIProj/" + line[:-1] + "/" + line[:-1] + "_config.json") as f:
                        config = json.load(f)
                        mat = mpr.get_structure_by_material_id(line[:-1])
                        plot.plot(config, mat, options=["EDOS", "BS"], spins=True)
                """try:
                    with open("data/TIProj/" + line[:-1] + "/" + line[:-1] + "_config.json") as f:
                        config = json.load(f)
                        mat = mpr.get_structure_by_material_id(line[:-1])
                        plot.plot(config, mat, options=["EDOS", "BS"], spins=True)
                except Exception:
                    print(line[:-1] + " could not be plotted.")"""
        completed.close()


if __name__ == "__main__":
    main()
