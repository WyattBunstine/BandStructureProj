import sys
import os
from warnings import *
import json
import pymatgen
from mp_api.client import MPRester


def get_MP_Mat(MPID):
    API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
    # inputs
    with MPRester(API_KEY) as mpr:
        data = mpr.summary.search(material_ids=[MPID])
        # save the relevent data to a dict
        return data[0]

def main():
    #check if the config file is given
    if len(sys.argv) < 1:
        warn("no args given")
        return -1
    config = sys.argv[0]
    if not os.path.exists(config):
        warn("config file, "+str(config)+", does not exist")
        return -1
    with open(config) as f:
        config = json.load(f)
    mat = 0
    if config["MatType"] == "MP":
        mat = get_MP_Mat(config["MatLoc"])
    elif config["MatType"] == "cif":
        if not os.path.exists(config["MatLoc"]):
            warn("specified cif does not exist")
            return -1
        mat = pymatgen.core.structure.Structure.from_file(config["MatLoc"])
    else:
        warn("MatType not properly specified")
        return -1
    if not os.path.isdir("data/" + config["MatID"]):
        os.mkdir("data/" + config["MatID"])

    return 0




if __name__ == "__main__":
    main()