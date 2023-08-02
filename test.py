import json
from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen
from pymatgen.symmetry.kpath import *
import numpy as np



MPID = "mp-22862"
path = "data/" + MPID + "/"
ROUNDING = 4
numBandPoints = 400
BtoA = 0.52917721090 #1 bohr radius is 0.52 angstroms

file = open(path+MPID+".json", 'r')
mat = json.load(file)
file.close()
'''
API_KEY = "cZPQqY0nH2aOGBqCGBfbibyF00XJZXWh"
# inputs
with MPRester(API_KEY) as mpr:
    data = mpr.summary.search(material_ids=[MPID])[0]
    # save the relevent data to a dict
    print(data.structure.lattice)
'''


'''
structure = pymatgen.core.structure.Structure.from_dict(mat["structure"])
kpoints = KPathSetyawanCurtarolo(structure).get_kpoints()
print(kpoints)
points = []
for point in range(len(kpoints[0])):
    if not kpoints[1][point] == '':
        coords = []
        for coord in kpoints[0][point]:
            coords.append(round(coord, ROUNDING))
        if len(points) == 0 or not points[-1][0] == kpoints[1][point]:
            latticeCoords = []
            for coord in structure.lattice.get_fractional_coords(coords):
                latticeCoords.append(round(coord, ROUNDING))
            print(latticeCoords)


            points.append([kpoints[1][point],latticeCoords])

'''

structure = pymatgen.core.structure.Structure.from_file("data/MITHRIL_Data/LiZn2Pt.cif")
sga = SpacegroupAnalyzer(structure)
sga.get_primitive_standard_structure()
kpoints = KPathSetyawanCurtarolo(sga.get_primitive_standard_structure()).kpath
print(kpoints)
points = []
for point in kpoints['path'][0]:
    points.append([point, kpoints['kpoints'][point]])
print(points)
