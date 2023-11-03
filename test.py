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

ROUNDING = 4
BtoA = 0.52917721090

file = open("data/mp-19009/BAND.OUT")
