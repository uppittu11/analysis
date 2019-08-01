import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from .molecules import *

def load_system(traj, system):
    start = 0
    for mol in system:
        for index in range(mol[1]):
            traj.topology.residue(index+start).name = mol[0]
        start += mol[1]
    return traj

def is_lipid(resname):
    return (resname in molecule)

