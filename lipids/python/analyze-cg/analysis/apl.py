import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from  .molecules import *

def calc_apl(frame, leaflets=2):
    atom_idxs = [molecule[residue.name][0]+residue.atom(0).index for residue in frame.topology.residues if residue.name in molecule]
    midpoint = np.mean([frame.xyz[0,atom_idxs,2]])
    count_top = len([None for atom in atom_idxs if frame.xyz[0, atom, 2]>midpoint])
    apl_top = frame.unitcell_lengths[0,0]*frame.unitcell_lengths[0,1] / count_top
    count_bot = len([None for atom in atom_idxs if frame.xyz[0, atom, 2]<midpoint])
    apl_bot = frame.unitcell_lengths[0,0]*frame.unitcell_lengths[0,1] / count_bot
    apl = [apl_top, apl_bot]
    return apl


