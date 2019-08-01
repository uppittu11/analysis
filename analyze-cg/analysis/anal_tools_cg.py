import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from .molecules import *
from .angles import *
from .directors import *
from .s2 import *

def calc_all_directors(frame, masses):
    r = [residue for residue in frame.top.residues if residue.name in molecule]
    masses = np.array(masses)
    directors = [[0, 0, 0]]
    for residue in r:
        for atom_idxs in molecule[residue.name][2]:
            atoms = frame.top.select('resid {}'.format(residue.index))
            atoms = np.array(atoms).take(atom_idxs)
            res_coords = frame.xyz[0, atoms]
            res_mass = masses.take([atoms])
            com = calc_com(res_coords, res_mass)
            centered_coords = res_coords - com
            moi = calc_moi(centered_coords, res_mass)
            w, v = np.linalg.eig(moi)
            director = v[:,np.argmin(w)]
            directors.append(director)
    directors.pop(0)
    directors = np.array(directors)
    return directors

def calc_order_parameter(directors):
    Q = calc_q(directors)
    S2 = calc_s2(Q)

    return S2

def calc_tilt_angle(directors):
    vec1 = directors
    vec2 = [[0, 0, 1]] * len(vec1)
    tilt = calc_angle(vec1, vec2)
    return tilt

