import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET


def calc_com(coords, masses):
    coords, masses = np.array(coords), np.array(masses)
    tot_mass = np.sum(masses)
    com = np.sum(coords * masses.T[np.newaxis,:], axis=1)
    com /= tot_mass
    return com

def calc_moi(coords, masses):
    I = np.zeros((3, 3))
    I[0,0] = np.sum((coords[:, 1]**2 + coords[:, 2]**2 )*masses)
    I[1,1] = np.sum((coords[:, 0]**2 + coords[:, 2]**2 )*masses)
    I[2,2] = np.sum((coords[:, 0]**2 + coords[:, 1]**2 )*masses)
    I[0,1] = I[1,0] = np.sum((coords[:, 0]*coords[:, 1])*masses)
    I[0,2] = I[2,0] = np.sum((coords[:, 0]*coords[:, 2])*masses)
    I[1,2] = I[2,1] = np.sum((coords[:, 1]*coords[:, 2])*masses)
    return I

def calc_director(moi):
    w, v = np.linalg.eig(moi)
    director = v[:,np.argmin(w)]
    return director

