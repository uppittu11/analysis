import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

def calc_q(directors):
    directors = np.array(directors)
    Q = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Q[i,j] = np.sum(directors[:,i]*directors[:,j]*3)
            if i == j:
                Q[i,j] -= len(directors)
    Q /= 2 * len(directors)
    return Q

def calc_s2(q):
    w, v = np.linalg.eig(q)
    return max(w)

