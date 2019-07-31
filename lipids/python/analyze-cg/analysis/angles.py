import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET


def calc_angle(vec1, vec2):
    dot = np.sum(np.array(vec1)*np.array(vec2), axis=1)
    theta = np.arccos(dot)*180/np.pi
    for i, angle in enumerate(theta):
        if angle > 90:
            theta[i] = 180 - theta[i]
    return theta

def calc_direction_vector(coord1, coord2):
    vector = np.array(coord2)-np.array(coord1)
    magnitude = np.sqrt(np.sum(vector**2, axis=1))
    vector = vector/magnitude[:,np.newaxis]
    return vector

