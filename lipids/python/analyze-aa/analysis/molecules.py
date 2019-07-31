import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

# format for molecule:
# 'name':(head_atom, tail_atom, list_of_arrays_of_tail_atoms, number_of_atoms)

global molecule
molecule = {'ucer2':(0, 125, [np.arange(74), np.arange(81, 129)], 129),
            'ecer2':(0, 101, [np.arange(50), np.arange(57, 105)], 105),
            'ucer3':(0, 129, [np.arange(74), np.arange(81, 131)], 132),
            'ecer3':(0, 105, [np.arange(50), np.arange(57, 107)], 107),
            'ucer6':(0, 129, [np.arange(74), np.arange(81, 131)], 132),
            'ecer6':(0, 105, [np.arange(50), np.arange(51, 107)], 108),
            'chol':(0, 73, [np.arange(74)], 74),
            'ffa6':(5, 0, [np.arange(6)[::-1]], 20),
            'ffa16':(19, 0, [np.arange(16)[::-1]], 50),
            'ffa18':(21, 0, [np.arange(18)[::-1]], 56),
            'ffa20':(23, 0, [np.arange(20)[::-1]], 62),
            'ffa22':(25, 0, [np.arange(22)[::-1]], 68),
            'ffa24':(27, 0, [np.arange(24)[::-1]], 74)
           }


