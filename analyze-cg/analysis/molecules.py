import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

global molecule
molecule = {'ucer2':(9, [np.arange(8), np.arange(10, 15)], 17),
            'ecer2':(6, [np.arange(5), np.arange(7, 12)], 14),
            'ucer3':(9, [np.arange(8), np.arange(10, 15)], 18),
            'ecer3':(6, [np.arange(5), np.arange(7, 12)], 15),
            'ucer6':(9, [np.arange(3), np.arange(10, 15)], 19),
            'ecer6':(6, [np.arange(5), np.arange(7, 12)], 16),
            'chol':(0, [np.arange(9)], 9),
            'ffa6':(0, [np.arange(3)[::-1]], 3),
            'ffa16':(0, [np.arange(6)[::-1]], 6),
            'ffa18':(0, [np.arange(7)[::-1]], 7),
            'ffa20':(0, [np.arange(8)[::-1]], 8),
            'ffa22':(0, [np.arange(9)[::-1]], 9),
            'ffa24':(8, [np.arange(9)[::-1]], 9)
           }


