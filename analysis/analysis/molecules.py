import numpy as np

def collect_molecules(cg):
    if cg:
        molecule = {
            'ucer2':(9, [np.arange(8), np.arange(10, 15)], 17),
            'ecer2':(6, [np.arange(5), np.arange(7, 12)], 14),
            'ucer3':(9, [np.arange(8), np.arange(10, 15)], 18),
            'ecer3':(6, [np.arange(5), np.arange(7, 12)], 15),
            'ucer6':(9, [np.arange(3), np.arange(10, 15)], 19),
            'ecer6':(6, [np.arange(5), np.arange(7, 12)], 16),
            'chol' :(0, [np.arange(9)], 9),
            'ffa6' :(0, [np.arange(3)[::-1]], 3),
            'ffa16':(0, [np.arange(6)[::-1]], 6),
            'ffa18':(0, [np.arange(7)[::-1]], 7),
            'ffa20':(0, [np.arange(8)[::-1]], 8),
            'ffa22':(0, [np.arange(9)[::-1]], 9),
            'ffa24':(8, [np.arange(9)[::-1]], 9)
            }
    else:
        molecule = {
            'ucer2':(0, [np.arange(74), np.arange(81, 129)], 129),
            'ecer2':(0, [np.arange(50), np.arange(57, 105)], 105),
            'ucer3':(0, [np.arange(74), np.arange(81, 131)], 132),
            'ecer3':(0,  [np.arange(50), np.arange(57, 107)], 107),
            'ucer6':(0, [np.arange(74), np.arange(81, 131)], 132),
            'ecer6':(0, [np.arange(50), np.arange(51, 107)], 108),
            'chol':(0, [np.arange(74)], 74),
            'ffa6':(5, [np.arange(6)[::-1]], 20),
            'ffa16':(19, [np.arange(16)[::-1]], 50),
            'ffa18':(21, [np.arange(18)[::-1]], 56),
            'ffa20':(23, [np.arange(20)[::-1]], 62),
            'ffa22':(25, [np.arange(22)[::-1]], 68),
            'ffa24':(27, [np.arange(24)[::-1]], 74)
            }
        
    return molecule