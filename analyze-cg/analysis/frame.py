import numpy as np

class Frame(object):
    def __init__(self, xyz=[], unitcell_lengths=[]):
        self.xyz = xyz
        self.unitcell_lengths = unitcell_lengths

    @property
    def xyz(self):
        return self.xyz
    
    @xyz.setter
    def xyz(self, xyz):
        assert type(xyz) == np.ndarray
        self.xyz = xyz
    
    @property
    def unitcell_lengths(self):
        return unitcell_lengths

    @unitcell_lengths.setter
    def unitcell_lengths(self, unitcell_lengths):
        assert type(unitcell_lengths) == np.ndarray
        self.unitcell_lengths = unitcell_lengths