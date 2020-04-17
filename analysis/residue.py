import numpy as np

__all__ = ["Residue"]


class Residue(object):
    """
    Class for a single residue

    This residue class is meant to replace the mdtraj Topology object
    with only the vital information. Reduces the amount of memory
    usage.

    Parameters
    ----------
    name : string, optional, default='RES'
        The name of the compound. Should match the names in the mdtraj
        Topology
    tails : list, optional, default=[]
        The indices for each tail of the molecule. This is typically
        taken from the molecules dictionary
    la_regions : list, optional, default=[]
        The indices of the molecule used for performing "linear analysis",
        i.e., the directory, S2, and tilt angle.
    """

    def __init__(self, name="RES", tails=[], la_regions=[]):
        self._name = name
        self._tails = tails
        self._la_regions = la_regions

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def tails(self):
        return self._tails

    @property
    def la_regions(self):
        return self._la_regions
    
    @tails.setter
    def tails(self, tails):
        self._tails = tails

    def add_tail(self, tail):
        tail = np.array(tail)
        self._tails.append(tail)

    @la_regions.setter
    def la_regions(self, la_regions):
        self._la_regions = la_regions
    
    def add_la_regions(self, la_region):
        la_region = np.array(la_region)
        self._la_regions.append(la_region)
