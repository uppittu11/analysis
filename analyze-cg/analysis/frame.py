import numpy as np

class Frame(object):
    """
    Class for a single frame

    Each frame contains the coordinates of each atom, unitcell lengths
    and the mass of each atom

    Parameters
    ----------
    xyz : np.ndarray, optional, default=[]
        The coordinates of each atom. N x 3 array
    unitcell_lengths : np.ndarray, optional, default=[]
        The box lengths for the frame. 3 element array
    masses : np.ndarray, optional, default=[]
        The masses of each atom. N element array
    residuelist : list, optional, default=[]
        List of Residue objects
    n_leaflets : int, optional, default=2
        Number of leaflets in the frame

    """

    def __init__(self, xyz=[], unitcell_lengths=[], masses=[], residuelist=[], n_leaflets=2):
        self.xyz = xyz
        self.unitcell_lengths = unitcell_lengths
        self.masses = masses
        self.residuelist = residuelist

    @property
    def xyz(self):
        return self.xyz
    
    @xyz.setter
    def xyz(self, xyz):
        assert type(xyz) == np.ndarray
        self.xyz = xyz
    
    @property
    def unitcell_lengths(self):
        return self.unitcell_lengths

    @unitcell_lengths.setter
    def unitcell_lengths(self, unitcell_lengths):
        assert type(unitcell_lengths) == np.ndarray
        self.unitcell_lengths = unitcell_lengths

    @property
    def masses(self):
        return self.masses
    
    @masses.setter
    def masses(self, masses):
        assert type(masses) == np.ndarray
        self.masses = masses

    @property
    def residuelist(self):
        return self.residuelist
    
    @residuelist.setter
    def residuelist(self, residuelist):
        assert type(residuelist) == np.ndarray
        self.residuelist = residuelist

    @property
    def n_leaflets(self):
        return self.n_leaflets
    
    @n_leaflets.setter
    def n_leaflets(self, n_leaflets):
        n_leaflets == int(n_leaflets)
        self.n_leaflets = n_leaflets
    
    def validate_frame(self):
        self.masses = np.array(self.masses)
        self.unitcell_lengths = np.array(self.unitcell_lengths)
        self.xyz = np.array(self.xyz)

        assert self.masses.shape[0] == self.unitcell_lengths.shape[0]
        assert self.unitcell_lengths.shape[0] == self.xyz.shape[0]
        assert self.xyz.shape[0] == len(self.residuelist)