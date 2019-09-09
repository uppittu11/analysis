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

    def __init__(self, xyz=[], unitcell_lengths=[], masses=[], 
                    residuelist=[], atomnames=[], n_leaflets=2,
                    cg=False):
        self._xyz = xyz
        self._unitcell_lengths = unitcell_lengths
        self._masses = masses
        self._residuelist = residuelist
        self._n_leaflets = n_leaflets
        self._atomnames = atomnames
        self._cg = cg

    @property
    def xyz(self):
        return self._xyz
    
    @xyz.setter
    def xyz(self, xyz):
        self._xyz = xyz
    
    @property
    def unitcell_lengths(self):
        return self._unitcell_lengths

    @unitcell_lengths.setter
    def unitcell_lengths(self, unitcell_lengths):
        self._unitcell_lengths = unitcell_lengths

    @property
    def masses(self):
        return self._masses
    
    @masses.setter
    def masses(self, masses):
        self._masses = masses

    @property
    def residuelist(self):
        return self._residuelist
    
    @residuelist.setter
    def residuelist(self, residuelist):
        self._residuelist = residuelist
    
    @property
    def atomnames(self):
        return self._atomnames
    
    @atomnames.setter
    def atomnames(self, atomnames):
        self._atomnames = atomnames

    @property
    def n_leaflets(self):
        return self._n_leaflets
    
    @n_leaflets.setter
    def n_leaflets(self, n_leaflets):
        n_leaflets == int(n_leaflets)
        self._n_leaflets = n_leaflets
    
    @property
    def cg(self):
        return self._cg
    
    @cg.setter
    def cg(self, cg):
        self._cg = cg

    def __repr__(self):
        return "<Frame with {} residues and {} atoms>".format(
                    len(self._residuelist), len(self.masses))
    
    def validate_frame(self):
        """ Ensure that this frame can be analyzed.
        Verifies that the number of coordinates, masses, and residues 
        match and that the dimensions/type of unitcell lengths and 
        n_leaflets are the correct.

        Notes
        -----
        Raises an assertion error in the case that there is a mismatch
        in dimensions/types
        """
        self._masses = np.array(self._masses)
        self._unitcell_lengths = np.array(self._unitcell_lengths)
        self._xyz = np.array(self._xyz)
        self._atomnames = np.array(self.atomnames)

        assert self._masses.shape[0] == self._xyz.shape[0]
        assert self._masses.shape[0] == self._atomnames.shape[0]
        assert self._unitcell_lengths.shape[0] == 3
        assert type(self._n_leaflets) == int
    
    def select(self, names=None, mass_range=None):
        if names:
            return np.array([index for index, atom in 
                                enumerate(self._atomnames) if 
                                atom in set(names)])
        elif mass_range:
            mass_range = np.array(mass_range)
            assert mass_range.shape[-1] == 2
            return np.array([index for index, mass in 
                                enumerate(self.masses) if
                                mass_range[0] < 
                                mass
                                < mass_range[1]])