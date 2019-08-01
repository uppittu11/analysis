import mdtraj as md
from .molecules import *

def load_system(traj):
    """ load in the correct residue names. Also makes sure that
    the residues are in the list of analyzeable molecules

    Parameters:
    -----------
    traj : mdtraj.Trajectory
        trajectory to load data into

    Returns:
    --------
    traj : mdtraj.Trajectory
    """

    for res in traj.top.residues:
        if not is_lipid(res.name):
            print("residue {} is not in the database".format(res.name))
            assert 1 == 0
    return traj

def is_lipid(resname):
    """ checks to see if the residue is in the list of analyzable
    molecules

    Parameters:
    -----------
    resname : string
        name of residue

    Returns:
    --------
    bool
    """

    return (resname in molecule)

