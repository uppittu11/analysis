import mdtraj as md
import numpy as np
from .molecules import molecule
from .residue import Residue

def _is_lipid(resname):
    """ Determine whether a residue name is an analyzeable lipid

    Parameters:
    -----------
    resname : string
        name of the residue to check

    Returns:
    --------
    boolean
        True if resnam is in molecule. False otherwise
    """
    return (resname in molecule)

def get_standard_topology(traj):
    """ Load system information into a trajectory.
    Based on the number and name of CG beads in a residue, load in the
    correct residue name

    Parameters:
    -----------
    traj : md.Trajectory
        Trajectory in which system information is loaded.

    Returns:
    --------
    traj : md.Trajectory
        Trajectory with information loaded.
    """
    for i, residue in enumerate(traj.top.residues):
        atoms = set([atom.name for atom in residue.atoms])
        if 'chead' in atoms:
            name = 'chol'
        elif 'head' in atoms:
            name = 'ffa24'
        elif 'oh4' in atoms:
            if 'ter2' in atoms:
                name = 'ucer6'
            else:
                name = 'ecer6'
        elif 'oh3' in atoms:
            if 'ter2' in atoms:
                name = 'ucer3'
            else:
                name = 'ecer3'
        elif 'oh2' in atoms:
            if 'ter2' in atoms:
                name = 'ucer2'
            else:
                name = 'ecer2'
        elif 'water' in atoms:
            name = 'water'
        traj.top.residue(i).name = name

    return traj

def to_residuelist(topology):
    """ Convert a topology into a list of Residue objects.
    Residue objects are intended to reduce memory consumption by
    eliminating unnecessary information.

    Parameters:
    -----------
    topology : md.Topology
        Topology to be convert

    Returns:
    --------
    residuelist : list
        list of Residue objects
    """
    assert type(topology) == md.Topology
    residuelist = []
    for residue in topology:
        if not _is_lipid(residue.name):
            continue
        res_idx = topology.select('residue {}'.format(residue.index))
        tails = []
        for tail_idx in molecule[residue.name][1]:
            tails.append(np.array(res_idx).take(tail_idx))
        new_residue = Residue(name=residue.name, tails=tails)
        residuelist.append(new_residue)
    return residuelist