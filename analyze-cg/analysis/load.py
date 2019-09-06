import mdtraj as md
import numpy as np
from .molecules import collect_molecules
from .residue import Residue
from xml.etree import cElementTree as ET

def _is_lipid(resname, cg):
    molecule = collect_molecules(cg)
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

def get_cg_residuename(residue):
    name = None
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
    if name == None:
        raise KeyError("Residue {} ".format(residue.name) +
                        "is not in the database")
    return name

def get_standard_topology(traj, cg):
    """ Load system information into a trajectory.
    CG systems: Based on the number and name of CG beads in a residue, 
    load in the correct residue name

    AA systems: Validate that all residue names are in the molecule
    list.

    Parameters:
    -----------
    traj : md.Trajectory
        Trajectory in which system information is loaded.
    cg :  boolean
        True if the system is coarse-grained. False for atomistic

    Returns:
    --------
    traj : md.Trajectory
        Trajectory with information loaded.
    """

    if cg:
        for i, residue in enumerate(traj.top.residues):
            name = get_cg_residuename(residue)
            traj.top.residue(i).name = name
    else:
        for res in traj.top.residues:
            if not _is_lipid(res.name, cg):
                raise KeyError("Residue {} ".format(res.name) +
                                "is not in the database")

    return traj

def load_masses(cg, topology=None, topfile=None):
    if cg:
        if topfile == None: 
            raise ValueError("topfile is a required arguement " +
                                "if cg=True")
        tree = ET.parse(topfile)
        root = tree.getroot()
        masses = np.fromstring(root[0].find('mass').text, sep='\n')
    else:
        if topology == None: 
            raise ValueError("topology is a required arguement " +
                                "if cg=True")
        masses = []
        for atom in topology.atoms:
            if atom.element.symbol == "C":
                masses.append(12.0107)
            elif atom.element.symbol == "H":
                masses.append(1.00794)
            elif atom.element.symbol == "N":
                masses.append(14.0067)
            elif atom.element.symbol == "O":
                masses.append(15.999)
    return masses

def to_residuelist(topology, cg):
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
    molecule = collect_molecules(cg)
    assert type(topology) == md.Topology
    residuelist = []
    for residue in topology.residues:
        if not _is_lipid(residue.name, cg=cg):
            continue
        res_idx = topology.select('residue {}'.format(residue.index))
        tails = []
        for tail_idx in molecule[residue.name][1]:
            tails.append(np.array(res_idx).take(tail_idx))
        new_residue = Residue(name=residue.name, tails=tails)
        residuelist.append(new_residue)
    return residuelist