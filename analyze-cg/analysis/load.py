import mdtraj as md
from .molecules import *

def load_system(traj, system):
    start = 0
    for mol in system:
        for index in range(mol[1]):
            traj.topology.residue(index+start).name = mol[0]
        start += mol[1]
    return traj

def is_lipid(resname):
    return (resname in molecule)

def get_standard_topology(traj):
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
