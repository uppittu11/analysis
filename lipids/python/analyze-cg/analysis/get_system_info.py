import mdtraj as md
import numpy as np
import sys
import pickle

topfile = sys.argv[1]
outfile = sys.argv[2]

traj = md.load(topfile)
system = []
prev_name = None
counter = 0

for residue in traj.top.residues:
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

    if prev_name == name:
        counter += 1
    else:
        counter += 1
        system += [[prev_name, counter]]
        counter = 0
    prev_name = name

counter += 1
system += [[prev_name, counter]]

system.pop(0)
print(system)
with open(outfile, 'wb') as f:
    pickle.dump(system, f)
