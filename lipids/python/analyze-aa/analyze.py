from analysis import *
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import multiprocessing as mp
import pickle
import sys

def analyze_all(frame, masses, n_leaflets, bilayer=False):
    print('imaframe')
    residues = [residue for residue in frame.top.residues if residue.name in molecule]
    masses = np.array(masses)
    directors = [[0, 0, 0]]

    if bilayer:
        midpoint = np.mean([frame.xyz[0, r.atom(molecule[r.name][0]).index, 2] for r in residues])
        top = [r for r in residues if frame.xyz[0, r.atom(molecule[r.name][0]).index, 2] > midpoint]
        bottom = [r for r in residues if frame.xyz[0, r.atom(molecule[r.name][0]).index, 2] < midpoint]

        for residue in top:
            for atom_idxs in molecule[residue.name][2]:
                atoms = frame.top.select('resid {}'.format(residue.index))
                atoms = np.array(atoms).take(atom_idxs)
                res_coords = frame.xyz[0, atoms]
                res_mass = masses.take([atoms])
                com = calc_com(res_coords, res_mass)
                centered_coords = res_coords - com
                moi = calc_moi(centered_coords, res_mass)
                w, v = np.linalg.eig(moi)
                director = v[:,np.argmin(w)]
                directors.append(director)
        directors.pop(0)
        directors = np.array(directors)
        # Calculate Tilt Angles
        tilt_top = calc_tilt_angle(directors)
        # Calculate Nematic Order Parameter
        s2_top = calc_order_parameter(directors)
        apl_top = frame.unitcell_lengths[0, 0]**2 / len(top)

        directors = [[0, 0, 0]]
        for residue in bottom:
            for atom_idxs in molecule[residue.name][2]:
                atoms = frame.top.select('resid {}'.format(residue.index))
                atoms = np.array(atoms).take(atom_idxs)
                res_coords = frame.xyz[0, atoms]
                res_mass = masses.take([atoms])
                com = calc_com(res_coords, res_mass)
                centered_coords = res_coords - com
                moi = calc_moi(centered_coords, res_mass)
                w, v = np.linalg.eig(moi)
                director = v[:,np.argmin(w)]
                directors.append(director)
        directors.pop(0)
        directors = np.array(directors)

        # Calculate Tilt Angles
        tilt_bot = calc_tilt_angle(directors)
        # Calculate Nematic Order Parameter
        s2_bot = calc_order_parameter(directors)
        apl_bot = frame.unitcell_lengths[0, 0]**2 / len(bottom)
        atomselection = '(resname ucer2 ecer2 ucer3 ecer3 ucer6 ecer6 and name N1) or \
                                        (resname chol and name C1) or \
                                        (resname ffa24 and name O27)'
        height = calc_height(frame, atomselection, 1, 5, n_leaflets/2+1, masses)

        return [np.mean(tilt_top), np.std(tilt_top), np.mean(tilt_bot), np.std(tilt_bot),
                np.mean(s2_top), np.mean(s2_bot), apl_top, apl_bot, height]

    else:
        for residue in residues:
            for atom_idxs in molecule[residue.name][2]:
                atoms = frame.top.select('resid {}'.format(residue.index))
                atoms = np.array(atoms).take(atom_idxs)
                res_coords = frame.xyz[0, atoms]
                res_mass = masses.take([atoms])
                com = calc_com(res_coords, res_mass)
                centered_coords = res_coords - com
                moi = calc_moi(centered_coords, res_mass)
                w, v = np.linalg.eig(moi)
                director = v[:,np.argmin(w)]
                directors.append(director)
        directors.pop(0)
        directors = np.array(directors)
        # Calculate Tilt Angles
        tilt = calc_tilt_angle(directors)
        s2 = calc_order_parameter(directors)
        apl = frame.unitcell_lengths[0, 0]**2 / len(residues) * n_leaflets

        atomselection = '(resname ucer2 ecer2 ucer3 ecer3 ucer6 ecer6 and name N1) or \
                                        (resname chol and name C1) or \
                                        (resname ffa24 and name O27)'
        height = calc_height(frame, atomselection, 1, 5, n_leaflets/2+1, masses)
        return [np.mean(tilt), np.std(tilt), np.mean(s2), apl, height]

def main():
    ## PARSING INPUTS
    trajfile = sys.argv[1]
    topfile  = sys.argv[2]
    outputdir = sys.argv[3]
    try:
        n_leaflets = int(sys.argv[4])
    except:
        n_leaflets = 2

    try:
        bilayer = eval(sys.argv[5])
    except:
        bilayer = False

    ## LOADING TRAJECTORIES
    # If previous traj exists:
    try:
        traj = md.load('{}/traj.h5'.format(outputdir))
        print("Loaded trajectory from {}/traj.h5".format(outputdir))

    # If previous traj isn't there load the files inputted
    except:
        try:
            traj = md.load(trajfile, top=topfile)
        except:
            traj = md.load(trajfile)

        # Load system information
        traj = load_system(traj)

        # keep only the lipids
        orig_top = traj.top
        select_atoms = traj.top.select("not resname tip3p")
        traj.atom_slice(select_atoms, inplace=True)

        # generate bonds - nearest neighbor atoms with a residue are bonded
        bonds = []

        for residue in traj.top.residues:
            pairs = [[i.index, j.index] for i in residue.atoms for j in residue.atoms
                    if i.index < j.index and {i.element.name, j.element.name} != {'hydrogen', 'hydrogen'}]
            distances = md.compute_distances(traj[0], pairs, periodic=False)
            cutoff = distances<.18
            take = np.extract(cutoff, np.arange(len(cutoff[0])))
            bonds += np.take(pairs, take, axis=0).tolist()

        for pair in bonds:
            traj.top.add_bond(traj.top.atom(pair[0]), traj.top.atom(pair[1]))

        traj.save('{}/traj.h5'.format(outputdir))

    # Set number of frames
    n_frames = traj.n_frames
    print('Loaded trajectory with {} frames'.format(n_frames))

    # Get masses from element names
    masses = []
    for atom in traj.top.atoms:
        if atom.element.symbol == "C":
            masses.append(12.0107)
        elif atom.element.symbol == "H":
            masses.append(1.00794)
        elif atom.element.symbol == "N":
            masses.append(14.0067)
        elif atom.element.symbol == "O":
            masses.append(15.999)
    print('Loaded masses')

    # Extract atoms within a specified z range
    '''
    z_max = 4.0
    z_min = 3.2
    selected_atoms = [[atom.index for atom in residue.atoms]
            for residue in traj.top.residues
            if residue.name in molecule
            and z_min < np.mean(traj.xyz[:,residue.atom(molecule[residue.name][0]).index,2]) < z_max]
    selected_atoms = [atom for residue in selected_atoms for atom in residue]
    selected_atoms = np.array(selected_atoms)
    traj = traj.atom_slice(selected_atoms)
    masses = masses.take(selected_atoms)
    '''

    # Get parallel processes
    pool = mp.Pool(mp.cpu_count())
    inputs = zip(traj, [masses]*len(traj), [n_leaflets]*len(traj), [bilayer]*len(traj))
    chunksize = int(len(traj)/mp.cpu_count()) + 1
    results = pool.starmap(analyze_all, inputs, chunksize=chunksize)
    results = np.array(results)

    # Dump pickle file of results
    with open('{}/results.p'.format(outputdir), 'wb') as f:
        pickle.dump(results, f)

if __name__ == "__main__": main()
