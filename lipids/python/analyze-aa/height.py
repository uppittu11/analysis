from analysis import *
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from xml.etree import cElementTree as ET
import multiprocessing as mp
import pickle
import sys
import matplotlib.pyplot as plt

def analyze_my_shit(frame, masses):
    # Calculates parameters for top and bottom if True
    n_leaflets = 4
    bilayer = False

    residues = [residue for residue in frame.top.residues if residue.name in molecule]
    # Calculate directors for tails
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

        atomselection = 'name mhead2 head chead'
        height = calc_height(frame, atomselection, 1, 5, n_leaflets/2+1, masses)
        return [np.mean(tilt_top), np.std(tilt_top), np.mean(tilt_bot), np.std(tilt_bot),
                np.mean(s2_top), np.std(s2_top), np.mean(s2_bot), np.std(s2_bot), apl_top, apl_bot, height]

    else:
        """
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
        tilt_all = calc_tilt_angle(directors)
        # Calculate Nematic Order Parameter
        s2_all = calc_order_parameter(directors)
        apl_all = frame.unitcell_lengths[0, 0]**2 / len(residues) * n_leaflets
        """
        atomselection = 'name mhead2 head chead'
        peaks = calc_peaks(frame, atomselection, 1, 5, n_leaflets/2+1, masses)
        """
        return [np.mean(tilt_all), np.std(tilt_all),
                np.mean(s2_all), np.std(s2_all), apl_all, peaks]
        """
        return peaks
def main():
    trajfile = sys.argv[1]
    try:
        topfile  = sys.argv[2]
        outputdir = sys.argv[3]
    except:
        topfile = sys.argv[1]
        outputdir = sys.argv[2]

    # Loading trajectories
    # If previous traj exists:
    try:
        with open("{}/traj.p".format(outputdir),"rb") as file:
            traj = pickle.load(file)
    # If previous traj isn't there
    except:
        if trajfile != topfile:
            traj = md.load(trajfile, top=topfile)[1:]
        else:
            traj = md.load(trajfile)

        #anchor_molecules = [mol for mol in traj.top.find_molecules() if len(mol) > 2]
        #traj.image_molecules(inplace=True, anchor_molecules=anchor_molecules)
        """
        # Load system information
        with open(systeminfo, 'rb') as f:
            system = pickle.load(f)
        """
        traj = get_standard_topology(traj)

        with open("{}/traj.p".format(outputdir),"wb") as file:
            pickle.dump(traj, file)

    n_frames = traj.n_frames
    print('Loaded trajectory with {} frames'.format(n_frames))

    # Get masses from hoomdxml
    tree = ET.parse(topfile)
    root = tree.getroot()
    masses = np.fromstring(root[0][3].text, sep='\n')
    #masses = [1.0 for i in range(traj.n_atoms)]
    #select_atoms = traj.top.select("resname ucer2 ecer2 ucer3 ecer3 ucer6 ecer6 chol ffa24")
    #traj.atom_slice(select_atoms, inplace=True)
    #masses = np.take(masses, select_atoms)
    '''
    # Extract atoms within a range
    selected_atoms = [[atom.index for atom in residue.atoms]
            for residue in traj.top.residues
            if residue.name in molecule
            and 3.2 < np.mean(traj.xyz[:,residue.atom(molecule[residue.name][0]).index,2]) < 4]

    selected_atoms = np.array(selected_atoms).reshape(-1)
    traj = traj.atom_slice(selected_atoms)
    masses = masses.take(selected_atoms)
    '''

    # Get parallel processes
    pool = mp.Pool(mp.cpu_count())
    inputs = zip(traj, [masses]*len(traj))
    chunksize = int(len(traj)/mp.cpu_count()) + 1
    results = pool.starmap(analyze_my_shit, inputs, chunksize=chunksize)

    # Tidy up the queues
    results = np.array(results)
    #apl = results[:, 4] * 36
    #tilt = results[:, :2]
    #s2 = results[:, 2:4]
    #height = results[:, 5]

    height = np.array(results)

    print('Printing output files')
    # Save all my shit
    #np.savetxt('{}/apl.txt'.format(outputdir), apl, header='all', fmt='%.4f', delimiter='\t')
    #np.savetxt('{}/tilt.txt'.format(outputdir), tilt, header='all', fmt='%.4f', delimiter='\t')
    #np.savetxt('{}/s2.txt'.format(outputdir), s2, header='all', fmt='%.4f', delimiter='\t')
    np.savetxt('{}/height.txt'.format(outputdir), height, header='all', fmt='%.4f', delimiter='\t')

    #tilt = tilt[:,0]
    #s2 = s2[:,0]

    #with open('{}/summary.txt'.format(outputdir), 'w') as f:
        #f.write('OVERALL\n')
        #f.write('Area per Lipid:\t{0:.3f} +/- {1:.3f}\n'.format(np.mean(apl), np.std(apl)))
        #f.write('Tilt Angle:\t{0:.3f} +/- {1:.3f}\n'.format(np.mean(tilt), np.std(tilt)))
        #f.write('S2:\t\t{0:.3f} +/- {1:.3f}\n'.format(np.mean(s2), np.std(s2)))
        #f.write('Height:\t{0:.3f} +/- {1:.3f}\n'.format(np.mean(height), np.std(height)))
        #f.write("\n")
if __name__ == "__main__": main()
