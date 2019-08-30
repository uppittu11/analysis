import analysis.utils
from analysis.molecules import molecule
from xml.etree import cElementTree as ET
import mdtraj as md
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
import pickle
import sys

## TO USE THIS SCRIPT ##
# python analyze.py {trajectory_file} {topology_file} {output_directory} {N_leaflets} 


def analyze_all(frame, masses, n_leaflets):
    # Prints phase to terminal for each frame. Can be piped to a file and used to
    # track progress
    print('imaframe')

    # Gets the residue Topologies
    residues = [residue for residue in frame.top.residues]

    # Sanitize inputs
    masses = np.array(masses)

    # Note: if you want to calculate properties for a particular layer, slice it
    # out prior to running this function
    # Calculates directors for a given set of residues
    directors = analysis.utils.calc_all_directors(frame, masses, residues)

    # Calculate Tilt Angles
    tilt = analysis.utils.calc_tilt_angle(directors)

    # Calculate Nematic Order Parameter
    s2 = analysis.utils.calc_order_parameter(directors)

    # Calculate Area per Lipid: cross section / n_lipids
    apl = frame.unitcell_lengths[0, 0]**2 / len(residues) * n_leaflets

    # Calculate the height -- uses the "head" atoms specified below
    atomselection = 'name mhead2 oh1 oh2 oh3 oh4 oh5 amide chead head'
    height = analysis.utils.calc_height(frame, atomselection,int(n_leaflets/2+1), masses)
    return [np.array(tilt), np.mean(s2), apl, np.array(height)]

def main():
    ## PARSING INPUTS
    trajfile = sys.argv[1]
    topfile  = sys.argv[2]
    outputdir = sys.argv[3]
    try:
        n_leaflets = int(sys.argv[4])
    except:
        n_leaflets = 1

    ## LOADING TRAJECTORIES
    # If previous traj exists:
    try:
        traj = md.load('{}/traj.h5'.format(outputdir))
        print("Loading trajectory from {}/traj.h5".format(outputdir))

    # If previous traj isn't there load the files inputted via
    # command line
    except:
        try:
            traj = md.load(trajfile, top=topfile)
            print("Loading trajectory from {} and topology from {}".format(trajfile, topfile))
        except:
            traj = md.load(trajfile)
            print("Loading trajectory from {}".format(trajfile))

        # keep only the lipids
        select_atoms = traj.top.select("not name water")
        traj.atom_slice(select_atoms, inplace=True)

        # Load system information
        traj = analysis.load.get_standard_topology(traj)

        traj.save('{}/traj.h5'.format(outputdir))

    # Get number of frames
    n_frames = traj.n_frames
    print('Loaded trajectory with {} frames'.format(n_frames))

    # Get masses from hoomdxml
    tree = ET.parse(topfile)
    root = tree.getroot()
    masses = np.fromstring(root[0].find('mass').text, sep='\n')
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
    print('Starting {} parallel threads'.format(mp.cpu_count()))
    pool = mp.Pool(mp.cpu_count())
    inputs = zip(traj, [masses]*len(traj), [n_leaflets]*len(traj))
    chunksize = int(len(traj)/mp.cpu_count()) + 1
    results = pool.starmap(analyze_all, inputs, chunksize=chunksize)

    print('Cleaning up results')
    results = np.array(results)


    # Dump pickle file of results
    with open('{}/results.p'.format(outputdir), 'wb') as f:
        pickle.dump(results, f)
    print('Finished!')

if __name__ == "__main__": main()
