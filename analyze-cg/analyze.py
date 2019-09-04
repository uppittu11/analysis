import multiprocessing as mp
import pickle
import sys
from optparse import OptionParser
from xml.etree import cElementTree as ET
import numpy as np
import scipy.stats as stats
import mdtraj as md
import analysis
from analysis.frame import Frame
from analysis.molecules import molecule
import copy as cp

## TO USE THIS SCRIPT ##
# python analyze.py {trajectory_file} {topology_file} {output_directory} {N_leaflets} 


def analyze_all(frame):
    # Prints frame number to terminal for each frame. Can be piped to a file and used to
    # track progress
    print('imaframe')

    # Note: if you want to calculate properties for a particular layer, slice it
    # out prior to running this function

    # Unpack inputs
    frame.validate_frame()
    
    # Calculates directors for a given set of residues
    directors = analysis.utils.calc_all_directors(frame.xyz, frame.masses, frame.residuelist)

    # Calculate Tilt Angles
    tilt = analysis.utils.calc_tilt_angle(directors)

    # Calculate Nematic Order Parameter
    s2 = analysis.utils.calc_order_parameter(directors)

    # Calculate Area per Lipid: cross section / n_lipids
    apl = frame.unitcell_lengths[0] * frame.unitcell_lengths[1] / len(frame.residuelist) * frame.n_leaflets

    # Calculate the height -- uses the "head" atoms specified below
    atomselection = 'name mhead2 oh1 oh2 oh3 oh4 oh5 amide chead head'
    height = analysis.density.calc_height(frame.xyz, atomselection, int(frame.n_leaflets/2+1), frame.masses)
    results = {'tilt' :  np.array(tilt),
                's2' : s2,
                'apl' : apl,
                'height' : np.array(height)}
    return results

def main():
    ## PARSING INPUTS
    parser = OptionParser()
    parser.add_option("-f", "--file", action="store", type="string", dest="trajfile")
    parser.add_option("-c", "--conf", action="store", type="string", dest="topfile")
    parser.add_option("-o", "--output", action="store", type="string", dest="outputdir", default="./")
    parser.add_option("-n", "--nleaflets", action="store", type="int", dest="n_leaflets",  default=1)
    (options, _) = parser.parse_args()

    trajfile = options.trajfile
    topfile  = options.topfile
    outputdir = options.outputdir
    n_leaflets = options.n_leaflets

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
    # Convert to Frame/residuelist format
    residuelist = cp.deepcopy(analysis.load.to_residuelist(traj.top))
    frames = []
    for i in range(n_frames):
        frame = Frame(xyz=np.squeeze(traj.xyz[i,:,:]),
                unitcell_lengths=np.squeeze(traj.unitcell_lengths[i,:]),
                masses=masses, residuelist=residuelist, 
                n_leaflets=n_leaflets)
        frames.append(cp.deepcopy(frame))
    
    # Purge the old trajectory from memory
    del traj

    # Get parallel processes
    print('Starting {} parallel threads'.format(mp.cpu_count()))
    pool = mp.Pool(mp.cpu_count())
    chunksize = int(len(frames)/mp.cpu_count()) + 1
    results = pool.starmap(analyze_all, frames, chunksize=chunksize)

    print('Cleaning up results')
    results = np.array(results)


    # Dump pickle file of results
    with open('{}/results.p'.format(outputdir), 'wb') as f:
        pickle.dump(results, f)
    print('Finished!')

if __name__ == "__main__": main()
