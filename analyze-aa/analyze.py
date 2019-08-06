import analysis.utils
from analysis.molecules import molecule
import mdtraj as md
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
import pickle
import sys

## TO USE THIS SCRIPT ##
# python analyze.py {trajectory} {topology} {output directory} {N leaflets} {is bilayer?}


def analyze_all(frame, masses, n_leaflets, bilayer=False):
    # Prints phase to terminal for each frame. Can be piped to a file and used to
    # track progress
    print('imaframe')

    # Gets the residue Topologies
    residues = [residue for residue in frame.top.residues]

    # Sanitize inputs
    masses = np.array(masses)

    # If the system is a bilayer, calculate properties for top and bottom leaflets
    # separately
    if bilayer:
        # Gets residues from top and bottom leaflets; assumes a symmetrical system
        midpoint = np.mean([frame.xyz[0, :, 2] for r in residues])
        top = [r for r in residues if frame.xyz[0, r.atom(molecule[r.name][0]).index, 2] > midpoint]
        bottom = [r for r in residues if frame.xyz[0, r.atom(molecule[r.name][0]).index, 2] < midpoint]

        ## For the top leaflet:
        # Calculates directors for a given set of residues
        directors = analysis.utils.calc_all_directors(frame, masses, top)

        # Calculate Tilt Angles
        tilt_top = analysis.utils.calc_tilt_angle(directors)

        # Calculate Nematic Order Parameter
        s2_top = analysis.utils.calc_order_parameter(directors)

        # Calculate Area per Lipid: cross section / n_lipids
        apl_top = frame.unitcell_lengths[0, 0]**2 / len(top)


        ## For the top leaflet:
        # Calculates directors for a given set of residues
        directors = analysis.utils.calc_all_directors(frame, masses, bottom)

        # Calculate Tilt Angles
        tilt_bottom = analysis.utils.calc_tilt_angle(directors)

        # Calculate Nematic Order Parameter
        s2_bottom = analysis.utils.calc_order_parameter(directors)

        # Calculate Area per Lipid: cross section / n_lipids
        apl_bottom = frame.unitcell_lengths[0, 0]**2 / len(bottom)

        # Calculate the height -- uses the "head" atoms specified below
        atomselection = '(resname ucer2 ecer2 ucer3 ecer3 ucer6 ecer6 and name N1) or \
                                        (resname chol and name C1) or \
                                        (resname ffa24 and name O27)'
        height = analysis.utils.calc_height(frame, atomselection, int(n_leaflets/2+1), masses)

        return [np.mean(tilt_top), stats.sem(tilt_top), np.mean(tilt_bottom), stats.sem(tilt_bottom),
                np.mean(s2_top), np.mean(s2_bottom), apl_top, apl_bottom, height]

    # If the system is a multilayer, do not attempt to differentiate layers here.
    # Note: if you want to calculate properties for a particular layer, slice it
    # out prior to running this function
    else:
        # Calculates directors for a given set of residues
        directors = analysis.utils.calc_all_directors(frame, masses, residues)

        # Calculate Tilt Angles
        tilt = analysis.utils.calc_tilt_angle(directors)

        # Calculate Nematic Order Parameter
        s2 = analysis.utils.calc_order_parameter(directors)

        # Calculate Area per Lipid: cross section / n_lipids
        apl = frame.unitcell_lengths[0, 0]**2 / len(residues) * n_leaflets

        # Calculate the height -- uses the "head" atoms specified below
        atomselection = '(resname ucer2 ecer2 ucer3 ecer3 ucer6 ecer6 and name N1) or \
                                        (resname chol and name C1) or \
                                        (resname ffa24 and name O27)'
        height = analysis.utils.calc_height(frame, atomselection,int(n_leaflets/2+1), masses)
        return [np.mean(tilt), stats.sem(tilt), np.mean(s2), apl, np.mean(height)]

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

    # If previous traj isn't there load the files inputted via
    # command line
    except:
        try:
            traj = md.load(trajfile, top=topfile)
        except:
            traj = md.load(trajfile)

        # keep only the lipids
        orig_top = traj.top
        select_atoms = traj.top.select("not resname tip3p")
        traj.atom_slice(select_atoms, inplace=True)

        # Load system information
        traj = analysis.load.load_system(traj)

        # generate bonds - atoms the same residue within 1.8 A are bonded
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
    # TODO: add n_electrons for height calculation
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
    print('Starting {} parallel threads'.format(mp.cpu_count()))
    pool = mp.Pool(mp.cpu_count())
    inputs = zip(traj, [masses]*len(traj), [n_leaflets]*len(traj), [bilayer]*len(traj))
    chunksize = int(len(traj)/mp.cpu_count()) + 1
    results = pool.starmap(analyze_all, inputs, chunksize=chunksize)

    print('Cleaning up results')
    results = np.array(results)


    # Dump pickle file of results
    with open('{}/results.p'.format(outputdir), 'wb') as f:
        pickle.dump(results, f)
    print('Finished!')

if __name__ == "__main__": main()
