import numpy as np
import mdtraj as md
import pickle as p
import multiprocessing as mp
from analysis import *


def analyze_my_shit(frame, midpoint):

    heads = frame.top.select('name mhead2')
    midpoint = np.mean(frame.xyz[-1, heads, 2])
    heads = [i for i in heads if midpoint-15/60 < frame.xyz[-1,i,2] < midpoint+15/60]
    trip = [[frame.top.atom(i).residue.atom(molecule[frame.top.atom(i).residue.name][2][0][-1]).index,
            i, frame.top.atom(i).residue.atom(molecule[frame.top.atom(i).residue.name][2][1][-1]).index] for i in heads]
    angles = md.compute_angles(frame, trip)
    angles = angles.reshape(-1)
    hairpin = [1  if -90*np.pi/180.0 < i < 90*np.pi/180.0 else 0 for i in angles]
    return np.mean([hairpin])

def main():
    for i in [1, 2, 3]:
        with open('../6-{}-6/traj.p'.format(i), 'rb') as f:
            traj = p.load(f)

        traj = traj[-20:]
        n_frames = traj.n_frames
        print('Loaded trajectory with {} frames'.format(n_frames))

        peaks = np.loadtxt('../6-{}-6/height.txt'.format(i))
        midpoint = np.mean(peaks[:,1])

        # Get parallel processes
        pool = mp.Pool(mp.cpu_count())
        inputs = zip(traj, [midpoint]*len(traj))
        chunksize = int(len(traj)/mp.cpu_count()) + 1
        results = pool.starmap(analyze_my_shit, inputs, chunksize=chunksize)

        # Tidy up the queues
        hairpin = np.array(results)
        print(np.mean(hairpin))

        a = peaks
        b = (a[:,2]-a[:,1])/2*6
        print(np.mean(b))
        print('Printing output files')
        # Save all my shit
        np.savetxt('../6-{}-6/hairpin.txt'.format(i), hairpin, header='all', fmt='%.4f', delimiter='\t')

if __name__ == "__main__": main()

