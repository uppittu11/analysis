import mdtraj as md
import numpy as np
from .molecules import molecule
from .smoothing import savitzky_golay
from scipy.signal import find_peaks

def calc_peaks(frame, atomselection, window=41):
    n_layers = int(frame.n_leaflets/2+1)
    atoms = frame.select(atomselection) # which atoms to plot
    box_length = frame.unitcell_lengths[2]

    z = frame.xyz[atoms, 2].reshape(-1) - np.mean(frame.xyz[atoms, 2])
    z_range = [-box_length*0.5-.01,box_length*0.5+.01]

    hist, edges = np.histogram(z, weights=frame.masses.take(atoms),
                                range=z_range, bins=400)
    bins = (edges[1:]+edges[:-1]) * 0.5
    hist = savitzky_golay(hist, window, 5)
    peaks, _ = find_peaks(hist, prominence=np.max(hist)*0.25)
    peaks = np.sort(peaks)
    peaks = bins[peaks]

    if len(peaks) != n_layers:
        print("There is an unequal number of peaks ({}) and layers ({})".format(len(peaks), n_layers))
        if len(peaks) > n_layers:
            peaks = peaks[:n_layers]
        else:
            for _ in range(n_layers-len(peaks)):
                peaks.append(2*peaks[-1]-peaks[-2])
    

    return peaks

def calc_height(frame, atomselection, window=41):
    peaks = calc_peaks(frame, atomselection, window)
    peaks = np.sort(peaks)
    height = peaks[1:] - peaks[:-1]
    return height

