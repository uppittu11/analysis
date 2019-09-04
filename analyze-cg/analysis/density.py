import mdtraj as md
import numpy as np
from .molecules import molecule
from .smoothing import savitzky_golay
from scipy.signal import find_peaks

def calc_peaks(frame, atomselection, n_frames, window, n_layers, masses):
    atoms = frame.top.select(atomselection) # which atoms to plot
    box_length = np.mean([frame.unitcell_lengths[0,2]])
    hist, edges = np.histogram(frame.xyz[0, atoms, 2].reshape(-1), weights=np.tile(masses.take(atoms), n_frames),
                               range=[-.01,box_length+.01], bins=400)
    bins = (edges[1:]+edges[:-1]) / 2.0
    hist = savitzky_golay(hist, 40, 5)
    peaks, _ = find_peaks(hist, prominence=np.max(hist)*0.25)
    assert len(peaks) == n_layers
    peaks = bins[peaks]

    return peaks

def calc_height(frame, atomselection, n_frames, window, n_layers, masses):
    peaks = calc_peaks(frame, atomselection, n_frames, window, n_layers, masses)
    peaks = np.sort(peaks)
    height = peaks[1:] - peaks[:-1]
    return height

