import mdtraj as md
import numpy as np
from math import factorial
from  .molecules import *

def smoothing(y, window_size, order, deriv=0, rate=1):
    """
    Parameters
    ----------
    y : list
        data to be smoothed
    window_size : int
        window size for spline smoothing calc. Must positive and odd
    order : int
        polynomial order. Must be greater than window_size + 2
    deriv : int
        order of derivative. Must be less than or equal to order
    rate : float
        rate coefficient
    Returns
    -------
    y : list
        smoothed list
    """

    order_range = range(order+1)
    half_window = (window_size - 1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window,
                                                           half_window +1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

def calc_height(frame, atomselection, n_layers, masses):
    """
    Calculates the heights of the bilayers present in the system.
    These are generally "head-to-head" heights depending on the
    atomselection.

    Parameters
    ----------
    frame : mdtraj.Trajectory
        simulation frame
    atomselection : list
        list of selected atom indices to use for calculating height.
        These are typically the indices of head atoms
    window : int
        window_size for savitzky golay smoothing filter
    n_layers : int
        number of layers in bilayer. Note that this is the number of
        layer not the number of leaflets
    masses : list
        masses of each of the atoms in the frame.

    Returns
    -------
    heights : list
        list of bilayer heights for each bilayer in the system. These
        are not in any particular order
    """

    atoms = frame.top.select(atomselection) # which atoms to plot
    box_length = np.mean([frame.unitcell_lengths[0,2]])
    hist, edges = np.histogram(frame.xyz[0, atoms, 2].reshape(-1), weights=masses.take(atoms),
                               range=[-.01,box_length+.01], bins=400)
    bins = (edges[1:]+edges[:-1]) / 2.0
    hist = smoothing(hist, 9, 2)
    hist = np.array([(hist[i], bins[i]) for i in range(len(hist))])
    hist = hist[hist[:,0].argsort()][::-1]
    peaks = []

    for i in range(n_layers):
        peaks.append(hist[0])
        hist = hist[np.abs( hist[:,1] - peaks[i][1] ) > 0.5]

    peaks = np.sort(np.array(peaks)[:,1])
    heights = peaks[1:] - peaks[:-1]
    return heights
