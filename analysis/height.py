import mdtraj as md
import numpy as np
from .smoothing import savitzky_golay
from scipy.signal import find_peaks

__all__ = ['calc_peaks', 'calc_height']

def calc_peaks(frame, atoms, window=41):
    """ Calculate the locations of peaks in 1-D mass density
    Calculates a mass-weighted density histogram along the z-dimension

    Parameters:
    -----------
    frame : analysis.Frame
        The frame to be analyzed
    atoms : list
        A list of atom indices. These are used to create the
        mass density histogram
    window : int
        Window size for Savizky-Golay smoothing. Must be odd, positive

    Returns:
    --------
    peaks : list
        list of z-coordinates at which there are peaks in the mass
        density histogram
    """

    atoms = np.array(atoms)
    # Heuristic for getting the n_layers from n_leaflets
    n_layers = int(frame.n_leaflets/2+1)
    box_length = frame.unitcell_lengths[2]

    # Collect centered z coordinates and box dimensions
    z = frame.xyz[atoms, 2].reshape(-1) - np.mean(frame.xyz[atoms, 2])
    z_range = [-box_length*0.5-.01,box_length*0.5+.01]

    # Create histogram
    hist, edges = np.histogram(z, weights=frame.masses.take(atoms),
                                range=z_range, bins=400)
    bins = (edges[1:]+edges[:-1]) * 0.5

    # Smoothing via Savitzky-Golay
    hist = savitzky_golay(hist, window, 5)

    # Gets peak indices
    # Prominance: https://en.wikipedia.org/wiki/Topographic_prominence
    peaks, _ = find_peaks(hist, prominence=np.max(hist)*0.25)
    peaks = np.sort(peaks)
    peaks = bins[peaks]

    # Warns if there is an unequal number of peaks and layers
    if len(peaks) != n_layers:
        print("There is an unequal number of peaks " +
                "({}) and layers ({})".format(len(peaks), n_layers))
        # remove the last few peaks if there are too many
        if len(peaks) > n_layers:
            peaks = peaks[:n_layers]
        # adds peaks via linear interpolation if there are too few
        else:
            for _ in range(n_layers-len(peaks)):
                try:
                    np.append(peaks, 2*peaks[-1]-peaks[-2])
                except IndexError:
                    np.append(peaks, peaks[-1])

    return peaks

def calc_height(frame, atoms, window=41):
    """ Calculate the height of layers in frame
    Obtains peak locations the calc_peaks function and takes the
    difference in adjacent peak locations to get the heights.

    Parameters:
    -----------
    frame : analysis.Frame
        The frame to be analyzed
    atoms : list
        A list of atom indices. These are used to create the
        mass density histogram
    window : int
        Window size for Savizky-Golay smoothing. Must be odd, positive

    Returns:
    --------
    height : list
        list of heights for each layer (see above for n_layers)
    """

    peaks = calc_peaks(frame, atoms, window)
    peaks = np.sort(peaks)
    height = peaks[1:] - peaks[:-1]
    return height

