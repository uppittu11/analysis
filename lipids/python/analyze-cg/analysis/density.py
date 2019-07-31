import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from  .molecules import *

def smoothing(hist, bins, window):
    bins_new = [np.mean(bins[i+int(window/2):i+int(3*window/2)]) for i in range(len(bins)-window)]
    hist_new = [np.mean(hist[i+int(window/2):i+int(3*window/2)]) for i in range(len(bins)-window)]
    err_new = [np.std(hist[i+int(window/2):i+int(3*window/2)]) for i in range(len(bins)-window)]
    return np.array(hist_new), np.array(bins_new), np.array(err_new)

def density_plot(traj, n_frames, window, atomselection=None, atoms=None):
    if not atoms:
        atoms = traj.top.select(atomselection) # which atoms to plot

    box_length = np.mean([traj.unitcell_lengths[-(n_frames):,2]])

    hist, edges = np.histogram(traj.xyz[-(n_frames):, atoms, 2].reshape(-1), weights=np.tile(masses.take(atoms), n_frames),
                               range=[-.001,box_length+.001], bins=400)
    bins = (edges[1:]+edges[:-1]) / 2.0
    print(bins[1]-bins[0])
    hist /= 6.02e23 * (10e-9)**3 * 1000 * n_frames
    hist /= np.mean(traj.unitcell_lengths[-(n_frames):,0]**2) * (bins[1]-bins[0]) # divide by volume of slice
    hist, bins, err = smoothing(hist, bins, window)
    return hist, bins, err

def scattering_length_density_plot(traj, n_frames, window, atomselection=None, atoms=None, d=False):
    if not atoms:
        atoms = traj.top.select(atomselection) # which atoms to plot

    box_length = np.mean([traj.unitcell_lengths[-(n_frames):,2]])

    if d:
        s = scattering_length_d
    else:
        s = scattering_length

    hist, edges = np.histogram(traj.xyz[-(n_frames):, atoms, 2].reshape(-1), weights=np.tile(s.take(atoms), n_frames),
                               range=[-.001,box_length+.001], bins=400)
    bins = (edges[1:]+edges[:-1]) / 2.0
    err = bins[1]-bins[0]
    hist /= 6.02e23 * (10e-9)**3 * 1000 * n_frames
    hist /= np.mean(traj.unitcell_lengths[-(n_frames):,0]**2) * (bins[1]-bins[0]) # divide by volume of slice
    hist, bins, err = smoothing(hist, bins, window)
    bins -= box_length/2
    return hist, bins, err

def calc_peaks(frame, atomselection, n_frames, window, n_layers, masses):
    atoms = frame.top.select(atomselection) # which atoms to plot
    box_length = np.mean([frame.unitcell_lengths[0,2]])
    hist, edges = np.histogram(frame.xyz[0, atoms, 2].reshape(-1), weights=np.tile(masses.take(atoms), n_frames),
                               range=[-.01,box_length+.01], bins=400)
    bins = (edges[1:]+edges[:-1]) / 2.0
    hist *= 0.55409730 / 1000 / n_frames
    hist /= np.mean(frame.unitcell_lengths[0,0]**2) * (bins[1]-bins[0]) # divide by volume of slice
    bins *= 6
    hist, bins, err = smoothing(hist, bins, window)
    peaks = [(hist[i], bins[i]) for i in range(len(hist))]
    peaks.sort(reverse=True)
    points_sorted = peaks
    peaks = []

    for i in points_sorted:
        if len(peaks) == 0:
            peaks += [i[1]]
        else:
            add_peak = True
            for j in peaks:
                if np.abs(j-i[1]) < 3:
                    add_peak = False
            if add_peak:
                peaks += [i[1]]
        if len(peaks) == n_layers:
            break

    peaks.sort()
    peaks = np.array(peaks)
    return peaks

def calc_height(frame, atomselection, n_frames, window, n_layers, masses):
    peaks = calc_peaks(frame, atomselection, n_frames, window, n_layers, masses)
    height = np.mean(peaks[1:] - peaks[:-1])
    return height

