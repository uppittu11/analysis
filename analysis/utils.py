import mdtraj as md
import numpy as np

# imports all python files
from analysis.angles import calc_angle
from analysis.directors import calc_com, calc_moi, calc_director
from analysis.s2 import calc_q, calc_s2

__all__ = ["calc_all_directors", "calc_tilt_angle", "calc_order_parameter"]


def calc_all_directors(xyz, masses, residues, return_coms=False):
    """ Calculates directors for all residues in a frame. This is
    a wrapper for the calc_director function which only works for
    a single residue

    Parameters:
    -----------
    frame : mdtraj.Trajectory
        frame to analyze
    masses : list
        list of masses corresponding to each bead in the frame
    com : boolean
        Returns the COM for each tail if true

    Returns:
    --------
    directors : list
        list of directors
    """
    masses = np.array(masses)

    n_tails = sum([len(residue.tails) for residue in residues])
    coms = np.zeros((n_tails, 3), dtype=np.float)
    directors = np.zeros((n_tails, 3), dtype=np.float)

    tail_num = 0
    for residue in residues:
        for tail_atoms in residue.tails:

            res_coords = xyz[tail_atoms]
            res_mass = masses.take([tail_atoms])

            com = calc_com(res_coords, res_mass)
            centered_coords = res_coords - com
            moi = calc_moi(centered_coords, res_mass)

            director = calc_director(moi)

            coms[tail_num] = com
            directors[tail_num] = director

            tail_num = tail_num + 1

    assert tail_num == n_tails

    results_dict = {"directors": directors}
    if return_coms:
        results_dict.update({"coms": coms})

    return results_dict

    def tail_worker(atoms):
        """ worker function for calculating a director. This allows for
        list comprehension

        Parameters:
        -----------
        atom_idxs : list
            list of indices corresponding to a the tail

        Returns:
        --------
        director : list
            returns a single director for a tail
        """

        res_coords = xyz[atoms]
        res_mass = masses.take([atoms])
        com = calc_com(res_coords, res_mass)
        centered_coords = res_coords - com
        moi = calc_moi(centered_coords, res_mass)
        director = calc_director(moi)
        if return_com:
            return [director, com]
        else:
            return [director]

    tail_idxs = [tail for residue in residues for tail in residue.tails]
    results = [tail_worker(atom_indices) for atom_indices in tail_idxs]
    if com:
        directors = [result[0] for result in results]
        directors = np.array(directors)
        coms = [result[1] for result in results]
        coms = np.array(coms)
        return directors, coms
    return directors


def calc_order_parameter(directors):
    Q = calc_q(directors)
    S2 = calc_s2(Q)
    return S2


def calc_tilt_angle(directors):
    vec1 = directors
    vec2 = [[0, 0, 1]] * len(vec1)
    tilt = calc_angle(vec1, vec2)
    return tilt
