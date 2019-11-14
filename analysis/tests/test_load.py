import pytest
import numpy as np
import analysis
from .conf_test import ConfTest


class TestLoad(ConfTest):
    def test_load_cg(self):
        topology = "./include/test_cg.hoomdxml"
        trajectory = "./include/test_cg.dcd"
        traj = analysis.load.load_from_trajectory(trajectory, topology)
        assert traj.n_frames == 11
        assert traj.n_atoms == 747
        assert traj.n_residues == 182

    def test_load_cg_notraj(self):
        topology = "./include/test_cg.hoomdxml"
        traj = analysis.load.load_from_trajectory(topology)
        assert traj.n_frames == 1
        assert traj.n_atoms == 747
        assert traj.n_residues == 182



