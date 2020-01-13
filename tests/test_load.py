import pytest
import numpy as np
from os.path import join
from os.path import dirname
import analysis
from .conf_test import ConfTest


class TestLoad(ConfTest):
    def test_load_cg(self):
        topology = join(dirname(__file__), "include/test_cg.hoomdxml")
        trajectory = join(dirname(__file__), "include/test_cg.dcd")
        traj = analysis.load.load_from_trajectory(trajectory, topology)
        assert traj.n_frames == 11
        assert traj.n_atoms == 747
        assert traj.n_residues == 182

    def test_load_cg_notraj(self):
        topology = join(dirname(__file__), "include/test_cg.hoomdxml")
        traj = analysis.load.load_from_trajectory(topology)
        assert traj.n_frames == 1
        assert traj.n_atoms == 747
        assert traj.n_residues == 182

    def test_load_aa(self):
        topology = join(dirname(__file__), "include/test_aa.gro")
        trajectory = join(dirname(__file__), "include/test_aa.dcd")
        traj = analysis.load.load_from_trajectory(trajectory, topology)
        assert traj.n_frames == 10
        assert traj.n_atoms == 2264
        assert traj.n_residues == 153

    def test_load_aa_notraj(self):
        topology = join(dirname(__file__), "include/test_aa.gro")
        traj = analysis.load.load_from_trajectory(topology)
        assert traj.n_frames == 1
        assert traj.n_atoms == 2264
        assert traj.n_residues == 153



