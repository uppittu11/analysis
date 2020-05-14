import pytest
import numpy as np
import analysis
from .conf_test import ConfTest
import os

class TestEntrypoints(ConfTest):
    def test_entrypoint_atomistic(self):
        traj = self.get_fn("test_aa.dcd")
        top = self.get_fn("test_aa.gro")
        output_dir = self.get_fn("data")

        exit_status = os.system(
                f"analyze -c {top} -f {traj} -o {output_dir} --reload"
                f"> {output_dir+'/log'}")

        assert exit_status == 0

    def test_entrypoint_cg(self):
        traj = self.get_fn("test_cg.dcd")
        top = self.get_fn("test_cg.hoomdxml")
        output_dir = self.get_fn("data")

        exit_status = os.system(
                f"analyze -c {top} -f {traj} -o {output_dir} --reload --cg "
                f"> {output_dir+'/log'}")

        assert exit_status == 0

