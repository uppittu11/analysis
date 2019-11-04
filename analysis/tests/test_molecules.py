import pytest
from .conf_test import ConfTest


class TestMolecules(ConfTest):
    def test_load_cg_defaults(self):
        from analysis import collect_molecules

        mols = collect_molecules("cg")
        assert len(mols) > 0

    def test_load_atomistic_defaults(self):
        from analysis import collect_molecules

        mols = collect_molecules("atomistic")
        assert len(mols) > 0

    def test_molecule_validator(self):
        from analysis.molecules.molecules import Molecule
        molecule = Molecule()
        molecule.name = "Ceramide"
        molecule.add_tail([1, 2, 3])
        molecule.head = 0
        with pytest.raises(TypeError):
            molecule.name = [1, 2, 3]

