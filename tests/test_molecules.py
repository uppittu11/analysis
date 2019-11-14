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

    def test_name(self):
        from analysis.molecules.molecules import Molecule
        molecule = Molecule()
        molecule.name = "Ceramide"
        assert molecule.name == "Ceramide"
        with pytest.raises(TypeError):
            molecule.name = [1, 2, 3]

    def test_head(self):
        from analysis.molecules.molecules import Molecule
        molecule = Molecule()
        molecule.head = 5
        assert molecule.head == 5
        with pytest.raises(TypeError):
            molecule.head = 2.5

    def test_n_atoms(self):
        from analysis.molecules.molecules import Molecule
        molecule = Molecule()
        molecule.n_atoms = 3
        assert molecule.n_atoms == 3
        with pytest.raises(TypeError):
            molecule.n_atoms = "three"

    def test_tails(self):
        from analysis.molecules.molecules import Molecule
        molecule = Molecule()
        molecule.add_tail([1, 2, 3])
        molecule.add_tail((4, 5, 6))
        molecule.add_tail([1, 2, 3])
        molecule.remove_tail([1, 2, 3])
        molecule.remove_tail(1)
        assert molecule.tails == [(4, 5, 6)]
        with pytest.raises(TypeError):
            molecule.add_tail(5)

