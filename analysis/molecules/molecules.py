import numpy as np
import os
from pkg_resources import resource_filename
import glob
import json

__all__ = ["collect_molecules", "Molecule"]


def collect_molecules(library_dir, defaults=True):
    if defaults:
        library_dir = resource_filename("analysis", "molecules/{}/".format(library_dir))

    library = _load_jsons(library_dir)

    return library


def _load_json(fname):
    with open(fname, "r") as f:
        molecule = Molecule()
        molecule_dict = json.load(f)
        molecule.from_dict(molecule_dict)
    return {molecule.name : molecule}

def _load_jsons(library_dir):
    json_files = glob.glob(os.path.join(library_dir, "*.json"))
    library = dict()
    for fname in json_files:
        library.update(_load_json(fname))
    return library

class Molecule(object):
    def __init__(self, name="Molecule", head=0, tails=[], n_atoms=0):
        self._name = name
        self._head = head
        self._tails = tails
        self._n_atoms = n_atoms
        self._validate_molecule()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name
        self._validate_molecule()

    @property
    def head(self):
        return self._head

    @head.setter
    def head(self, head):
        self._head = head
        self._validate_molecule()

    @property
    def tails(self):
        return self._tails

    def add_tail(self, tail):
        self._tails.append(tail)
        self._validate_molecule()

    def remove_tail(self, tail):
        if isinstance(tail, int):
            del self._tails[tail]
        elif isinstance(tail, list):
            self._tails.remove(tail)

    @property
    def n_atoms(self):
        return self._n_atoms

    @name.setter
    def n_atoms(self, n_atoms):
        self._n_atoms = n_atoms
        self._validate_molecule()

    def from_dict(self, molecule_dict):
        try:
            assert isinstance(molecule_dict, dict)
        except:
            raise TypeError("Argument 'molecule_dict' must be dict type")

        for key in molecule_dict:
            setattr(self, '_'+str(key), molecule_dict[key])

        self._validate_molecule()

    def _validate_molecule(self):
        try:
            assert isinstance(self._name, str)
        except AssertionError:
            raise TypeError("Attribute 'name' must be str type")

        try:
            assert isinstance(self._head, int)
        except AssertionError:
            raise TypeError("Attribute 'head' must be int type")

        try:
            assert isinstance(self._tails, list)
        except AssertionError:
            raise TypeError("Attribute 'tails' must be list type")

        try:
            assert isinstance(self._n_atoms, int)
        except AssertionError:
            raise TypeError("Attribute 'n_atoms' must be int type")
