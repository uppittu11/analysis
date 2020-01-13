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
    def __init__(self, name="Molecule", head=0, tails=None, n_atoms=0):
        self._name = name
        self._head = head
        if tails == None:
            self._tails = []
        else:
            for index, tail in enumerate(tails):
                tails[index] = self._validate_tail(tail)
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
        tail = self._validate_tail(tail)
        self._tails.append(tail)
        self._validate_molecule()

    def remove_tail(self, tail):
        # if the tail is an int the tail at that index is removed
        if isinstance(tail, int):
            del self._tails[tail]
        # it the tail is iterable, the first matching tail is removed
        else:
            tail = self._validate_tail(tail)
            self._tails.remove(tail)

    def _validate_tail(self, tail):
        # check if tail is an iterable
        try:
            iter(tail)
            tail = tuple(tail)
        except TypeError:
            raise TypeError("tail must be iterable")

        # check if indices are ints
        try:
            for index in tail: assert isinstance(index, int)
        except TypeError:
            raise TypeError("Indices must be ints")

        return tail

    @property
    def n_atoms(self):
        return self._n_atoms

    @n_atoms.setter
    def n_atoms(self, n_atoms):
        self._n_atoms = n_atoms
        self._validate_molecule()

    def from_dict(self, molecule_dict):
        try:
            assert isinstance(molecule_dict, dict)
        except:
            raise TypeError("Argument 'molecule_dict' must be dict type")

        for key in molecule_dict:
            if key == "name":
                self.name = molecule_dict[key]
            elif key == "n_atoms":
                self.n_atoms = molecule_dict[key]
            elif key == "tails":
                for tail in molecule_dict[key]:
                    self.add_tail(tail)
            elif key == "head":
                self.head = molecule_dict[key]
            else:
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
            for tail in self._tails:
                assert isinstance(tail, tuple)
        except AssertionError:
            raise TypeError("Attribute 'tails' must be list type")

        try:
            assert isinstance(self._n_atoms, int)
        except AssertionError:
            raise TypeError("Attribute 'n_atoms' must be int type")

    def __repr__(self):
        return "<Molecule {}, {} atoms, id: {}>".format(self._name,
                self.n_atoms, id(self))
