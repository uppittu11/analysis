import numpy as np
import os
from pkg_resources import resource_filename
import glob
import json

__all__ = ["collect_molecules"]

# TODO: make a molecule class to replace dictionary

def collect_molecules(library_dir, defaults=True):
    if defaults:
        library_dir = resource_filename('analysis',
                'molecules/{}/'.format(library_dir))

    library = load_jsons(library_dir)

    return library

def load_json(fname):
    with open(fname, "r") as f:
        molecule_dict = json.load(f)
        return {molecule_dict["name"] : molecule_dict}


def load_jsons(library_dir):
    json_files = glob.glob(os.path.join(library_dir, "*.json"))
    library = dict()
    for fname in json_files:
        library.update(load_json(fname))
    return library
