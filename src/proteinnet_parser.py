"""
Text-based parser for ProteinNet Records.
"""

__author__ = "Mohammed AlQuraishi"
__copyright__ = "Copyright 2019, Harvard Medical School"
__license__ = "MIT"

# !/usr/bin/python

# imports
import sys
import re
import numpy as np
from itertools import groupby

# Constants
NUM_DIMENSIONS = 3

# Functions for conversion from Mathematica protein files to TFRecords
_aa_dict = {'A': '0', 'C': '1', 'D': '2', 'E': '3', 'F': '4', 'G': '5', 'H': '6', 'I': '7', 'K': '8', 'L': '9', 'M': '10', 'N': '11', 'P': '12', 'Q': '13', 'R': '14', 'S': '15', 'T': '16', 'V': '17', 'W': '18', 'Y': '19'}
_dssp_dict = {'L': '0', 'H': '1', 'B': '2', 'E': '3', 'G': '4', 'I': '5', 'T': '6', 'S': '7'}
_mask_dict = {'-': '0', '+': '1'}


def letter_to_num(string, dict_):
    """ Convert string of letters to list of ints """
    patt = re.compile('[' + ''.join(dict_.keys()) + ']')
    num_string = patt.sub(lambda m: dict_[m.group(0)] + ' ', string)
    return [int(i) for i in num_string.split()]


def yield_records_from_file(file, num_evo_entries: int = 20):
    def get_record(lines):
        entry = {"ID": lines[0].strip()}
        for i, line in enumerate(lines):
            if line == '[PRIMARY]' + '\n':
                primary = lines[i + 1].strip()
                entry.update({'primary': primary})
            elif line == '[EVOLUTIONARY]' + '\n':
                evolutionary = []
                for residue in range(num_evo_entries):
                    evolutionary.append(
                        [float(step) for step in lines[i + 1].strip().split()]
                    )
                entry.update({'evolutionary': np.array(evolutionary)})
            elif line == '[SECONDARY]' + '\n':
                secondary = letter_to_num(lines[i + 1].strip(), _dssp_dict)
                entry.update({'secondary': secondary})
            elif line == '[TERTIARY]' + '\n':
                tertiary = []
                for axis in range(NUM_DIMENSIONS):
                    tertiary.append([float(coord) for coord in lines[i + 1 + axis].strip().split()])
                entry.update({'tertiary': np.array(tertiary).T})
            elif line == '[MASK]' + '\n':
                mask = letter_to_num(lines[i + 1].strip(), _mask_dict)
                entry.update({'mask': mask})
            else:
                continue
        return entry

    for k, g in groupby(open(file, "r"), lambda x: x.startswith("[ID]")):
        if not k:
            yield get_record(list(g))


def clean_entry(entry, atom="ca"):
    sequence = "primary"
    mask = np.where(np.array(entry['mask']) == 1)[0]
    entry[sequence] = ''.join(entry[sequence][x] for x in mask)
    mask_3d = np.array([i for n in mask for i in range(n*3, n*3+3)]).astype(int)
    entry['tertiary'] = entry['tertiary'][mask_3d]
    if atom == "ca":
        index = 1
    elif atom == "n":
        index = 0
    elif atom == "cb":
        index = 2
    else:
        raise ValueError("atom must be one of n, ca, cb")
    entry['tertiary'] = entry['tertiary'][np.arange(index, entry['tertiary'].shape[0]+index, 3)] / 100
    assert entry['tertiary'].shape[0] == len(entry[sequence]), (entry['tertiary'].shape[0], len(entry[sequence]))
    return entry
