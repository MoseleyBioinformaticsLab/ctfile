#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

"""
ctfile.utils
~~~~~~~~~~~~

This module provides reusable utility functions and objects.
"""

from collections import OrderedDict
from collections import Counter
from collections import defaultdict
import ctfile


def d_colorize_mol(mol, max_d, isotope_resolved=True, stereo_resolved=True):
    """

    This is an improved implementation of the coloring from CASS, see doi:10.3389/fgene.2014.00237

    Every node in a graph can be assigned one or more 'd' colors which represent the local substructure
    around that node out to d edges away. The d color of a node can be expressed as a combination of the
    d-1 colors of its neighbors. Every node has a trivial d0 color that is its element type and charge.
    For example, an uncharged carbon has the d0-color C0

    If isotope resolved, the isotope symbol is added to the d0 color. e.g. 13C0 (only non-monoisotopic
    isotopes are added to the color).

    The dN color of a node is equal to its: d0 color + the dN-1 color of every neighbor + the order of the
    bond connecting that neighbor, which parantheses separating different components of the color.

    If stereo resolved, the stereo identifier of the bond is added to the bond order.

    Atoms with identical d3 colors are centered in the same substructure out to three bonds. They will
    have similar / identical J-coupling.

    This function will add a 'color' field to every atom in a mol file object and a 'color_groups' field
    to the mol file. The color_groups represent the sets of atoms with identical dN colors.

    :param mol: a molfile object from ctfile
    :param max_d: integer - the maximum d out to which to color
    :param isotope_resolved: boolean - if true, add non-monoisotope information to d0 colors
    :param stereo_resolved: boolean - if true, add stereo information to bonds when constructing colors.
    """
    if not isinstance(mol, ctfile.Molfile):
        raise Exception("can only colorize mol files")

    # make the isotope lookup table by atom index if there are isotopes and if the flag is set.
    if 'ISO' in mol['Ctab']['CtabPropertiesBlock'] and isotope_resolved:
        isotope_lookup = {entry['atom_number'] : entry['absolute_mass'] for entry in mol['Ctab']['CtabPropertiesBlock']['ISO']}
    else:
        isotope_lookup = {}

    bond_lookup = {}
    for bond in mol.bonds:
        bond_lookup[(bond['first_atom_number'], bond['second_atom_number'])] = bond
        bond_lookup[(bond['second_atom_number'], bond['first_atom_number'])] = bond

    # make the base d0 color for every atom
    for index, atom in enumerate(mol.atoms):
        d0_color = isotope_lookup[index] + atom["atom_symbol"] + atom["charge"] if index in isotope_lookup else atom["atom_symbol"] + atom["charge"]
        atom["colors"] = {0: d0_color}

    # now make the dN color for N = 1 to max_d
    for d in range(1, max_d):
        for n, atom in enumerate(mol.atoms):
            color_components = []
            for neighbor in atom.neighbors:
                connecting_bond = bond_lookup[(neighbor.atom_number, atom.atom_number)]
                if stereo_resolved:
                    color_components.append((neighbor["colors"][d - 1], connecting_bond['bond_type'] + connecting_bond['bond_stereo']))
                else:
                    color_components.append((neighbor["colors"][d - 1], connecting_bond['bond_type']))
            d_color = atom["colors"][0]
            for a, b in sorted(color_components):
                d_color += '(' + a + "," + b + ')'
            atom["colors"][d] = d_color

    # group all atoms by their colors for N = 0 to N = max_d.
    for d in range(0, max_d):
        color_groups = defaultdict(dict)
        for atom in mol.atoms:
            d_color = atom["colors"][d]
            if d_color not in color_groups[d]:
                color_groups[d][d_color] = []
            color_groups[d][d_color].append(atom)
    mol["color_groups"] = color_groups
    return mol

class OrderedCounter(Counter, OrderedDict):
    """Counter that remembers the order of elements that are seen first."""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, OrderedDict(self)

