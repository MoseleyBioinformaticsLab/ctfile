#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ctfile.tokenizer
~~~~~~~~~~~~~~~~

This module implements tokenizer that processes text file and yields 
tokens necessary to create ``CTfile`` objects.
"""

from __future__ import print_function, division, unicode_literals
from collections import deque
from collections import namedtuple


HeaderBlock = namedtuple('HeaderBlock', ['molecule_name', 'software', 'comment'])
MolfileStart = namedtuple('MolfileStart', [])
MolfileEnd = namedtuple('MolfileEnd', [])
CtabBlockStart = namedtuple('CtabBlockStart', [])
CtabBlockEnd = namedtuple('CtabBlockEnd', [])
DataBlockStart = namedtuple('DataBlockStart', [])
DataBlockEnd = namedtuple('DataBlockEnd', [])
EndOfFile = namedtuple('EndOfFile', [])

CtabCountsLine = namedtuple('CtabCountsLine', ['number_of_atoms', 'number_of_bonds', 'number_of_atom_lists',
                                               'not_used1', 'chiral_flag', 'number_of_stext_entries', 'not_used2',
                                               'not_used3', 'not_used4', 'not_used5', 'number_of_properties',
                                               'version'])

CtabAtomBlockLine = namedtuple('CtabAtomBlock', ['x', 'y', 'z', 'atom_symbol', 'mass_difference', 'charge',
                                                 'atom_stereo_parity', 'hydrogen_count', 'stereo_care_box',
                                                 'valence', 'h0designator', 'not_used1', 'not_used2',
                                                 'atom_atom_mapping_number', 'inversion_retention_flag',
                                                 'exact_change_flag'])

CtabBondBlockLine = namedtuple('CtabBondBlock', ['first_atom_number', 'second_atom_number', 'bond_type',
                                                 'bond_stereo', 'not_used1', 'bond_topology',
                                                 'reacting_center_status'])

CtabPropertiesBlockLine = namedtuple('CtabPropertiesBlock', ['name', 'line'])


DataHeader = namedtuple('DataHeader', ['header'])
DataItem = namedtuple('DataItem', ['data_item'])


def tokenizer(text):
    """A lexical analyzer for the `CTfile` formatted files.

    :param str text: `CTfile` formatted text.
    :return: Tuples of data.
    :rtype: py:class:`~collections.namedtuple`
    """
    for entry in text.split('$$$$\n'):
        if entry.rstrip():
            lines_stream = deque(entry.split('\n'))
        else:
            continue

        # yield from _molfile(stream=lines_stream)
        for token in _molfile(stream=lines_stream):
            yield token

        if len(lines_stream):
            # yield from _sdfile(stream=lines_stream)
            for token in _sdfile(stream=lines_stream):
                yield token

    yield EndOfFile()


def _molfile(stream):
    """Process ``Molfile``.
    
    :param stream: Queue containing lines of text.
    :type stream: :py:class:`collections.deque`
    :return: Tuples of data.
    """
    yield MolfileStart()
    yield HeaderBlock(stream.popleft().strip(), stream.popleft().strip(), stream.popleft().strip())
    # yield from _ctab(stream)
    for token in _ctab(stream):
        yield token
    yield MolfileEnd()


def _sdfile(stream):
    """Process ``SDfile``.
    
    :param stream: Queue containing lines of text.
    :type stream: :py:class:`collections.deque`
    :return: Tuples of data.
    """
    yield DataBlockStart()
    # yield from _data_block(stream=stream)
    for token in _data_block(stream=stream):
        yield token
    yield DataBlockEnd()


def _ctab(stream):
    """Process ``Ctab``.
    
    :param stream: Queue containing lines of text.
    :type stream: :py:class:`collections.deque`
    :return: Tuples of data.
    """
    yield CtabBlockStart()
    counts_line = stream.popleft()
    counts_line_values = [counts_line[i:i + 3].strip() for i in range(0, len(counts_line) - 6, 3)] + \
                         [counts_line[-6:len(counts_line)].strip()]
    ctab_counts_line = CtabCountsLine(*counts_line_values)
    yield ctab_counts_line

    number_of_atoms = ctab_counts_line.number_of_atoms
    number_of_bonds = ctab_counts_line.number_of_bonds

    # yield from _ctab_atom_bond_block(number_of_lines=number_of_atoms, block_type=CtabAtomBlockLine, stream=stream)
    for token in _ctab_atom_bond_block(number_of_lines=number_of_atoms, block_type=CtabAtomBlockLine, stream=stream):
        yield token

    # yield from _ctab_atom_bond_block(number_of_lines=number_of_bonds, block_type=CtabBondBlockLine, stream=stream)
    for token in _ctab_atom_bond_block(number_of_lines=number_of_bonds, block_type=CtabBondBlockLine, stream=stream):
        yield token

    # yield from _ctab_property_block(stream=stream)
    for token in _ctab_property_block(stream=stream):
        yield token

    yield CtabBlockEnd()


def _ctab_atom_bond_block(number_of_lines, block_type, stream):
    """Process atom and bond blocks of ``Ctab``.

    :param number_of_lines: Number of lines to process from stream.
    :param block_type: :py:class:`collections.namedtuple` to use for data processing.
    :type block_type: :class:`~ctfile.tokenizer.CtabAtomBlockLine` or :class:`~ctfile.tokenizer.CtabBondBlockLine`
    :param stream: Queue containing lines of text.
    :type stream: :py:class:`collections.deque`
    :return: Tuples of data.
    :rtype: :class:`~ctfile.tokenizer.CtabAtomBlockLine` or :class:`~ctfile.tokenizer.CtabBondBlockLine`
    """
    for _ in range(int(number_of_lines)):
        line = stream.popleft()
        yield block_type(*line.split())


def _ctab_property_block(stream):
    """Process properties block of ``Ctab``.

    :param stream: Queue containing lines of text.
    :type stream: :py:class:`collections.deque`
    :return: Tuples of data.
    :rtype: :class:`~ctfile.tokenizer.CtabPropertiesBlockLine`
    """
    line = stream.popleft()
    while line != 'M  END':
        name = line.split()[1]
        yield CtabPropertiesBlockLine(name, line)
        line = stream.popleft()


def _data_block(stream):
    """Process data block of ``CTfile``.
    
    :param stream: Queue containing lines of text.
    :type stream: :py:class:`collections.deque`
    :return: Tuples of data.
    :rtype: :class:`~ctfile.tokenizer.DataHeader` or :class:`~ctfile.tokenizer.DataItem`
    """
    while len(stream) > 0:
        line = stream.popleft()

        if line.startswith('>'):
            yield DataHeader(line[1:].strip())
        else:
            data_item = line.strip()
            if data_item:
                yield DataItem(line)
            else:
                continue
