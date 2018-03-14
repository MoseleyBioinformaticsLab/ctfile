from __future__ import print_function, division, unicode_literals
from collections import deque
from collections import namedtuple


HeaderBlock = namedtuple('HeaderBlock', ['molecule_name', 'software', 'comment'])

CtabBlock = namedtuple('CtabBlock', [])

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
    for entry in text.split('$$$$'):
        if entry:
            lines_stream = deque(entry.split('\n'))
        else:
            continue

        yield HeaderBlock(lines_stream.popleft().strip(), lines_stream.popleft().strip(), lines_stream.popleft().strip())

        counts_line = lines_stream.popleft()
        counts_line_values = [counts_line[i:i + 3].strip() for i in range(0, len(counts_line) - 6, 3)] + [counts_line[-6:len(counts_line)].strip()]
        ctab_counts_line = CtabCountsLine(*counts_line_values)
        number_of_atoms = ctab_counts_line.number_of_atoms
        number_of_bonds = ctab_counts_line.number_of_bonds

        yield CtabBlock()
        yield ctab_counts_line

        yield from _atom_bond_block(number_of_lines=number_of_atoms, block_type=CtabAtomBlockLine, stream=lines_stream)
        yield from _atom_bond_block(number_of_lines=number_of_bonds, block_type=CtabBondBlockLine, stream=lines_stream)
        yield from _property_block(stream=lines_stream)
        yield from _data_block(stream=lines_stream)


def _atom_bond_block(number_of_lines, block_type, stream):
    """

    :param number_of_lines:
    :param block_type:
    :param stream:
    :return:
    """
    for _ in range(int(number_of_lines)):
        line = stream.popleft()
        yield block_type(*line.split())


def _property_block(stream):
    """

    :param stream:
    :return:
    """
    line = ''
    while line != 'M  END':
        line = stream.popleft()
        name = line.split()[1]

        if name == 'END':
            continue
        else:
            yield CtabPropertiesBlockLine(name, line)

def _data_block(stream):
    """
    
    :param stream: 
    :return: 
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
