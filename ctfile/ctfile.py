#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

"""
ctfile.ctfile
~~~~~~~~~~~~~

This module implements core objects to represent ``CTfile`` file format.
"""

from __future__ import print_function, division, unicode_literals
import sys
import json
import io
from collections import OrderedDict

import more_itertools

from .tokenizer import tokenizer
from .conf import ctab_properties_conf
from .utils import OrderedCounter


class CTfile(OrderedDict):
    """Base class to represent collection of Chemical table file (``CTfile``) 
    formats, e.g. ``Molfile``, ``SDfile``."""

    ctab_properties = ctab_properties_conf

    def __init__(self, *args, **kwargs):
        """CTfile initializer.
        
        :param lexer: instance of the ``CTfile`` format tokenizer.
        :type lexer: :func:`~ctfile.tokenizer.tokenizer`.
        """
        super(CTfile, self).__init__(*args, **kwargs)

    @classmethod
    def load(cls, filehandle):
        """Load data into ``CTfile`` object from file-like object.

        :param filehandle: File-like object.
        :return: Instance of ``CTfile``.
        """
        return cls.loadstr(filehandle.read())

    @classmethod
    def loadstr(cls, string):
        """Load data into ``CTfile`` object from string.
        
        :param str string: String containing data in ``CTfile`` format.
        :return: Instance of ``CTfile``.
        """
        try:
            string = string.decode('utf-8')
        except AttributeError:
            pass

        lexer = tokenizer(string)

        if cls.is_molfile(string):
            molfile = Molfile()

            try:
                molfile.update(json.loads(string))
            except ValueError:
                molfile._build(lexer)
            return molfile

        elif cls.is_sdfile(string):
            sdfile = SDfile()

            try:
                sdfile.update(json.loads(string))
            except ValueError:
                sdfile._build(lexer)
            return sdfile

        else:
            raise ValueError('Cannot determine the format of string.')

    def write(self, filehandle, file_format):
        """Write :class:`~ctfile.ctfile.CTfile` data into file. 

        :param filehandle: File-like object.
        :param str file_format: Format to use to write data: ``ctfile`` or ``json``.
        :return: None.
        :rtype: :py:obj:`None`.
        """
        try:
            filehandle.write(self.writestr(file_format=file_format))
        except IOError:
            raise IOError('"filehandle" parameter must be writable.')

    def writestr(self, file_format):
        """Write :class:`~ctfile.ctfile.CTfile` data into string.
        
        :param str file_format: Format to use to write data: ``ctfile`` or ``json``.
        :return: String representing the :class:`~ctfile.ctfile.CTfile` instance.
        :rtype: :py:class:`str`.
        """
        if file_format == 'json':
            repr_str = self._to_json()
        elif file_format == 'ctfile':
            repr_str = self._to_ctfile()
        else:
            raise ValueError('Invalid "file_format" argument: "{}"'.format(file_format))
        return repr_str

    def print_file(self, file_format='ctfile', f=sys.stdout):
        """Print representation of :class:`~ctfile.ctfile.CTfile`.

        :param str file_format: Format to use: ``ctfile`` or ``json``.
        :param f: Print to file or stdout.
        :type f: File-like 
        :return: None.
        :rtype: :py:obj:`None`.
        """
        print(self.writestr(file_format=file_format), file=f)

    def _build(self, lexer):
        """Build :class:`~ctfile.ctfile.CTfile` instance.

        :return: :class:`~ctfile.ctfile.CTfile` instance.
        :rtype: :class:`~ctfile.ctfile.CTfile`.
        """
        raise NotImplementedError('Subclass must implement abstract method')

    def _to_json(self, sort_keys=False, indent=4):
        """Convert :class:`~ctfile.ctfile.CTfile` into JSON string.

        :return: ``JSON`` formatted string.
        :rtype: :py:class:`str`.
        """
        return json.dumps(self, sort_keys=sort_keys, indent=indent, cls=CtabAtomBondEncoder)

    def _to_ctfile(self):
        """Convert :class:`~ctfile.ctfile.CTfile` into `CTfile` formatted string.
        
        :return: ``CTfile`` formatted string.
        :rtype: :py:class:`str`.
        """
        raise NotImplementedError('Subclass must implement abstract method')

    @staticmethod
    def is_molfile(string):
        """Test if input string is in ``Molfile`` format.

        :param string: Input string.
        :type string: :py:class:`str` or :py:class:`bytes`.
        :return: Input string if in ``Molfile`` format or False otherwise.
        :rtype: :py:class:`str` or :py:obj:`False`.
        """
        try:
            string = string.decode('utf-8')
        except AttributeError:
            pass

        if '$$$$\n' in string:
            return False
        return True

    @staticmethod
    def is_sdfile(string):
        """Test if input string is in ``SDfile`` format.

        :param string: Input string.
        :type string: :py:class:`str` or :py:class:`bytes`.
        :return: Input string if in ``SDfile`` format or False otherwise.
        :rtype: :py:class:`str` or :py:obj:`False`.
        """
        try:
            string = string.decode('utf-8')
        except AttributeError:
            pass

        if '$$$$\n' in string:
            return True
        return False


class Ctab(CTfile):
    """Ctab - connection table contains information describing the structural relationships
    and properties of a collection of atoms.
    
    --------------------
    | CTab             |
    |                  |
    | Counts line      |
    | Atom block       |
    | Bond block       |
    | Properties block |
    |                  |
    --------------------
    
    * Counts line: specifies the number of atoms, bonds, Sgroups, 3D constituents, as well as
      the chiral flag setting, and the regno.
    * Atom block: specifies an atomic symbol and any mass difference, charge, stereochemistry,
      and associated hydrogens for each atom.
    * Bond block: Specifies the two atoms connected by the bond, the bond type, and any bond
      stereochemistry and topology (chain or ring properties) for each bond.
    * Properties block: specifies additional properties.
    
    counts line format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    where:
        aaa = number of atoms
        bbb = number of bonds
        lll = number of atom lists
        fff = (obsolete)
        ccc = chiral flag: 0=not chiral, 1=chiral
        sss = number of stext entries
        xxx = (obsolete)
        rrr = (obsolete)
        ppp = (obsolete)
        iii = (obsolete)
        mmm = number of lines of additional properties
    
    atom block format: xxxxxxxxxxyyyyyyyyyyzzzzzzzzzzaaaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    where:
        xxxxxxxxxx = atom x coordinate
        yyyyyyyyyy = atom y coordinate
        zzzzzzzzzz = atom z coordinate
        aaa        = atom symbol
        dd         = mass difference: -3, -2, -1, 0, 1, 2, 3, 4, (0 if value beyond these limits)
        ccc        = charge: 0=uncharged or value other than these, 1=+3, 2=+2, 3=+1, 
                     4=doublet radical, 5=-1, 6=-2, 7=-3
        sss        = atom stereo parity: 0=not stereo, 1=odd, 2=even, 3=either or unmarked stereo center
        hhh        = hydrogen count + 1: 1=H0, 2=H1, 3=H2, 4=H3, 5=H4
        bbb        = stereo care box: 0=ignore stereo configuration of this double bond atom, 
                     1=stereo configuration of double bond atom must match
        vvv        = valence: 0=no marking (default), (1 to 14)=(1 to 14), 15=zero valence
        HHH        = H0 designator: 0=not specified, 1=no H atoms allowed
        rrr        = (obsolete)
        iii        = (obsolete)
        mmm        = atom-atom mapping number: 1=number of atoms
        nnn        = inversion/retention flag: 0=property not applied 1=configuration is inverted, 
                     2=configuration is retained
        eee        = exact change flag: 0=property not applied, 1=change on atom must be exactly as shown
    
    bond block format: 111222tttsssxxxrrrccc
    where:
        111 = first atom number: 1=number of atoms
        222 = second atom number: 1=number of atoms
        ttt = bond type: 1=Single, 2=Double, 3=Triple, 4=Aromatic, 5=Single or Double, 6=Single or Aromatic, 
              7=Double or Aromatic, 8=Any
        sss = bond stereo: Single bonds: 0=not stereo, 1=Up, 4=Either, 6=Down; 
              Double bonds: 0=Use x-, y-, z-coords from atom block to determine cis or trans, 
              3=Cis or trans (either) double bond
        xxx = (obsolete)
        rrr = bond topology: 0=Either, 1=Ring, 2=Chain
        ccc = reacting center status: 0=unmarked, 1=a center, -1=not a center; 
              Additional: 2=no change, 4=bond made/broken, 8=bond order changes 12=4+8 (both made/broken and changes); 
              5=(4 + 1), 9=(8 + 1), and 13=(12 + 1) are also possible
    
    properties block: 
    where:
        Most lines in the properties block are identified by a prefix of the form "M  XXX" where two spaces 
        separate the M and XXX.
        The prefix: "M  END" terminates the properties block.
    """
    counts_line_format = 'aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv'
    atom_block_format = 'xxxxxxxxxxyyyyyyyyyyzzzzzzzzzzaaaaddcccssshhhbbbvvvHHHrrriiimmmnnneee'
    bond_block_format = '111222tttsssxxxrrrccc'

    def __init__(self):
        """Ctab initializer."""
        super(Ctab, self).__init__()
        self['CtabCountsLine'] = OrderedDict()
        self['CtabAtomBlock'] = []
        self['CtabBondBlock'] = []
        self['CtabPropertiesBlock'] = OrderedDict()

        self.edges = OrderedDict()

    def _build(self, lexer):
        """Build :class:`~ctfile.ctfile.Ctab` instance.
        
        :return: :class:`~ctfile.ctfile.Ctab` instance.
        :rtype: :class:`~ctfile.ctfile.Ctab`.
        """
        atom_number = 1
        while True:
            token = next(lexer)
            key = token.__class__.__name__

            if key == 'CtabCountsLine':
                self[key].update(token._asdict())

            elif key == 'CtabAtomBlock':
                self[key].append(Atom(atom_number=str(atom_number), **token._asdict()))
                atom_number += 1

            elif key == 'CtabBondBlock':
                first_atom_number, second_atom_number, bond_type, bond_stereo, \
                not_used1, bond_topology, reacting_center_status = token

                first_atom = self['CtabAtomBlock'][int(first_atom_number) - 1]
                second_atom = self['CtabAtomBlock'][int(second_atom_number) - 1]

                first_atom.neighbors.append(second_atom)
                second_atom.neighbors.append(first_atom)

                bond = Bond(first_atom=first_atom, second_atom=second_atom, bond_type=bond_type,
                            bond_stereo=bond_stereo, not_used1=not_used1, bond_topology=bond_topology,
                            reacting_center_status=reacting_center_status)
                self[key].append(bond)

            elif key == 'CtabPropertiesBlock':
                property_name = token.name
                single_entry_keys = self.ctab_properties[self.version][property_name]['values']
                self[key].setdefault(property_name, [])
                ctab_properties = more_itertools.sliced(token.line.split()[3:], len(single_entry_keys))

                for ctab_property in ctab_properties:
                    self[key][property_name].append(OrderedDict(zip(single_entry_keys, ctab_property)))

            elif key == 'CtabBlockEnd':
                break

            else:
                raise KeyError('Ctab object does not supposed to have any other information: "{}".'.format(key))

    def _to_ctfile(self):
        """Convert :class:`~ctfile.ctfile.CTfile` into `CTfile` formatted string.

        :return: `CTfile` formatted string.
        :rtype: :py:class:`str`.
        """
        output = io.StringIO()

        for key in self:

            if key == 'CtabCountsLine':
                counter = OrderedCounter(self.counts_line_format)
                counts_line = ''.join([str(value).rjust(spacing) for value, spacing
                                       in zip(self[key].values(), counter.values())])
                output.write(counts_line)
                output.write('\n')

            elif key == 'CtabAtomBlock':
                counter = OrderedCounter(Atom.atom_block_format)
                for atom in self[key]:
                    atom_line = ''.join([str(value).rjust(spacing) for value, spacing
                                         in zip(atom._ctab_data.values(), counter.values())])
                    output.write(atom_line)
                    output.write('\n')

            elif key == 'CtabBondBlock':
                counter = OrderedCounter(Bond.bond_block_format)
                for bond in self[key]:
                    bond_line = ''.join([str(value).rjust(spacing) for value, spacing
                                         in zip(bond._ctab_data.values(), counter.values())])
                    output.write(bond_line)
                    output.write('\n')

            elif key == 'CtabPropertiesBlock':
                for property_name in self[key]:
                    ctab_property_identifier = self.ctab_properties[self.version][property_name]['fmt']

                    for entry in self[key][property_name]:
                        ctab_property_line = '{}  {}{}'.format(ctab_property_identifier, 1, ''.join([str(value).rjust(4) for value in entry.values()]))
                        output.write('{}'.format(ctab_property_line))
                        output.write('\n')
                output.write('{}'.format(self.ctab_properties[self.version]['END']['fmt']))
                output.write('\n')

            else:
                raise KeyError('Ctab object does not supposed to have any other information: "{}".'.format(key))

        return output.getvalue()

    @property
    def version(self):
        """Version of the `CTfile` formatting.
        
        :return: Version of the `CTfile`.
        :rtype: :py:class:`str`.
        """
        return self['CtabCountsLine']['version']

    @property
    def atoms(self):
        """List of atoms.

        :return: List of atoms.
        :rtype: :py:class:`list`.
        """
        return self['CtabAtomBlock']

    @property
    def bonds(self):
        """List of atoms.

        :return: List of atoms.
        :rtype: :py:class:`list`.
        """
        return self['CtabBondBlock']

    @property
    def positions(self):
        """List of positions of atoms in atoms block starting from 1.
        
        :return: List of positions of atoms in atoms block starting from 1. 
        :rtype: :py:class:`list`.
        """
        return [str(i) for i in range(1, len(self.atoms)+1)]

    @property
    def iso(self, property_specifier='ISO'):
        """Return list of isotopic properties per each atom position.
        
        :return: List of isotopic properties per each atom position.
        :rtype: :py:class:`list`.
        """
        isotopes = []

        if property_specifier in self['CtabPropertiesBlock']:
            atoms_by_position = dict(zip(self.positions, self.atoms))

            for iso_property in self['CtabPropertiesBlock'][property_specifier]:
                atom = atoms_by_position[iso_property['atom_number']]
                isotopes.append({'atom_symbol': atom['atom_symbol'],
                                 'isotope': iso_property['absolute_mass'],
                                 'atom_number': iso_property['atom_number']})
        return isotopes

    def atoms_by_symbol(self, atom_symbol):
        """Access all atoms of specified type.

        :param str atom_symbol: Atom symbol.
        :return: List of atoms.
        :rtype: :py:class:`list`.
        """
        return [atom for atom in self['CtabAtomBlock'] if atom['atom_symbol'] == atom_symbol]

    @property
    def carbon_atoms(self, atom_symbol='C'):
        """Access all carbon atoms within ``Ctab`` atom block.
        
        :param str atom_symbol: Carbon atom symbol.
        :return: List of carbon atoms.
        :rtype: :py:class:`list`.
        """
        return self.atoms_by_symbol(atom_symbol=atom_symbol)

    @property
    def hydrogen_atoms(self, atom_symbol='H'):
        """Access all hydrogen atoms within ``Ctab`` atom block.

        :param str atom_symbol: Hydrogen atom symbol.
        :return: List of hydrogen atoms.
        :rtype: :py:class:`list`.
        """
        return self.atoms_by_symbol(atom_symbol=atom_symbol)

    def add_ctab_property(self, ctab_property_name, values):
        """Add new values to existing ``CtabPropertiesBlock``.

        :param str ctab_property_name: Name of the ``Ctab`` property.
        :param values: Sequence of values. 
        :return: None.
        :rtype: :py:obj:`None`.
        """
        if ctab_property_name.upper() not in ctab_properties_conf[self.version]:
            raise ValueError('Unknown property: "{}".\n'
                             'Available Ctab properties are: {}'.format(ctab_property_name, ", ".join('"{}"'.format(ctab_property) for ctab_property in ctab_properties_conf[self.version])))
        else:
            self['CtabPropertiesBlock'].setdefault(ctab_property_name, [])

            single_entry_keys = self.ctab_properties[self.version][ctab_property_name]['values']
            for value in values:
                annotated_value = OrderedDict(zip(single_entry_keys, value))
                if annotated_value not in self['CtabPropertiesBlock'][ctab_property_name]:
                    self['CtabPropertiesBlock'][ctab_property_name].append(annotated_value)

    def replace_ctab_property(self, ctab_property_name, values):
        """Replace with new values ``CtabPropertiesBlock``.

        :param str ctab_property_name: Name of the ``Ctab`` property.
        :param values: Sequence of values. 
        :return: None.
        :rtype: :py:obj:`None`.
        """
        if ctab_property_name.upper() not in ctab_properties_conf[self.version]:
            raise ValueError('Unknown property: "{}".\n'
                             'Available Ctab properties are: {}'.format(ctab_property_name, ', '.join('"{}"'.format(ctab_property) for ctab_property in ctab_properties_conf[self.version])))
        else:
            self['CtabPropertiesBlock'][ctab_property_name] = []
            single_entry_keys = self.ctab_properties[self.version][ctab_property_name]['values']
            for value in values:
                annotated_value = OrderedDict(zip(single_entry_keys, value))
                if annotated_value not in self['CtabPropertiesBlock'][ctab_property_name]:
                    self['CtabPropertiesBlock'][ctab_property_name].append(annotated_value)

    def add_charge(self, atom_number, atom_symbol, charge):
        """Add charge to neutral "N", "O", or "S".
        
        :param str atom_number: Atom id in order of appearance in ``CTfile``.
        :param str atom_symbol: Atom symbol.
        :param str charge: Charge: "+1" or "-1".
        :return: None.
        :rtype: :py:obj:`None`.
        """
        atom = self['CtabAtomBlock'][int(atom_number)-1]

        if atom.atom_symbol != atom_symbol:
            raise ValueError('Mismatch between atom symbol and atom number.')

        if atom.atom_symbol not in {'N', 'O', 'S'}:
            raise ValueError('Cannot ionize atom: "{}"'.format(atom.atom_symbol))

        if int(charge) == 1:
            hydrogen = Atom(atom_number=str(len(self['CtabAtomBlock'])+1),
                            atom_symbol='H',
                            x='0.0000',
                            y='0.0000',
                            z='0.0000',
                            mass_difference='0',
                            charge='0',
                            atom_stereo_parity='0',
                            hydrogen_count='0',
                            stereo_care_box='0',
                            valence='0',
                            h0designator='0',
                            not_used1='0',
                            not_used2='0',
                            atom_atom_mapping_number='0',
                            inversion_retention_flag='0',
                            exact_change_flag='0')

            bond = Bond(first_atom=atom,
                        second_atom=hydrogen,
                        bond_type='1',
                        bond_stereo='0',
                        not_used1='0',
                        bond_topology='0',
                        reacting_center_status='0')

            self['CtabAtomBlock'].append(hydrogen)
            self['CtabBondBlock'].append(bond)

        elif int(charge) == -1:
            terminal_hydrogens = atom.neighbor_hydrogen_atoms

            if len(terminal_hydrogens) == 1:
                hydrogen = terminal_hydrogens[0]
                remove_index = []

                # update atom numbers
                for atom in self['CtabAtomBlock']:
                    if int(atom.atom_number) > int(hydrogen.atom_number):
                        atom.atom_number = str(int(atom.atom_number)-1)

                # find index of a bond to remove and update ctab data dict with new atom numbers
                for index, bond in enumerate(self['CtabBondBlock']):
                    if hydrogen.atom_number in {bond.first_atom_number, bond.second_atom_number}:
                        remove_index.append(index)
                    bond.update_atom_numbers()

                if len(remove_index) != 1:
                    raise ValueError('Removing more than one bond! Operation is not permitted.')
                else:
                    # remove atom from neighbors list
                    for atom in self['CtabAtomBlock']:
                        if hydrogen in atom.neighbors:
                            atom.neighbors.remove(hydrogen)

                    # remove atom and bond from Ctab
                    self['CtabAtomBlock'].pop(int(hydrogen.atom_number)-1)
                    self['CtabBondBlock'].pop(remove_index[0])
            else:
                raise ValueError('Atom "{}{}" has incorrect number of terminal hydrogens.'.format(atom.atom_symbol, atom.atom_number))

        else:
            raise ValueError('Can only handle charges "+1" and "-1".')

        # update atom and bond counts
        self['CtabCountsLine']['number_of_atoms'] = str(len(self['CtabAtomBlock']))
        self['CtabCountsLine']['number_of_bonds'] = str(len(self['CtabBondBlock']))
        self.add_ctab_property(ctab_property_name='CHG', values=[(str(atom_number), str(charge))])


class Molfile(CTfile):
    """Molfile - each molfile describes a single molecular structure which can
    contain disjoint fragments.
    
    --------------------
    | molfile          |
    |                  |
    |                  |
    |   ---------      |
    |   | Ctab  |      |
    |   ---------      |
    |                  |
    --------------------
    """
    def __init__(self):
        """Molfile initializer."""

        super(Molfile, self).__init__()
        self['HeaderBlock'] = OrderedDict()
        self['Ctab'] = OrderedDict()

    def _build(self, lexer):
        """Build :class:`~ctfile.ctfile.Molfile` instance.

        :return: :class:`~ctfile.ctfile.Molfile` instance.
        :rtype: :class:`~ctfile.ctfile.Molfile`.
        """
        key = ''
        while key != 'EndOfFile':

            token = next(lexer)
            key = token.__class__.__name__

            if key == 'HeaderBlock':
                self[key].update(token._asdict())

            elif key == 'CtabBlockStart':
                ctab = Ctab()
                ctab._build(lexer)
                self['Ctab'] = ctab

            elif key == 'MolfileStart':
                pass

            elif key in ('MolfileEnd', 'EndOfFile'):
                break

            else:
                raise KeyError('Molfile object does not supposed to have any other information: "{}".'.format(key))

        return self

    def _to_ctfile(self):
        """Convert :class:`~ctfile.ctfile.CTfile` into `CTfile` formatted string.

        :return: ``CTfile`` formatted string.
        :rtype: :py:class:`str`.
        """
        output = io.StringIO()

        for key in self:
            if key == 'HeaderBlock':
                for line in self[key].values():
                    output.write(line)
                    output.write('\n')

            elif key == 'Ctab':
                ctab_str = self[key]._to_ctfile()
                output.write(ctab_str)

            else:
                raise KeyError('Molfile object does not supposed to have any other information: "{}".'.format(key))

        return output.getvalue()

    @property
    def version(self):
        """Version of the `CTfile` formatting.

        :return: Version of the `CTfile`.
        :rtype: :py:class:`str`.
        """
        return self['Ctab'].version

    @property
    def atoms(self):
        """List of atoms.

        :return: List of atoms.
        :rtype: :py:class:`list`.
        """
        return self['Ctab'].atoms

    @property
    def bonds(self):
        """List of bonds.

        :return: List of bonds.
        :rtype: :py:class:`list`.
        """
        return self['Ctab'].bonds

    @property
    def positions(self):
        """List of positions of atoms in atoms block starting from 1.
        
        :return: List of positions of atoms in atoms block starting from 1. 
        :rtype: :py:class:`list`.
        """
        return self['Ctab'].positions

    @property
    def iso(self):
        """Return list of isotopic properties per each atom position.
        
        :return: List of isotopic properties per each atom position.
        :rtype: :py:class:`list`.
        """
        return self['Ctab'].iso

    @property
    def molfiles(self):
        """Create list of ``Molfile`` instances (list of 1 in this case).

        :return: List of ``Molfile`` instances.
        :rtype: :py:class:`list`.
        """
        return [self]

    @property
    def carbon_atoms(self):
        """Access all carbon atoms within ``Ctab`` atom block.

        :return: List of carbon atoms.
        :rtype: :py:class:`list`.
        """
        return self['Ctab'].carbon_atoms

    @property
    def hydrogen_atoms(self):
        """Access all hydrogen atoms within ``Ctab`` atom block.

        :return: List of hydrogen atoms.
        :rtype: :py:class:`list`.
        """
        return self['Ctab'].hydrogen_atoms

    def as_sdfile(self, data=None):
        """Create ``SDfile`` from ``Molfile``.
        
        :param dict data: Data associated with ``Molfile``. 
        :return: New ``SDfile`` instance.
        :rtype: :class:`~ctfile.ctfile.SDfile`.
        """
        return SDfile.from_molfile(molfile=self, data=data)

    def add_ctab_property(self, ctab_property_name, values):
        """Add new values to existing ``CtabPropertiesBlock``.

        :param str ctab_property_name: Name of the ``Ctab`` property.
        :param values: Sequence of values. 
        :return: None.
        :rtype: :py:obj:`None`.
        """
        self['Ctab'].add_ctab_property(ctab_property_name=ctab_property_name, values=values)

    def replace_ctab_property(self, ctab_property_name, values):
        """Replace with new values ``CtabPropertiesBlock``.

        :param str ctab_property_name: Name of the ``Ctab`` property.
        :param values: Sequence of values. 
        :return: None.
        :rtype: :py:obj:`None`.
        """
        self['Ctab'].replace_ctab_property(ctab_property_name=ctab_property_name, values=values)

    def add_charge(self, atom_number, atom_symbol, charge):
        """Add charge to neutral "N", "O", or "S".

        :param str atom_number: Atom id in order of appearance in ``CTfile``.
        :param str atom_symbol: Atom symbol.
        :param str charge: Charge: "+1" or "-1".
        :return: None.
        :rtype: :py:obj:`None`.
        """
        self['Ctab'].add_charge(atom_number=atom_number, atom_symbol=atom_symbol, charge=charge)


class SDfile(CTfile):
    """SDfile - each structure-data file contains structures and data for any number
    of molecules.

    ---------------------
    | SDfile       .    |
    |            .      |
    |          .        |
    | ----------------- |
    | | ------------- | |
    | | | molfile   | | |
    | | | or RGfile | | |
    | | ------------- | |
    | | ------------- | |
    | | | data      | | |
    | | | block     | | |
    | | ------------- | |
    | ----------------- |
    ---------------------
    """
    def __init__(self):
        """SDfile initializer."""
        super(SDfile, self).__init__()

    @classmethod
    def from_molfile(cls, molfile, data=None):
        """Construct new ``SDfile`` object from ``Molfile`` object.
        
        :param molfile: ``Molfile`` object.
        :type molfile: :class:`~ctfile.ctfile.Molfile`.
        :return: ``SDfile`` object.
        :rtype: :class:`~ctfile.ctfile.SDfile`.
        """
        if not data:
            data = OrderedDict()

        if not isinstance(molfile, Molfile):
            raise ValueError('Not a Molfile type: "{}"'.format(type(molfile)))

        if not isinstance(data, dict):
            raise ValueError('Not a dict type: "{}"'.format(type(data)))

        sdfile = cls()
        sdfile['1'] = OrderedDict()
        sdfile['1']['molfile'] = molfile
        sdfile['1']['data'] = data
        return sdfile

    def add_data(self, id, key, value):
        """Add new data item.
        
        :param str id: Entry id within ``SDfile``.
        :param str key: Data item key.
        :param str value: Data item value.
        :return: None.
        :rtype: :py:obj:`None`.
        """
        self[str(id)]['data'].setdefault(key, [])
        self[str(id)]['data'][key].append(value)

    def add_molfile(self, molfile, data):
        """Add ``Molfile`` and data to ``SDfile`` object.
        
        :param molfile: ``Molfile`` instance.
        :type molfile: :class:`~ctfile.ctfile.Molfile`.
        :param dict data: Data associated with ``Molfile`` instance. 
        :return: None.
        :rtype: :py:obj:`None`.
        """
        if not isinstance(molfile, Molfile):
            raise ValueError('Not a Molfile type: "{}"'.format(type(molfile)))

        if not isinstance(data, dict):
            raise ValueError('Not a dict type: "{}"'.format(type(data)))

        entry_ids = sorted(self.keys(), key=lambda x: int(x))
        if entry_ids:
            last_entry_id = str(entry_ids[-1])
        else:
            last_entry_id = '0'

        new_entry_id = str(int(last_entry_id) + 1)
        self[new_entry_id] = OrderedDict()
        self[new_entry_id]['molfile'] = molfile
        self[new_entry_id]['data'] = data

    def add_sdfile(self, sdfile):
        """Add new ``SDfile`` to current ``SDfile``.
        
        :param sdfile: ``SDfile`` instance. 
        :return: None.
        :rtype: :py:obj:`None`.
        """
        if not isinstance(sdfile, SDfile):
            raise ValueError('Not a SDfile type: "{}"'.format(type(sdfile)))

        for entry_id in sdfile:
            self.add_molfile(molfile=sdfile[entry_id]['molfile'],
                             data=sdfile[entry_id]['data'])

    def _build(self, lexer):
        """Build :class:`~ctfile.ctfile.SDfile` instance.

        :return: :class:`~ctfile.ctfile.SDfile` instance.
        :rtype: :class:`~ctfile.ctfile.SDfile`.
        """
        current_entry_id = 0

        while True:
            token = next(lexer)
            key = token.__class__.__name__

            if key == 'MolfileStart':
                current_entry_id += 1
                molfile = Molfile()
                molfile._build(lexer)
                self[str(current_entry_id)] = OrderedDict(molfile=molfile, data=OrderedDict())

            elif key == 'DataBlockStart':
                data_block = self._build_data_block(lexer)
                self[str(current_entry_id)]['data'].update(data_block)

            elif key == 'EndOfFile':
                break

            else:
                raise KeyError('SDfile does not supposed to have any other information: "{}".'.format(key))

        return self

    def _build_data_block(self, lexer):
        """Build the data block of :class:`~ctfile.ctfile.SDfile` instance. 
        
        :return: Data block.
        :rtype: :py:class:`collections.OrderedDict`.
        """
        data_block = OrderedDict()
        header = ''

        while True:
            token = next(lexer)
            key = token.__class__.__name__

            if key == 'DataHeader':
                header = token.header[1:-1]
                data_block.setdefault(header, [])

            elif key == 'DataItem':
                data_block[header].append(token.data_item)

            elif key == 'DataBlockEnd':
                break

            else:
                raise KeyError('SDfile data block does not supposed to have any other information: "{}".'.format(key))

        return data_block

    def _to_ctfile(self):
        """Convert :class:`~ctfile.ctfile.CTfile` into `CTfile` formatted string.

        :return: ``CTfile`` formatted string.
        :rtype: :py:class:`str`.
        """
        output = io.StringIO()

        for entry in self.values():
            output.write(entry['molfile']._to_ctfile())

            for header, values in entry['data'].items():
                output.write('> <{}>\n'.format(header))
                output.write('\n'.join(values))
                output.write('\n')
            output.write('\n$$$$\n')

        return output.getvalue()

    @property
    def molfiles(self):
        """Create list of ``Molfile`` instances.
        
        :return: List of ``Molfile`` instances.
        :rtype: :py:class:`list`. 
        """
        return [entry['molfile'] for entry in self.values()]


class Atom(object):
    """Atom within ``Ctab`` block."""

    atom_block_format = 'xxxxxxxxxxyyyyyyyyyyzzzzzzzzzzaaaaddcccssshhhbbbvvvHHHrrriiimmmnnneee'

    def __init__(self, atom_number, atom_symbol, x, y, z, mass_difference, charge, atom_stereo_parity,
                 hydrogen_count, stereo_care_box, valence, h0designator, not_used1, not_used2,
                 atom_atom_mapping_number, inversion_retention_flag, exact_change_flag):
        """Atom initializer.
        
        :param str atom_number: Atom id in order of appearance in ``CTfile``.
        :param str atom_symbol: Atom symbol.
        :param str x: Atom x coordinate.
        :param str y: Atom y coordinate.
        :param str z: Atom z coordinate.
        :param str mass_difference: Atom mass difference.
        :param str charge: Atom charge.
        :param str atom_stereo_parity: Atom stereo parity.
        :param str hydrogen_count: Hydrogen count.
        :param str stereo_care_box: Atom stereo care.
        :param str valence: Atom valence.
        :param str h0designator: H0 designator.
        :param str not_used1: Unused field.
        :param str not_used2: Unused field.
        :param str atom_atom_mapping_number: Atom-atom mapping.
        :param str inversion_retention_flag: Inversion/retention flag.
        :param str exact_change_flag: Exact change flag.
        """
        self.atom_number = atom_number
        self.neighbors = []
        self._ctab_data = OrderedDict()

        self._ctab_data['x'] = x
        self._ctab_data['y'] = y
        self._ctab_data['z'] = z
        self._ctab_data['atom_symbol'] = atom_symbol
        self._ctab_data['mass_difference'] = mass_difference
        self._ctab_data['charge'] = charge
        self._ctab_data['atom_stereo_parity'] = atom_stereo_parity
        self._ctab_data['hydrogen_count'] = hydrogen_count
        self._ctab_data['stereo_care_box'] = stereo_care_box
        self._ctab_data['valence'] = valence
        self._ctab_data['h0designator'] = h0designator
        self._ctab_data['not_used1'] = not_used1
        self._ctab_data['not_used2'] = not_used2
        self._ctab_data['atom_atom_mapping_number'] = atom_atom_mapping_number
        self._ctab_data['inversion_retention_flag'] = inversion_retention_flag
        self._ctab_data['exact_change_flag'] = exact_change_flag

    def neighbor_atoms(self, atom_symbol=None):
        """Access neighbor atoms.
        
        :param str atom_symbol: Atom symbol.
        :return: List of neighbor atoms.
        :rtype: :py:class:`list`.
        """
        if not atom_symbol:
            return self.neighbors
        else:
            return [atom for atom in self.neighbors if atom['atom_symbol'] == atom_symbol]

    @property
    def neighbor_carbon_atoms(self, atom_symbol='C'):
        """Access neighbor carbon atoms.
        
        :param str atom_symbol: Atom symbol.
        :return: List of neighbor carbon atoms.
        :rtype: :py:class:`list`.
        """
        return self.neighbor_atoms(atom_symbol=atom_symbol)

    @property
    def neighbor_hydrogen_atoms(self, atom_symbol='H'):
        """Access neighbor hydrogen atoms.
        
        :param str atom_symbol: Atom symbol.
        :return: List of neighbor hydrogen atoms.
        :rtype: :py:class:`list`.
        """
        return self.neighbor_atoms(atom_symbol=atom_symbol)

    def __getitem__(self, item):
        """Provide dict-like item access to atom ``Ctab`` data."""
        return self._ctab_data[item]

    def __setitem__(self, key, value):
        """Provide dict-like item setting to atom ``Ctab`` data."""
        self._ctab_data[key] = value

    def __getattr__(self, item):
        """Provide dot item access to atom ``Ctab`` data."""
        return self._ctab_data[item]

    def __str__(self):
        """String representation of atom ``Ctab`` data."""
        return str(self._ctab_data)

    def __repr__(self):
        """Representation of atom ``Ctab`` data."""
        return str(self._ctab_data)


class Bond(object):
    """Bond that connects two atoms within ``Ctab`` block."""

    bond_block_format = '111222tttsssxxxrrrccc'

    def __init__(self, first_atom, second_atom, bond_type, bond_stereo,
                 not_used1, bond_topology, reacting_center_status):
        """Bond initializer.
        
        :param first_atom: Atom object.
        :type first_atom: :class:`~ctfile.ctfile.Atom`
        :param second_atom: Atom object.
        :type second_atom: :class:`~ctfile.ctfile.Atom`
        :param str bond_type: Bond type.
        :param str bond_stereo: Bond stereo.
        :param str not_used1: Unused field.
        :param str bond_topology: Bond topology.
        :param str reacting_center_status: Reacting center status.
        """
        self.first_atom = first_atom
        self.second_atom = second_atom
        self._ctab_data = OrderedDict()

        self._ctab_data['first_atom_number'] = first_atom.atom_number
        self._ctab_data['second_atom_number'] = second_atom.atom_number
        self._ctab_data['bond_type'] = bond_type
        self._ctab_data['bond_stereo'] = bond_stereo
        self._ctab_data['not_used1'] = not_used1
        self._ctab_data['bond_topology'] = bond_topology
        self._ctab_data['reacting_center_status'] = reacting_center_status

    def update_atom_numbers(self):
        """Update links "first_atom_number" -> "second_atom_number"
        
        :return: None.
        :rtype: :py:obj:`None`.
        """
        self._ctab_data['first_atom_number'] = self.first_atom.atom_number
        self._ctab_data['second_atom_number'] = self.second_atom.atom_number

    def __getitem__(self, item):
        """Provide dict-like item access to bond ``Ctab`` data."""
        return self._ctab_data[item]

    def __setitem__(self, key, value):
        """Provide dict-like item setting to bond ``Ctab`` data."""
        self._ctab_data[key] = value

    def __getattr__(self, item):
        """Provide dot item access to bond ``Ctab`` data."""
        return self._ctab_data[item]

    def __str__(self):
        """String representation of bond ``Ctab`` data."""
        return str(self._ctab_data)

    def __repr__(self):
        """Representation of bond ``Ctab`` data."""
        return str(self._ctab_data)


class CtabAtomBondEncoder(json.JSONEncoder):
    """Custom serializer for Atom and Bond objects."""

    def default(self, o):
        """Default encoder.

        :param o: Atom or Bond instance.
        :type o: :class:`~ctfile.ctfile.Atom` or :class:`~ctfile.ctfile.Bond`.
        :return: Dictionary that contains information required for atom and bond block of ``Ctab``.
        :rtype: :py:class:`collections.OrderedDict`
        """
        if isinstance(o, Atom) or isinstance(o, Bond):
            return o._ctab_data
        else:
            return o.__dict__
