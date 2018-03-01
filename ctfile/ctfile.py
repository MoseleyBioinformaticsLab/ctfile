#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

from __future__ import print_function, division, unicode_literals
import json
import os
from collections import OrderedDict

from .utils import OrderedCounter


class CTfile(OrderedDict):
    """"""

    def __init__(self, lexer):
        """CTfile initializer.
        
        :param lexer: 
        """
        super(CTfile, self).__init__()
        self.lexer = lexer

        this_directory = os.path.dirname(__file__)
        template_path = "{}/conf/{}_template.json".format(this_directory, self.__class__.__name__)

        with open(template_path, "r") as infile:
            self.update(json.load(infile, object_pairs_hook=OrderedDict))

        self._build()

    def read(self, filehandle):
        """

        :param filehandle:
        :return:
        """
        input_str = filehandle.read()

    def write(self, filehandle, file_format):
        """

        :param filehandle:
        :param file_format:
        :return:
        """
        raise NotImplementedError("Abstract class")

    def _build(self):
        """

        :return:
        """
        raise NotImplementedError("Abstract class")

    def _to_json(self):
        """Convert :class:`~ctfile.ctfile.CTfile` into JSON string.
        
        :return: JSON string.
        :rtype: :py:class:`str`
        """
        return json.dumps(self, sort_keys=False, indent=4)

    def _to_ctfile(self):
        """Convert :class:`~ctfile.ctfile.CTfile` into `CTfile` string.
        
        :return: CTfile string.
        :rtype: :py:class:`str`
        """
        pass

    @staticmethod
    def _is_molfile_str(string):
        """Test if input string is in ``Molfile`` format.

        :param string: Input string.
        :type string: :py:class:`str` or :py:class:`bytes`
        :return: Input string if in ``Molfile`` format or False otherwise.
        :rtype: :py:class:`str` or :py:obj:`False`
        """
        pass

    @staticmethod
    def _is_sdfile_str(string):
        """

        :param str string:
        :return:
        """
        pass


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
    counts_line_format = "aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv"
    atom_block_format = "xxxxxxxxxxyyyyyyyyyyzzzzzzzzzzaaaaddcccssshhhbbbvvvHHHrrriiimmmnnneee"
    bond_block_format = "111222tttsssxxxrrrccc"

    def _build(self):
        for token in self.lexer:
            key = token.__class__.__name__

            if key == "CtabCountsLine":
                self[key].update(token._asdict())

            elif key in ("CtabAtomBlock", "CtabBondBlock"):
                self[key].append(token._asdict())

            elif key == "CtabPropertiesBlock":
                self[key].append(token._asdict())

    def write(self, filehandle, file_format):
        """
        
        :param filehandle: 
        :param file_format: 
        :return: 
        """
        for key in self:

            if key == "CtabCountsLine":
                counter = OrderedCounter(self.counts_line_format)
                counts_line = "".join([str(value).rjust(spacing) for value, spacing in zip(self[key].values(), counter.values())])
                filehandle.write(bytes(counts_line, "utf-8"))
                filehandle.write(bytes("\n", "utf-8"))

            elif key == "CtabAtomBlock":
                counter = OrderedCounter(self.atom_block_format)

                for i in self[key]:
                    atom_line = "".join([str(value).rjust(spacing) for value, spacing in zip(i.values(), counter.values())])
                    filehandle.write(bytes(atom_line, "utf-8"))
                    filehandle.write(bytes("\n", "utf-8"))

            elif key == "CtabBondBlock":
                counter = OrderedCounter(self.bond_block_format)

                for i in self[key]:
                    bond_line = "".join([str(value).rjust(spacing) for value, spacing in zip(i.values(), counter.values())])
                    filehandle.write(bytes(bond_line, "utf-8"))
                    filehandle.write(bytes("\n", "utf-8"))

            elif key == "CtabPropertiesBlock":
                for i in self[key]:
                    filehandle.write(bytes(i["property_line"], "utf-8"))
                    filehandle.write(bytes("\n", "utf-8"))

            else:
                raise KeyError("Ctab object does not supposed to have any other keys.")


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
    def _build(self):
        """
        
        :return: 
        """
        for token in self.lexer:
            key = token.__class__.__name__

            if key == "HeaderBlock":
                self[key].update(token._asdict())

            elif key == "CtabBlock":
                ctab = Ctab(lexer=self.lexer)
                self["Ctab"] = ctab

            else:
                pass

    def write(self, filehandle, file_format):
        """
        
        :param filehandle: 
        :param file_format: 
        :return: 
        """
        for key in self:

            if key == "HeaderBlock":
                for line in self[key].values():
                    filehandle.write(bytes(line, "utf-8"))
                    filehandle.write(bytes("\n", "utf-8"))

            elif key == "Ctab":
                self[key].write(fileahandle, file_format)

            else:
                raise KeyError("Molfile object does not supposed to have any other keys.")


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
    pass