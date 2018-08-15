#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ctfile.api
~~~~~~~~~~

This module provides core API functions and objects.
"""

from filehandles import filehandles
from .ctfile import CTfile
from .ctfile import Ctab
from .ctfile import Molfile
from .ctfile import SDfile
from .ctfile import Atom
from .ctfile import Bond


load = CTfile.load
loadstr = CTfile.loadstr


def read_file(path, verbose=False):
    """Read single ``CTfile`` formatted file.
    
    :param str path: Path to ``CTfile``. 
    :return: Subclass of :class:`~ctfile.ctfile.CTfile`.
    :rtype: :class:`~ctfile.ctfile.CTfile`
    """
    for filehandle in filehandles(path, verbose=verbose):
        return CTfile.load(filehandle)


def read_files(path, verbose=True):
    """Read multiple ``CTfile`` formatted files.
    
    :param str path: Path to ``CTfile``. 
    :return: Subclass of :class:`~ctfile.ctfile.CTfile`.
    :rtype: :class:`~ctfile.ctfile.CTfile`
    """
    for filehandle in filehandles(path, verbose=verbose):
        yield CTfile.load(filehandle)