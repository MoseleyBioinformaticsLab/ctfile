#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from filehandles import filehandles
from .ctfile import CTfile


load = CTfile.load
loadstr = CTfile.loadstr


def read_file(path, verbose=False):
    """Read single ``CTfile`` formatted file.
    
    :param str path: Path to ``CTfile``. 
    :return: Subclass of :class:`~ctfile.ctfile.CTfile`.
    :rtype: :class:`~ctfile.ctfile.CTfile`
    """
    for filehandle in filehandles(path, verbose=verbose):
        try:
            return CTfile.load(filehandle)
        except:
            raise ValueError('Cannot load file from path: "{}"'.format(path))


def read_files(path, verbose=True):
    """Read multiple ``CTfile`` formatted files.
    
    :param str path: Path to ``CTfile``. 
    :return: Subclass of :class:`~ctfile.ctfile.CTfile`.
    :rtype: :class:`~ctfile.ctfile.CTfile`
    """
    for filehandle in filehandles(path, verbose=verbose):
        yield CTfile.load(filehandle)