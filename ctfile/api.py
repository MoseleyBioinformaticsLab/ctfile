#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ctfile.api
~~~~~~~~~~

This module provides core API functions and objects.
"""

from filehandles import filehandles
from .ctfile import CTfile


def read_file(path, verbose=False):
    """Read single ``CTfile`` formatted file.
    
    :param str path: Path to ``CTfile``.
    :param verbose: Print what files are processing.
    :type verbose: :py:obj:`True` or :py:obj:`False`.
    :return: Subclass of :class:`~ctfile.ctfile.CTfile`.
    :rtype: :class:`~ctfile.ctfile.CTfile`
    """
    for filehandle in filehandles(path, verbose=verbose):
        return CTfile.load(filehandle)


def read_files(path, verbose=True):
    """Read multiple ``CTfile`` formatted files.
    
    :param str path: Path to ``CTfile``.
    :param verbose: Print what files are processing.
    :type verbose: :py:obj:`True` or :py:obj:`False`.
    :return: Subclass of :class:`~ctfile.ctfile.CTfile`.
    :rtype: :class:`~ctfile.ctfile.CTfile`
    """
    for filehandle in filehandles(path, verbose=verbose):
        yield CTfile.load(filehandle)
