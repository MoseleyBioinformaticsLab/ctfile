#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""The ``ctfile`` package is a Python library that facilitates 
reading and writing of ``CTfile`` formats. CTfile stands for 
Chemical Table file which is a family of text-based chemical 
file formats and is used to describe chemical molecules and 
reactions.

This package includes the following modules:

``api``
    This modules provides core API functions and objects.
   
``conf``
    This module reads in configuration files.

``ctfile``
    This module implements core objects to represent ``CTfile`` file format.
   
``tokenizer``
    This module implements tokenizer that processes text file and yields 
    tokens necessary to create ``CTfile`` objects.

``utils``
    This module provides reusable utility functions and objects.
"""

import logging

from .api import read_file
from .api import read_files
from .ctfile import CTfile
from .ctfile import Molfile
from .ctfile import SDfile
from .ctfile import Atom
from .ctfile import Bond


load = CTfile.load
loadstr = CTfile.loadstr


__version__ = '0.1.7'


try:  # Python 2/3 compatibility code
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass


# Setting default logging handler
logging.getLogger(__name__).addHandler(NullHandler())
