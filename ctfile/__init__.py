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

from .api import load
from .api import loadstr
from .api import read_file
from .api import read_files
from .api import Ctab
from .api import Molfile
from .api import SDfile
from .api import Atom
from .api import Bond


__version__ = '0.1.6'


try:  # Python 2/3 compatibility code
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass


# Setting default logging handler
logging.getLogger(__name__).addHandler(NullHandler())
