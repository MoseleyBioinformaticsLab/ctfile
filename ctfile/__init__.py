#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging

from .api import load
from .api import loadstr
from .api import read_file
from .api import read_files

__version__ = '0.1.4'


try:  # Python 2/3 compatibility code
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass


# Setting default logging handler
logging.getLogger(__name__).addHandler(NullHandler())
