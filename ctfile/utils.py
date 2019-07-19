#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

"""
ctfile.utils
~~~~~~~~~~~~

This module provides reusable utility functions and objects.
"""

from collections import OrderedDict
from collections import Counter


class OrderedCounter(Counter, OrderedDict):
    """Counter that remembers the order of elements that are seen first."""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, OrderedDict(self)
