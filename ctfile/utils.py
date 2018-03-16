#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

from collections import OrderedDict
from collections import Counter


class OrderedCounter(Counter, OrderedDict):
    """Counter that remembers the order of elements are first seen."""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, OrderedDict(self)
