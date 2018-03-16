#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

import os
import json
from collections import OrderedDict
from collections import Counter


this_directory = os.path.abspath(os.path.dirname(__file__))
ctab_properties_path = '{}/conf/Ctab_properties.json'.format(this_directory)
with open(ctab_properties_path, 'r') as infile:
    ctab_properties_conf = json.load(infile)


class OrderedCounter(Counter, OrderedDict):
    """Counter that remembers the order of elements are first seen."""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, OrderedDict(self)
