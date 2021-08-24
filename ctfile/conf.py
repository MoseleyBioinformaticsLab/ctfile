#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

"""
ctfile.conf
~~~~~~~~~~~

This module reads in configuration file.
"""

import os
import json


this_directory = os.path.abspath(os.path.dirname(__file__))
ctab_properties_path = '{}/config_files/Ctab_properties.json'.format(this_directory)
charge_index_path = '{}/config_files/Charge_index.json'.format(this_directory)
index_charge_path = '{}/config_files/Index_charge.json'.format(this_directory)

with open(ctab_properties_path, 'r') as infile:
    ctab_properties_conf = json.load(infile)


with open(charge_index_path, 'r') as infile:
    charge_index = json.load(infile)


with open(index_charge_path, 'r') as infile:
    index_charge = json.load(infile)
