#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
from setuptools import setup, find_packages


if sys.argv[-1] == 'publish':
    os.system('python3 setup.py sdist')
    os.system('twine upload dist/*')
    sys.exit()


def readme():
    with open('README.rst') as readme_file:
        return readme_file.read()


def find_version():
    with open('ctfile/__init__.py', 'r') as fd:
        version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                            fd.read(), re.MULTILINE).group(1)
    if not version:
        raise RuntimeError('Cannot find version information')
    return version


REQUIRES = [
    'docopt >= 0.6.2',
    'filehandles >= 0.3.1',
    'more_itertools >= 4.1.0'
]


setup(
    name='ctfile',
    version=find_version(),
    author='Andrey Smelter',
    author_email='andrey.smelter@gmail.com',
    description='Parser for chemical table file (CTfile) a family of text-based chemical file formats '
                '(e.g. Molfile, SDfile) that describe molecules and chemical reactions.',
    keywords='CTfile SDfile Molfile Ctab',
    license='BSD',
    url='https://github.com/MoseleyBioinformaticsLab/ctfile',
    packages=find_packages(),
    package_data={'ctfile': ['config_files/*.json']},
    platforms='any',
    long_description=readme(),
    install_requires=REQUIRES,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)