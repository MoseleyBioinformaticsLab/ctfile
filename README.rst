ctfile
======

.. image:: https://img.shields.io/pypi/l/ctfile.svg
   :target: https://choosealicense.com/licenses/bsd-3-clause-clear/
   :alt: License information

.. image:: https://img.shields.io/pypi/v/ctfile.svg
   :target: https://pypi.org/project/ctfile
   :alt: Current library version

.. image:: https://img.shields.io/pypi/pyversions/ctfile.svg
   :target: https://pypi.org/project/ctfile
   :alt: Supported Python versions


The ``ctfile`` package is a Python library that facilitates reading and writing
of ``CTfile`` formats.

``CTfile`` stands for `Chemical Table file`_ which is a family of text-based chemical
file formats and is used to describe chemical molecules and reactions.

The ``ctfile`` package provides ability to parse various chemical table file formats
(``CTfile`` formats), currently ``molfile`` and ``SDfile``.


Links
~~~~~

   * ctfile @ GitHub_
   * ctfile @ PyPI_


Installation
~~~~~~~~~~~~

The ``ctfile`` package runs under Python 2.7 and Python 3.4+. Use pip_ to install.

Install on Linux, Mac OS X
--------------------------

.. code:: bash

   python3 -m pip install ctfile

Install on Windows
------------------

.. code:: bash

   py -3 -m pip install ctfile


Quickstart
~~~~~~~~~~

.. code:: python

   >>> import ctfile
   >>>
   >>> with open('path/to/ctfile', 'r') as infile:
   >>>     ctf = ctfile.load(infile)
   >>>

TODO
~~~~

* Add support for version ``V3000`` reading  of ``CTfile`` formatted files.
* Add command-line interface.
* Add ability to convert between ``V2000`` and ``V3000`` syntax.
* Add proper package documentation.


License
~~~~~~~

This package is distributed under the BSD_ `license`.



.. _Chemical Table file: https://en.wikipedia.org/wiki/Chemical_table_file
.. _pip: https://pip.pypa.io
.. _PyPI: https://pypi.org/project/ctfile
.. _GitHub: https://github.com/MoseleyBioinformaticsLab/ctfile
.. _BSD: https://choosealicense.com/licenses/bsd-3-clause-clear/