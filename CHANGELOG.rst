.. :changelog:

Release History
===============

0.1.7 (2019-07-17)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Added methods to add/remove isotope property of atoms.
- Added methods to add/remove charge property of atoms.
- Add atom coloring algorithm.


0.1.6 (2018-08-20)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Added new functionality to add charge to "N", "O", and "S" atoms
  within neutral molecules.

**Bugfixes**

- Fixed bug where numeration within ``SDfile`` starts
  from "2" instead of "1".


0.1.5 (2018-08-15)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Added "Atom" object as part of ``Ctab`` block that contains information
  about neighbor atoms and additional functionality.
- Added "Bond" object as part of ``Ctab`` block.
- Added custom serializer for "Atom" and "Bond" objects.
- Added module-level documentation.
- Improved representation of ``Ctab`` properties block.
- Improved ``Ctab`` properties configuration file.


0.1.4 (2018-04-18)
~~~~~~~~~~~~~~~~~~

**Bugfixes**

- Fixed Python 2.7 and Python 3.4 compatibility.


0.1.2 - 0.1.3 (2018-04-17)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bugfixes**

- Fixed bug of not including configuration file into source distribution.


0.1.1 (2018-04-16)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Added ability to access list of ``Molfiles`` from ``SDfile`` instances.


0.1.0 (2018-04-04)
~~~~~~~~~~~~~~~~~~

- Initial public release.