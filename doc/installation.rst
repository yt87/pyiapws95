Installation
============

Prerequisites
-------------

* Required
   - Python (3.8 or later)
   - `numpy <https://numpy.org/doc/stable/index.html>`__

* Optional, but highly recommended
   - `numba <https://numba.pydata.org/numba-doc/latest/index.html>`__
   - `pint <https://pint.readthedocs.io/en/stable/>`__

* Testing
   - `pytest <https://github.com/pytest-dev/pytest>`__

Instructions
------------

From source::

  $ python setup.py install
  $ pytest tests

With conda::

  $ conda install -c yt87 pyiapws95
