========
Overview
========

Post-analysis and visualization of metage2metabo's results

* Free software: GNU Lesser General Public License v3 (LGPLv3)

Installation
============

::

    pip install m2m-postaviz

You can also install the in-development version with::

    pip install git+ssh://git@gitlab.inria.fr/postaviz/m2m-postaviz.git@main

Documentation
=============


https://m2m-postaviz.readthedocs.io/


Development
===========

To run all the tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
