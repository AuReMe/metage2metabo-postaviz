Development
===========

To run all the tests run::

    tox run -e clean,(pyXXX),report

Note, to combine the coverage data from all the tox environments run:

* Windows

.. code-block:: bash

    set PYTEST_ADDOPTS=--cov-append
    tox

* Other

.. code-block:: bash

    PYTEST_ADDOPTS=--cov-append tox
