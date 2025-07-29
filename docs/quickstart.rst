Quickstart
==========

To test the application with example data:

.. code-block:: bash

    m2m_postaviz --test

It takes a few seconds to launch because data needs to be uncompressed and processed in a temporary directory. Shiny will launch automatically afterwards.

A dataset with test data is available in this repository in ``postaviz_test_data`` and can be used to test the main functionalities of the tool.

Once on the homepage you're free to explore the test data.

.. note::
   Metacyc database information related to the ontology of metabolites and pathways is not included in test option.

.. warning::
   We assume that you arrive at this step having installed the tool first (see :doc:`installation`), for instance in a Python virtual environment, or conda (mamba) environment.

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
