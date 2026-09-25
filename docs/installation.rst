============
Installation
============

M2M-PostAViz is tested with Python version 3.12, 3.13 and 3.14.

At the command line::

    pip install m2m-postaviz

To install the latest development version from source::

    git clone https://gitlab.inria.fr/postaviz/m2m-postaviz.git
    cd m2m-postaviz
    pip install -r requirements.txt
    pip install .

Dependencies
============

M2M-PostAViz dependencies (installed automatically with pip):

- pandas
- padmet
- scipy
- skbio
- plotly
- scikit-bio
- shiny
- shinywidgets
- faicons
- matplotlib
- polars
- pyarrow
- seaborn
- statsmodels

To install dependencies manually::

    pip install -r requirements.txt

Troubleshooting
===============

Polars compatibility (``polars``)
---------------------------------

If installation succeeds but running ``m2m_postaviz`` fails because ``polars`` is unavailable or incompatible on your machine, install the long-term-support CPU wheel manually::

    pip install polars-lts-cpu

macOS ARM binary mismatch (``biom-format`` / ``scikit-bio``)
----------------------------------------------------------------

On Apple Silicon (arm64), you may occasionally hit an error like:

- ``ImportError: ... biom/_filter...so ... incompatible architecture (have 'x86_64', need 'arm64')``

This means a transitive binary extension (from ``biom-format``, required by ``scikit-bio``) was installed for the wrong architecture.

In most cases, reinstalling ``biom-format`` with arm64 build flags is enough::

    pip uninstall -y biom-format
    ARCHFLAGS="-arch arm64" pip install --no-binary=biom-format biom-format==2.1.17

Then verify with::

    m2m_postaviz --test


Citation
========

If you use the application for research, do not forget to cite the works associated to those dependencies.

Quickstart
==========

To test the application with example data:

.. code-block:: bash

    m2m_postaviz --test

It takes a few seconds to launch because data needs to be uncompressed and processed in a temporary directory. Shiny will launch automatically afterward.

A dataset with test data is available in this repository in ``postaviz_test_data`` and can be used to test the main functionalities of the tool.

Once on the homepage you're free to explore the test data.

.. note::
    The ``--test`` option does not load MetaCyc ontology data; therefore,
    metabolite-family exploration is unavailable in this mode. The bundled
    file ``metacyc28_5.padmet`` is a derived fixture used solely by the
    automated test suite. It was generated from MetaCyc release 28.5 using
    the PADMet library and subsequently filtered by the authors to remove
    reaction and pathway records, retaining only the compound annotations,
    cross-references, and classification relationships required by the
    tests. It is not intended as a general-purpose MetaCyc distribution or
    as input for user analyses. Users requiring MetaCyc data should obtain
    them directly from MetaCyc and comply with the applicable licence.

.. warning::
   We assume that you arrive at this step having installed the tool first (see above), for instance in a Python virtual environment, or conda (mamba) environment.
