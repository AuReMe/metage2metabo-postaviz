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

Dependencies
============


- pandas
- padmet
- scipy
- skbio
- plotly
- scikit-bio
- shiny
- shinywidgets
- pyarrow

::

    pip install -r requirement.txt


Documentation
=============


https://m2m-postaviz.readthedocs.io/


Development
===========

To run all the tests run::

    tox run -e clean,(pyXXX),report

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

Important
===========

In metadata file the first column must be the sample identification. Preferably named "smplID".

In taxonomy file the first column must be the metagenomes (mgs). Preferably named "mgs".

Utilisation
===========

m2m_postaviz can be run in two ways :

::

    m2m_postaviz -d Metage2metabo/samples/scopes/directory/path
                -m metadata/file/path
                -a abundance/file/path
                -t taxonomy/file/path
                -o save/path

This way is required as least one time to produce all dataframe and save them in -o save/path.

Once the dataframes are produced. Shiny will automatically run from the save/path given in -o option.
You can interrupt the process if you want and run postaviz with -l load option.

::

    m2m_postaviz -l save/directory/path

Which will directly launch shiny and skip dataprocessing.