============
Installation
============

M2M-PostAViz is tested with Python version 3.12 and 3.13.

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
- pyarrow
- seaborn

To install dependencies manually::

    pip install -r requirements.txt

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
   Metacyc database information related to the ontology of metabolites and pathways is not included in test option.

.. warning::
   We assume that you arrive at this step having installed the tool first (see above), for instance in a Python virtual environment, or conda (mamba) environment.
