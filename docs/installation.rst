============
Installation
============

M2M-PostAViz is tested with Python version 3.10, 3.11 and 3.12.

At the command line::

    pip install m2m-postaviz

To install the latest development version from source::

    git clone https://gitlab.inria.fr/postaviz/m2m-postaviz.git
    cd m2m-postaviz
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
