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
- seaborn

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
                --no-metacyc (Optionnal)

This way is required as least one time to produce all dataframe and save them in -o save/path.

Once the dataframes are produced. Shiny will automatically run from the save/path given in -o option.
You can interrupt the process if you want and run postaviz with -l load option.

::

    m2m_postaviz -l save/directory/path

Which will directly launch shiny and skip dataprocessing.

Metadata tabulation
===========

Tabulation to observe the metadata given in CLI.
For future update we will use this tab to allow users to change metadata type directly inside the application.

PCOA tabulation
===========

Tab with 2 graph.
- The first one shows the Principal coordinates analysis done with all samples and saves as pcoa_dataframe_postaviz.tsv into the save directory.
You can use the input color to observe the PCOA in all its shapes. But the PCOA is not recalculated, this is just for observation.

- The second graph allow the calculate directly the PCOA from the samples filtered by the metadata input. The color can be used to regroup samples
by their metadata values.

Bins exploration tabulation
===========

Tab dedicated to the observation of the bins contained into each sample's cscope of Metage2metabo.
Some pattern in compounds production can be found by the taxonomic belonging of the bins.
If the taxonomy (-t option) is not provided, this tabulation will be disabled.


Input:

- Allow to choose between the taxonomic ranks, the individual metagenomes "mgs" or all bins with "all"

- The second input automatically update from the input above. It allow the selection of the specific group of bins in CATEGORY ???

- The third input allow a filtering to the samples level, all samples (and associated bins !) will be removed from the plots if excluded by this input.

- Updated from the third input liek the second input, allow a more precise selection ????

- Color grouping option for all plots.

- Use abundance, this options will use the "normalised abundance dataframe" instead of the "main dataframe". Instead of using 0,1 value for production, the abundance dataframe is multi with the abundace of each bins in their respective sample.

Plots :

- Plot 1 the sum of unique metabolites produced by the selected bins in each samples.

- Plot 2 is a boxplot of the unique metabolites production of each selected bins in their samples.

- Plot 3 show the abundance for each selected bins in their respective sample.

.. warning::
    The "all" option on all sample (No metadata filter applied) can be long to produce the plots. Also heavy plots will impact the performance of the application. 

.. note::
    A small text output under the Processing button show how many bins are selected to avoid large calculation. Also if only 
    one bin (mgs) is selected it will display how many samples have this specfic bin.


Compounds exploration tabulation
===========

Input :

IF METACYC ENABLED

- Compounds input divided into three sub input:
    - List of metacyc category ordered from the top to the bottom of the tree.
    - Any category selected above will update this input to a list of all sub-category.
    - Automatically filled with the compounds corresponding to the category / sub category selected in the input above.
    Allow the selection of compounds directly if none of the first input are used or if --no-metacyc option is used in CLI.

The plots generated will only take the compounds selected as input.

- Metadata filter and color

    - Metadata filter
    - Plot color and regroup

- Sample filter
    - All (no filter), Inlucde or Exclude from the plots the selection in filter.
    - Metadata column selection. REDUNDANCY WITH METADATA FILTER ?????
    - unique choice of the metadata column previously selected.
    - Automatically filled. Select sample corresponding to previous choices. NO CUSTOM SELECTION BY NAME POSSIBLE !! NONE option in metadata should be enable!!


- Enable row/columns clustering (Only for heatmap) will change column and or row order. Optionnal and independant from each other.

- Generate statistical dataframe / Should be enable by default / performance 

Plot

- Heatmap

    Heatmap displaying the number of bins producing the compound in the sample.
    Cscope / Iscope / Added Value

- Percentage of samples producing selected compounds

    Divided as Cscope / Iscope

- Boxplot of the production of compounds (Y axis) by sample regrouped by metadata (X axis)

- Stats tests dataframe.

    Wilcoxon / Mann-Whitney tests for factor/categoric data.
    Corrélation test for integers data.

    Tested pair are determined by the metadata input.

    Sample filtering is not applied here.