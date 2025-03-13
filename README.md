# M2M-PostAViz

## Overview

M2M-PostAViz (_M2M Post-Analysis and Visualization_) is an interactive exploration platform for the exploration of metabolic potential predictions performed by [Metage2Metabo (M2M)](https://github.com/AuReMe/metage2metabo/tree/main). M2M predicts reachable metabolites by communities of microorganisms by computing a Boolean abstract abstraction (Network expansion) of metabolic activity. More precisely, it takes as input a collection of genome-scale metabolic networks (GSMNs) and a description of the nutritional environment, then first compute the molecules reachable (predicted to be producible in the environment) by each GSMN individually, and by the community of GSMNs. In the latter case, it takes into account mutualistic interactions (cross-feeding) that may occur within the community, and therefore increase its overall metabolic potential. As such, several outputs can be distinguished: what is produced by each member alone or _individual scope_, what becomes producible only through metabolic interactions _added-value_ and the union of both that is the community metabolic potential, or _community scope_. For each compound in one of the _scopes_, the information related to the producer(s) -- i.e. which GSMNs -- of the compound is retrieved by M2M and provided in dedicated tables. The overview of M2M pipeline is illustrated below:

<img src="./docs/pictures/m2m_overview.png" alt="General overview of Metage2Metabo's pipeline" width="50%"/>

M2M-PostAViz integrates and analyses all these outputs, especially in a context when multiple runs of M2M are performed, each one aiming at studying the metabolic potential of a microbial community. A typical use-case would be to run M2M for a cohort of microbiome samples, each described by a collection of GSMNs, for instance. In a cohort, the data comes with metadata that is used by M2M-PostAViz to analyse M2M's results and explore whether the predicted metabolic potential is statistically associated with metadata variables. 

<img src="./docs/pictures/postaviz_overview.png" alt="General overview of M2M-PostAViz" width="50%"/>


### License 
GNU Lesser General Public License v3 (LGPLv3)

## Installation

```pip install .```

    pip install m2m-postaviz

You can also install the in-development version with::

    pip install git+ssh://git@gitlab.inria.fr/postaviz/m2m-postaviz.git@main

### Dependencies


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

` pip install -r requirement.txt `

## Quick start

A dataset with test data is available in this repository and can be used to test the main functionalities of the tool.

##  Documentation



### Important

In metadata file the first column must be the sample identification. Preferably named "smplID".

In taxonomy file the first column must be the metagenomes (mgs). Preferably named "mgs".

### Usage

m2m_postaviz can be run in two ways :

```
    m2m_postaviz -d Metage2metabo/samples/scopes/directory/path
                -m metadata/file/path
                -a abundance/file/path
                -t taxonomy/file/path
                -o save/path
                --no-metacyc (Optionnal)
```

This way is required as least one time to produce all dataframe and save them in -o save/path.

Once the dataframes are produced. Shiny will automatically run from the save/path given in -o option.
You can interrupt the process if you want and run postaviz with -l load option.

```
    m2m_postaviz -l save/directory/path
````

Which will directly launch shiny and skip dataprocessing.

### Metadata tab

Tabulation to observe the metadata given in CLI.
For future update we will use this tab to allow users to change metadata type directly inside the application.

### PCOA tabulation

- Base PCOA shows the Principal coordinates analysis done with all samples and saves as pcoa_dataframe_postaviz.tsv into the save directory.

You can use the input color to observe the PCOA in all its shapes. But the PCOA is not recalculated, this is just for global observation.

- Customisable PCOA allow to calculate directly the PCOA from the samples filtered by the metadata input. The color can be used to regroup samples by their metadata values.

### Bins exploration tabulation


Tab dedicated to the observation of the bins contained into each sample's cscope of Metage2metabo.
Some pattern in compounds production can be found by the taxonomic belonging of the bins.
If the taxonomy (-t option) is not provided, this tabulation will be disabled.

Input:

- Allow to choose between the taxonomic ranks, the individual metagenomes "mgs" or all bins with "all"

- The second input automatically update from the selection above. It allow the selection of the specific group of bins in taxonomic rank selected.

- The third input allow a filtering to the samples level, all samples (and associated bins !) will be removed from the plots if excluded by this input.

- Updated from the third input like the second input, allow a more precise selection ????

- Color grouping option for all plots.

- Use abundance, this options will use the "normalised abundance dataframe" instead of the "main dataframe". Instead of using 0,1 value for production, the abundance dataframe is multi with the abundace of each bins in their respective sample.

Plots :

- Plot 1 the sum of unique metabolites produced by the selected bins in each samples.

- Plot 2 is a boxplot of the unique metabolites production of each selected bins in their samples.

- Plot 3 show the abundance for each selected bins in their respective sample.

>WARNING

>The "all" option on all sample (No metadata filter applied) can be long to produce the plots. Also heavy plots will impact the performance of the application. 

> [!NOTE]

> A small text output under the Processing button show how many bins are selected to avoid large calculation. Also if only one bin (mgs) is selected, it will display how many samples have this specfic bin.

#### Method

The bins exploration present some challenges, it need to retain the production by the bins inside of each samples.

The individual production of each bins depend of the bins present (community interaction highlighted by **Meta2metabo**) which also differ from the treatment/condition of the sample.

All of this information can scale poorly if a lot of sample are present in input which is the point of the application: being able observe from lots of angles large chunks of data.

In order to keep all this valuable data we choosed to use the hard drive and load/manipulate large dataframe by subset using Parquet format. That way RAM memory is not overload and we can pick what we need from the hard drive, using query similar to a SQL database.


### Compounds exploration tab

Input :

IF METACYC ENABLED

- Compounds input divided into three sub input:
    - List of metacyc category ordered from the top to the bottom of the tree.
    - Any category selected above will update this input to a list of all sub-category.
    - Automatically filled with the compounds corresponding to the category / sub category selected in the input above.
    Allow the selection of compounds directly if none of the first input are used or if **--no-metacyc** option is used in CLI.

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

#### Method

To produce the plots, m2m-postaviz use several dataframes made in processing part.

The producers_cscope_dataframe and producers_iscope_dataframe are the two main dataframe.
Both are produced by the sum of each rows (bins) in each sample's cscope/iscope and the concatenation of all Pandas Series produced.

That way we have the number of metagenomes producing the compound (columns) for each rows (samples).

The difference between these dataframe (production in community and individual) gives us the added_value dataframe




## Development

To run all the tests run::

    tox run -e clean,(pyXXX),report

Note, to combine the coverage data from all the tox environments run:

* Windows

```
    set PYTEST_ADDOPTS=--cov-append
    tox
```

* Other

```
    PYTEST_ADDOPTS=--cov-append tox
```

## Authors

Léonard Brindel and [Clémence Frioux](https://cfrioux.github.io) -- [Inria Pleiade team](https://team.inria.fr/pleiade/) 

### Acknowledgements

- David James Sherman
- Jean-Marc Frigerio
- Pleiade team members