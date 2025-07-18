[![PyPI version](https://img.shields.io/pypi/v/metage2metabo.svg)](https://pypi.org/project/Metage2Metabo/) [![GitHub license](https://img.shields.io/github/license/AuReMe/metage2metabo-postaviz.svg)](https://github.com/AuReMe/metage2metabo-postaviz/blob/main/LICENSE) [![Actions Status](https://github.com/AuReMe/metage2metabo/actions/workflows/pythonpackage.yml/badge.svg)](https://github.com/AuReMe/metage2metabo/actions/workflows/pythonpackage.yml) [![Documentation Status](https://readthedocs.org/projects/metage2metabo-postaviz/badge/?version=latest)](https://metage2metabo-postaviz.readthedocs.io/en/latest/?badge=latest)

# M2M-PostAViz

## Overview

M2M-PostAViz (_M2M Post-Analysis and Visualization_) is an interactive exploration platform for the exploration of metabolic potential predictions performed by [Metage2Metabo (M2M)](https://github.com/AuReMe/metage2metabo/tree/main). M2M predicts reachable metabolites by communities of microorganisms by computing a Boolean abstract abstraction (Network expansion) of metabolic activity. More precisely, it takes as input a collection of genome-scale metabolic networks (GSMNs) and a description of the nutritional environment, then first compute the molecules reachable (predicted to be producible in the environment) by each GSMN individually, and by the community of GSMNs. In the latter case, it takes into account mutualistic interactions (cross-feeding) that may occur within the community, and therefore increase its overall metabolic potential. As such, several outputs can be distinguished: what is produced by each member alone or _individual scope_, what becomes producible only through metabolic interactions _added-value_ and the union of both that is the community metabolic potential, or _community scope_. For each compound in one of the _scopes_, the information related to the producer(s) -- i.e. which GSMNs -- of the compound is retrieved by M2M and provided in dedicated tables. The overview of M2M pipeline is illustrated below:

<img src="./docs/pictures/m2m_overview.png" alt="General overview of Metage2Metabo's pipeline" width="70%"/>

M2M-PostAViz integrates and analyses all these outputs, especially in a context when multiple runs of M2M are performed, each one aiming at studying the metabolic potential of a microbial community. A typical use-case would be to run M2M for a cohort of microbiome samples, each described by a collection of GSMNs, for instance. In a cohort, the data comes with metadata that is used by M2M-PostAViz to analyse M2M's results and explore whether the predicted metabolic potential is statistically associated with metadata variables. 

<img src="./docs/pictures/postaviz_overview.png" alt="General overview of M2M-PostAViz" width="70%"/>


### License 
GNU Lesser General Public License v3 (LGPLv3)

## Installation

M2M-PostAViz is tested with Python version 3.10, 3.11 and 3.12.
You can install the application:

- By cloning and installing this repository for the latest version
    ```{sh}
    git clone git@gitlab.inria.fr:postaviz/m2m-postaviz.git # or git clone https://gitlab.inria.fr/postaviz/m2m-postaviz.git
    ```
    Then install the tool:
    ```
    pip install .
    ```

    
- Directly from the last release on Python Pypi with `pip`
    ```
    pip install m2m-postaviz
    ```
    Or as an alternative you can also directly install the in-development version with:
    ```
    pip install git+ssh://git@gitlab.inria.fr/postaviz/m2m-postaviz.git@main
    ```


### Dependencies

<details>
  <summary>Click to expand</summary>
M2M-PostAViz has a few dependencies that are listed below:

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

You can install them with
```
pip install -r requirement.txt
``` 

They should be installed automatically when installing the application though. 

If you use the application for research, do not forget to cite the works associated to those dependencies.
</details>

## Quick start

> **⚠️ Warning:** 
> 
> We assume that you arrive at this step having installed the tool first (see above section), for instance in a Python virtual environment, or conda (mamba) environment.

📁 A dataset with test data is available in this repository in `postaviz_test_data` and can be used to test the main functionalities of the tool.


The test can be run with the following command:

```
m2m_postaviz --test
```

It takes a few seconds to launch because data needs to be uncompressed and processed in a temporary directory.

Shiny will launch automatically afterward.

<img src="./docs/pictures/postaviz_first_tab.png" alt="Homepage of Postaviz" width="70%"/>

Once on the homepage you're free to explore the test data.

Metacyc database is not include in test option.

The differents tabulation will be explain below.

##  Documentation

### Input data 

The mandatory input data are the outputs of Metage2Metabo for each sample/microbial community, and the metadata associated to each of them. Additional facultative inputs are advised to gain the most out of the analysis: taxonomy of the genomes associated to the metabolic networks, abundance of these genomes in the samples/community. It is also possible to provide the Metacyc ontology of the metabolic compounds to analyse the predictions at the level of metabolite families. The latter is only relevant if the metabolic networks were obtained with PathwayTools, i.e. are made of compound identifiers that fit the Metacyc database. 

> **💡 Note:** Metage2Metabo has a first pipeline step dedicated to the reconstruction of metabolic networks with Pathway Tools.
>
> If you used `m2m recon`, your metabolic networks are compatible with the Metacyc database and PostAViz can use the Metacyc ontology of compound families. 

In practice, other input data can be provided, included precomputed M2M-PostAViz tables which allow for a much faster restart when rerunning the app on previously analyszed data.

As a summary, for a first run you can provide all individual inputs (📄 below), but once this is processed, M2M-PostAViz can save the pre-processed data for a fast startover later (🚀 below)

We detail below the input data:

- 📄 M2M output for each sample
  - It should be in the following format

     ```
    │   ├── sample_1
    │   ├── community_analysis
    │   │   ├── addedvalue.json
    │   │   ├── comm_scopes.json
    │   │   ├── contributions_of_microbes.json
    │   │   ├── mincom.json
    │   │   ├── rev_cscope.json
    │   │   ├── rev_cscope.tsv
    │   │   └── targets.sbml
    │   ├── indiv_scopes
    │   │   ├── indiv_scopes.json
    │   │   ├── rev_iscope.json
    │   │   ├── rev_iscope.tsv
    │   │   └── seeds_in_indiv_scopes.json
    │   ├── m2m_metacom.log
    │   └── producibility_targets.json
    └──  sample_2
         ├── community_analysis
         │   ├── addedvalue.json
         │   ├── comm_scopes.json
         │   ├── contributions_of_microbes.json
         │   ├── mincom.json
         │   ├── rev_cscope.json
         │   ├── rev_cscope.tsv
         │   └── targets.sbml
         ├── indiv_scopes
         │   ├── indiv_scopes.json
         │   ├── rev_iscope.json
         │   ├── rev_iscope.tsv
         │   └── seeds_in_indiv_scopes.json
         ├── m2m_metacom.log
         └── producibility_targets.json
    ```
- 📄 Metadata associated to samples
  - This file should be a table where the first column is the sample identifier matching the output of M2M. For instance as below

    | smplID    | Age | Country  |
    |-----------|----:|----------|
    | sample_1  |  2  | France   |
    | sample_2  |  30 | Canada   |
    | sample_3  |  68 | Germany  |
  - The expected format is a tabulated file

- 📄 Taxonomy of the MAGs/genomes corresponding to the metabolic networks used in the analysis.
  - This file should be a table where the first column is the identifier matching the IDs of the metabolic networks that were analyses by M2M. For instance as below

    | Genome_ID  | Domain      | Phylum          | Class          | Order           | Family            | Genus      | Species        |
    |------------|------------|----------------|---------------|----------------|-------------------|------------|---------------|
    | MAG_1      | Bacteria   | Proteobacteria  | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia | Escherichia coli |
    | Genome_1   | Bacteria   | Firmicutes      | Bacilli       | Lactobacillales | Lactobacillaceae  | Lactobacillus | Lactobacillus casei |
    | MAG_2      | Archaea    | Euryarchaeota   | Methanobacteria | Methanobacteriales | Methanobacteriaceae | Methanobacterium | Methanobacterium formicicum |

  - The expected format is a tabulated file

- 📄 Abundance of the MAGs/genomes in the samples/communities.
  - Tabulated file. 
  - It will be **normalised by column sum** during processing.
  - Data is used to normalise the Boolean prediction of metabolite producibility (0 non producible, 1 producible) by the relative abundance of the compound producer.

    | identifier | Sample_1 | Sample_2 | Sample_3 |
    |------------|----------|----------|----------|
    | MAG_1      | 12.5     | 8.3      | 15.2     |
    | Genome_1   | 5.8      | 10.1     | 7.6      |
    | MAG_2      | 20.3     | 14.7     | 18.9     |

  - The expected format is a tabulated file
  
- 📄 Metacyc database in padmet format
  - This optional input file is used to take advantage of the ontology of compound categories provided in the Metacyc database. 
  - This input is relevant only if the metabolic networks analysed in M2M were generated with the PathwayTools software.
  - The Metacyc database flat files can be downloaded provided subscription. Once downloaded, they can be integrated into a single file in the padmet format using the following command line 
  ```{sh}
  # install padmet
  pip install padmet
  # move to the directory were the Metacyc data is stored (unzipped). The data is in a directory named by the version of the database, for instance `28.5/data`
  padmet pgdb_to_padmet --pgdb 28.5/data --version 28.5 --db metacyc --output metacyc28_5.padmet
  ```

- 🚀 Precomputed data for M2M-PostAViz
  - #TODO This data can be stored when running the tool with the `-o` flag. It will be saved in the directory of your choice and can be loaded for future runs of M2M-PostAViz.

  ```{sh}
    # Run the tool once on the data and provide a path for saving the tables
    m2m_postaviz -d Metage2metabo/samples/scopes/directory/path
                -m metadata/file/path
                -a abundance/file/path
                -t taxonomy/file/path
                -o save/directory/path
    # for future runs, if data has not changed you can simply resume exploration with
    m2m_postaviz -l save/directory/path
    ```


<!-- ### Important

In metadata file the first column must be the sample identification. Preferably named "smplID".

In taxonomy file the first column must be the metagenomes (mgs). Preferably named "mgs". -->

### Usage

Based on the input listed above, `m2m_postaviz` can be run in two ways:

- 📄 🐢 by providing all input data. To avoid doing this at each run -- pre-processing of all data by the application can be quite lengthy if many samples are provided --, we advise users to use the `-o` flag and save the precomputed data for future runs. In that case, users can resume the analysis directly by loading the processed data (see item below). This procedure is needed at least once, for the first run with new datasets, or when datasets are altered. Note that metadata changes do not need to re-run the whole preprocessing: you can directly modify the file in the saved directory.
  
    ```
    m2m_postaviz -d Metage2metabo/samples/scopes/directory/path
                -m metadata/file/path
                -a abundance/file/path
                -t taxonomy/file/path
                -o save/path
                --no-metacyc (Optionnal)
    ```
- 🚀 by providing the preprocessed data

    ```
        m2m_postaviz -l save/directory/path
    ````

<!-- This way is required as least one time to produce all dataframe and save them in -o save/path.

Once the dataframes are produced. Shiny will automatically run from the save/path given in -o option.
You can interrupt the process if you want and run postaviz with -l load option.



Which will directly launch shiny and skip dataprocessing. -->

> **💡 Note:** The preprocessed dataset is stored in a directory in the form of dataframes and xxxxxx Parquet format. This format enables an efficient storage and data access by the application
> 
> Below is the structure of the preprocessed directory.
>
> #TODO add a tree of the preprocessed dir

## Application presentation

The application starts in a web browser and enables user to analyse metabolic potential predictions in the light of sample metadata, genome taxonomy and possibly taking genome abundance into account to weight producibility of metabolites. Several tabs dedicated to different analyses can be browsed. 

Users can modify the visualisations and the statistical analyses by selecting and filtering data and metadata. 

We detail below the contents of each tab and the analyses it enables to perform. 

### Metadata tab

Tabulation to observe the metadata given in CLI.

From this tab the metadata's column dtype can be changed.

This can be helpfull when you are not certain of your metadata types in CLI.

Sometime Plotly and seaborn do not treat numeric / not numeric columns the same when building plot's axes.

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

> **⚠️ Warning:**
>The "all" option on all sample (No metadata filter applied) can be long to produce the plots. Also heavy plots will impact the performance of the application. 


> **💡 Note:** A small text output under the Processing button show how many bins are selected to avoid large calculation. Also if only one bin (mgs) is selected, it will display how many samples have this specific bin.
e.

#### Method

The bins exploration present some challenges, it need to retain the production by the bins inside of each samples.

The individual production of each bins depend of the bins present (community interaction highlighted by **Meta2metabo**) which also differ from the treatment/condition of the sample.

All of this information can scale poorly if a lot of sample are present in input which is the point of the application: being able observe from lots of angles large chunks of data.

In order to keep all this valuable data we choosed to use the hard drive and load/manipulate large dataframe by subset using Parquet format. That way RAM memory is not overload, and we can pick what we need from the hard drive, using query similar to a SQL database.


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


- Enable row/columns clustering (Only for heatmap) will change column and or row order. Optional and independent from each other.

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