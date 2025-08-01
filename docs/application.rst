Application details
===================

The application starts in a web browser and enables users to analyse metabolic potential predictions in the light of sample metadata, genome taxonomy, and possibly taking genome abundance into account to weight producibility of metabolites. Several tabs dedicated to different analyses can be browsed.

Users can modify the visualisations and the statistical analyses by selecting and filtering data and metadata.

Tabs overview
-------------

The application is divided into several tabs, each dedicated to a specific type of analysis or data management:

*Overview* tab
~~~~~~~~~~~~

This tab summarises the dataset and permits a first exploration of it. 
A first summary presents the number of microorganism metabolic networks considered, the number of samples and the number of unique metabolites reached overall across the samples.

The tab displays several plots:

- **Number of unique metabolites reached in each community / sample.**  
  Metabolites can be either reached by individual community members or reached through interactions by the community. This plot can be refined by grouping the data using metadata variables.

- **Total number of compounds reached in samples.**  
  This plot is entirely customisable and aims at comparing groups of samples. It can take into account the relative abundance of each microorganism producer to weigh predictions of producibility instead of using {0,1} values (0: unproducible, 1: producible). Statistical results comparing metadata groups are provided as a table. Multiple test correction can be applied to the statistical tests. Note that statistical tests depend on the type of the variables. You can change those types in the Metadata Management tab if needed.

- **Principal Coordinate Analysis (PCoA) of compound production in samples.**  
  Two plots are provided. The first one is computed on the whole dataset, although one may show/hide certain samples according to metadata variables to ease visualisation. Initial matrix used for computation of this plot is the producibility (1) or non-producibility (0) of metabolites in samples considering that producibility may result from cross-feeding (metabolic complementarity) between populations. This relates to the so-called *community scopes* of Metage2Metabo.
  
  In the second PCoA plot, it is possible to filter data and perform de novo PCoA on groups of samples of interest, taking into account the relative abundance of taxa to weigh producibility or not. 

All data frames used for plot creation, as well as statistical tests and visualisations, can be exported by users.

*Taxonomy-based exploration* tab
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tab is dedicated to the exploration of metabolic networks matching the genomes (or bins) included in the samples. The tab proposes analyses performed at the scale of individual scopes (metabolites produced by individual populations, not considering microbial interactions) and community scopes (metabolites produced by the community, that may or may not depend on metabolic interactions).

Users can check the produced compounds for multiple taxonomic levels, according to the taxonomic metadata that is provided as input. In particular, it is possible to consider taxonomic groups, in which case, all metabolic networks corresponding to this group will be taken into account.

Customization of the plot:

- Choice of the taxonomic ranks: all, taxonomic groups (e.g. phylum), or individual genomes.
- Selection within the rank, e.g. one phylum of interest if phylum was selected above. The number of metabolic networks (bins) corresponding to the selection is displayed at the bottom of the box.
- "Filter" enables selecting samples based on a variable of interest from metadata. "Select" restricts the analysis on the samples matching the metadata variable values that are selected.
- "Color" permits grouping the samples on a metadata variable.

If checking the box "weigh the producibility value with relative abundance", the normalised abundance data frame will be used for computation instead of the main dataframe. This means that producibility values of compounds by microorganisms will no longer be 0 or 1 but values corresponding to the relative abundance of the microbe. For instance, 0.3 for a producer that represents 30% of the community's relative abundance.

Two plots are generated showing the number of metabolites selected taxa can produce:


- **Production of unique metabolites by taxa in samples.**  
  This plot shows the number of unique metabolites produced by the selected taxa in each sample. The number of metabolites is computed either considering the community scope (i.e. interactions between taxa) or the individual scope (i.e. no interactions). In the first case, selected taxa are the producers but they may depend on alternative populations that lead to the production of metabolite precursors. The second case considers only the taxa selected, meaning that they are the only producers of metabolites in the sample. The plot can be grouped by metadata variables.
- **Relative abundance of taxa in samples.**  
  This plot shows the relative abundance of the selected taxa in each sample. It can be grouped by metadata variables. The relative abundance is computed based on the normalised abundance data frame, which is generated by the application if relative abundance data is provided as input. If not, the plot will not be generated.

.. note::
   The "all" option on all samples (no metadata filter applied) can lead to long computational times when working with large datasets. Performances of the application may be impacted.

.. note::
   A text output under the "Go" button shows how many bins are selected to anticipate large calculations. Also, if only one metabolic network (genome, bin or MGS) is selected, the number of samples the bin is present in will be indicated.

*Metabolite-based exploration* tab
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tab is dedicated to the exploration of metabolites produced by the metabolic networks of microorganisms present in the samples. It allows users to explore the producibility of certain metabolites at the level of individual compounds or compound families.

Plot settings depend on whether Metacyc data containing the ontology of metabolite families is used for computation of results with Metage2Metabo and provided as input (see :doc:`input_data_details`) to the application or not. If Metacyc data is enabled, each metabolite is associated to an ontology of metabolite families, enabling to explore producibility at the level of these families.

Users may select groups of compounds based on the list of Metacyc categories ordered from the top to the bottom of the tree. Any category selected will update the following selection menu to a list of all sub-categories that can in turn be selected. The following field is automatically updated with compounds matching the (sub)-categories.

Without Metacyc input data, users may directly select their compounds of interest in this third field.

The plots generated will only consider the compounds selected as input.

- **Metadata filter and color**
    - Metadata filter
    - Plot color and regroup

- **Sample filtering** enables to target groups of samples of interest based on metadata variables.
    - Pick *all* (no filter, by default), *include* or *exclude*
    - Select a metadata variable
    - Select the values of the metadata variable to exclude or include.
    - Matching samples are automatically filled. You can further refine the selection by deleting some in the list.

- **Additional options:**
    - **Enable relative abundance weighting** will use the normalised abundance data frame to compute the producibility of metabolites instead of the main dataframe. This means that producibility values of compounds by microorganisms will no longer be 0 or 1 but values corresponding to the relative abundance of the microbe. For instance, 0.3 for a producer that represents 30% of the community's relative abundance.
    - **Add row/columns clustering** in the heatmaps. This will change column and/or row order based on producibility value similarity.
    - **Generate statistical dataframe** to compare groups of samples related to the boxplot of the tab.

Plots:

- **Heatmap:**  
  Heatmap displaying the number of bins producing the compound in the sample. Can be cusomtised with the options above: relative abundance weighting, sample filtering, row/column clustering. Three heatmaps are generated, each in a tab:
  - Community metabolic potential: takes into account cross-feeding and metabolic interactions between taxa.
  - Individual metabolic potential: does not take into account cross-feeding and metabolic interactions between taxa.
  - Added value: metabolites that are produced by the community but not by any of the individual taxa in the samples.

- **Percentage of samples producing selected compounds by groups:**  
  At least one producer in the sample has to produce the metabolite for producibility to be ensured. Two tabs are generated: community metabolic potential and individual metabolic potential. The percentage of samples producing the selected compounds is displayed, grouped by metadata variables.  

- **Boxplot of the number of producers by samples grouped by metadata:**
  Two metadata variables can be selected to group the samples. The boxplot shows the number of producers of the selected compounds in each sample, grouped by the metadata variables. A statistical test is performed to compare the groups of samples. The statistical test used depends on the type of metadata variable selected. For instance, if a categorical variable is selected, a Wilcoxon test is performed to compare the groups of samples. Results of the statistical tests are displayed in a table below the plot. Tested pairs are determined by the metadata input. Sample filtering is not applied here.

*Metadata management* tab
~~~~~~~~~~~~~~~~~~~~~~

This tab summarises the metadata that was provided as input to the tool and the data type that were associated to each variable. It is possible to change the type of each variable, which will impact the statistical tests performed in the application. For instance, if a variable is set as categorical, a Wilcoxon test will be performed to compare groups of samples. If it is set as numeric, a correlation test will be performed instead. Likewise, the plots generated across the application will adapt to the type of the variable.

Sometimes Plotly and seaborn do not treat numeric / non-numeric columns similarly when building plot's axes. You may want to pick *str* over *category* for factor variables. 

.. note::
   All data frames used for plot creation, as well as statistical tests and visualisations, can be exported by users.
