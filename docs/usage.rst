Usage
=====

Python usage
------------

To use m2m-postaviz in a project::

    import m2m_postaviz

Command-line usage
------------------

Based on the input listed in :doc:`input_data`, ``m2m_postaviz`` can be run in two ways:

- **By providing all input data (slow mode)** (first run or when datasets change):

  .. code-block:: bash

      m2m_postaviz -d Metage2metabo/samples/scopes/directory/path \
                   -m metadata/file/path \
                   -a abundance/file/path \
                   -t taxonomy/file/path \
                   -o save/path \
                   --no-metacyc  # (Optional)

- **By providing the preprocessed data (fast mode)** (for fast restarts):

  .. code-block:: bash

      m2m_postaviz -l save/directory/path

.. note::
   The preprocessed dataset is stored in a directory in the form of dataframes and files in Parquet format for efficient storage and access.

For more details on input data and directory structure, see below.

Input Data
==========

Summary of input files
----------------------

+-------------------+-------------------------------------------------------------+
| File              | Description                                                 |
+===================+=============================================================+
| M2M output        | Output directory for each sample from Metage2Metabo         |
+-------------------+-------------------------------------------------------------+
| Metadata          | Tabulated file, first column is sample identifier           |
+-------------------+-------------------------------------------------------------+
| Taxonomy          | Tabulated file, first column is genome/metabolic network ID |
+-------------------+-------------------------------------------------------------+
| Abundance         | Tabulated file, normalized by column sum                    |
+-------------------+-------------------------------------------------------------+
| Metacyc (optional)| Padmet format, for compound ontology                        |
+-------------------+-------------------------------------------------------------+
| Precomputed       | Directory with preprocessed dataframes (for fast restart)   |
+-------------------+-------------------------------------------------------------+

The mandatory input data are the outputs of Metage2Metabo for each sample/microbial community, and the metadata associated to each of them. Additional facultative inputs are advised to gain the most out of the analysis: taxonomy of the genomes associated to the metabolic networks, abundance of these genomes in the samples/community. It is also possible to provide the Metacyc ontology of the metabolic compounds to analyse the predictions at the level of metabolite families. The latter is only relevant if the metabolic networks were obtained with PathwayTools, i.e. are made of compound identifiers that fit the Metacyc database.

.. note::
   Metage2Metabo has a first pipeline step dedicated to the reconstruction of metabolic networks with Pathway Tools.
   If you used ``m2m recon``, your metabolic networks are compatible with the Metacyc database and PostAViz can use the Metacyc ontology of compound families.

In practice, other input data can be provided, including precomputed M2M-PostAViz tables which allow for a much faster restart when rerunning the app on previously analysed data.

.. include:: input_data_details.rst
