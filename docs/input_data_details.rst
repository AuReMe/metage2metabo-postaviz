.. This file is included from input_data.rst

Input data details
------------------

- **M2M output for each sample**:  
  Directory structure example:

  .. code-block:: text

      sample_1/
        community_analysis/
          addedvalue.json
          comm_scopes.json
          contributions_of_microbes.json
          mincom.json
          rev_cscope.json
          rev_cscope.tsv
          targets.sbml
        indiv_scopes/
          indiv_scopes.json
          rev_iscope.json
          rev_iscope.tsv
          seeds_in_indiv_scopes.json
        m2m_metacom.log
        producibility_targets.json
      sample_2/
        ...

- **📄 Metadata associated to samples**:  
  Tabulated file, first column is the sample identifier matching the output of M2M.

  +----------+-----+---------+
  | smplID   | Age | Country |
  +==========+=====+=========+
  | sample_1 |  2  | France  |
  +----------+-----+---------+
  | sample_2 | 30  | Canada  |
  +----------+-----+---------+
  | sample_3 | 68  | Germany |
  +----------+-----+---------+

- **📄 Taxonomy of the MAGs/genomes**:  
  Tabulated file, first column matches the IDs of the metabolic networks.

  +-----------+----------+-------------------+---------------------+---------------------+-----------------------+----------------+--------------------------+
  | Genome_ID | Domain   | Phylum            | Class               | Order               | Family                | Genus          | Species                  |
  +===========+==========+===================+=====================+=====================+=======================+================+==========================+
  | MAG_1     | Bacteria | Proteobacteria    | Gammaproteobacteria | Enterobacterales    | Enterobacteriaceae    | Escherichia    | Escherichia coli         |
  +-----------+----------+-------------------+---------------------+---------------------+-----------------------+----------------+--------------------------+
  | Genome_1  | Bacteria | Firmicutes        | Bacilli             | Lactobacillales     | Lactobacillaceae      | Lactobacillus  | Lactobacillus casei      |
  +-----------+----------+-------------------+---------------------+---------------------+-----------------------+----------------+--------------------------+
  | MAG_2     | Archaea  | Euryarchaeota     | Methanobacteria     | Methanobacteriales  | Methanobacteriaceae   | Methanobacterium | Methanobacterium formicicum |
  +-----------+----------+-------------------+---------------------+---------------------+-----------------------+----------------+--------------------------+

- **📊 Abundance of the MAGs/genomes in the samples/communities**:  
  Tabulated file, normalized by column sum during processing.

  +------------+----------+----------+----------+
  | identifier | Sample_1 | Sample_2 | Sample_3 |
  +============+==========+==========+==========+
  | MAG_1      |   12.5   |   8.3    |   15.2   |
  +------------+----------+----------+----------+
  | Genome_1   |   5.8    |  10.1    |   7.6    |
  +------------+----------+----------+----------+
  | MAG_2      |   20.3   |  14.7    |  18.9    |
  +------------+----------+----------+----------+

- **🚀 Precomputed data for M2M-PostAViz**:  
  Can be stored when running the tool with the ``-o`` flag and loaded for future runs.

  .. code-block:: bash

      m2m_postaviz -d Metage2metabo/samples/scopes/directory/path \
                   -m metadata/file/path \
                   -a abundance/file/path \
                   -t taxonomy/file/path \
                   -o save/directory/path

      # For future runs:
      m2m_postaviz -l save/directory/path

  The preprocessed dataset is stored in a directory in the form of dataframes and files in Parquet format. Example structure:

  .. code-block:: text

      saved_data_postaviz/
        abundance_file_normalised.tsv
        abundance_file.tsv
        ...
        sample_cscope_directory/
          Sample1.parquet.gzip
          ...
        sample_iscope_directory/
          Sample1.parquet.gzip
          ...
        ...
