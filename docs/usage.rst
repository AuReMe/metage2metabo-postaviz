=====
Usage
=====

Python usage
------------

To use m2m-postaviz in a project::

    import m2m_postaviz

Command-line usage
------------------

Based on the input listed in :doc:`input_data`, ``m2m_postaviz`` can be run in two ways:

- **By providing all input data** (first run or when datasets change):

  .. code-block:: bash

      m2m_postaviz -d Metage2metabo/samples/scopes/directory/path \
                   -m metadata/file/path \
                   -a abundance/file/path \
                   -t taxonomy/file/path \
                   -o save/path \
                   --no-metacyc  # (Optional)

- **By providing the preprocessed data** (for fast restarts):

  .. code-block:: bash

      m2m_postaviz -l save/directory/path

.. note::
   The preprocessed dataset is stored in a directory in the form of dataframes and files in Parquet format for efficient storage and access.

For more details on input data and directory structure, see :doc:`input_data`.
