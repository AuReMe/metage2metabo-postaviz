
Changelog
=========

1.0.0 (2026-11-20)
------------------

First stable release of M2M-PostAViz. This version is the one associated with
the publication of M2M-PostAViz in the Journal of Open Source Software (JOSS).

Licensing
^^^^^^^^^

* The project is now distributed under the GNU Affero General Public License v3
  (AGPLv3), replacing the GNU Lesser General Public License v3 (LGPLv3).

Application and analyses
^^^^^^^^^^^^^^^^^^^^^^^^

* Reworked statistical analyses, including multiple testing correction.
* Enhanced metabolite ontology management for the compound exploration, with the
  addition of an NPClassifier-based compound hierarchy.
* Extended taxonomy exploration: grouping of samples by metadata variable,
  additional plots, and export of the underlying data frames.
* Improved handling of datasets provided without abundance or taxonomic data.

Command-line interface and data management
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* The CLI supports an interactive mode to guide the selection of input files.
* The output directory is now a mandatory argument: precomputed data are stored
  on disk in the parquet format and reused in subsequent sessions via ``--load``.

Packaging, testing and documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Support and continuous integration for Python 3.12, 3.13 and 3.14, on Ubuntu
  and macOS, with scheduled runs.
* Dependency constraints adjusted to avoid failures caused by pyarrow and
  python-libsbml releases.
* Expanded documentation: tutorial with figures, installation troubleshooting,
  clarified application startup, and the origin of the bundled data files.

0.0.0 (2023-12-01)
------------------

* First release on PyPI.
