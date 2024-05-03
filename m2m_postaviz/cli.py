"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mm2m_postaviz` python will execute
    ``__main__.py`` as a script. That means there will not be any
    ``m2m_postaviz.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there"s no ``m2m_postaviz.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import argparse
import m2m_postaviz.data_utils as du
import m2m_postaviz.shiny_app as sh
from m2m_postaviz.data_struct import DataStorage
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="Directory containing the data")
parser.add_argument("-m", "--metadata", help="Tsv file containing metadata")
parser.add_argument("-t", "--taxonomy", help="Tsv file containing taxonomy data")
parser.add_argument("--test", help="Run postaviz with test files only", action="store_true")
parser.add_argument("--dev", help="Run postaviz for dev only", action="store_true")
parser.add_argument("-ut", help="Run postaviz unit test for dev only", action="store_true")

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(SRC_DIR))
TESTS_DIR = os.path.join(SRC_DIR, 'postaviz_test_data/')

data_table_filepath = os.path.join(TESTS_DIR, 'table_test_postaviz.tar.gz')


def main(args=None):
    arg_parser = parser.parse_args()

    if arg_parser.ut:
      du.unit_test_1()
      quit()
    if arg_parser.test:
      if not os.path.isdir(os.path.join(TESTS_DIR, 'data_test/')):
        print("No data_test/ directory found.")
        du.extract_tarfile(data_table_filepath, TESTS_DIR)
      data_test_dir = os.path.join(TESTS_DIR, 'data_test/')
      global_data, sample_data, abundance_data = du.build_test_data(data_test_dir)
      taxonomic_data = du.open_tsv(TESTS_DIR+"taxonomic_database.tsv")

    elif arg_parser.dev:
      # dir_path = TESTS_DIR+"metadata_ouput/"
      dir_path = "/home/lbrindel/output/palleja/"
      metadata = TESTS_DIR+"refined_palleja_metadata.tsv"
      abundance_path = "~/Downloads/matrix_palleja.tsv"
      taxonomic_data = du.open_tsv("~/Downloads/gtdbtk.summary_split.tsv")
      global_data, sample_data, norm_abundance_data, abundance_data = du.build_df(dir_path, metadata, abundance_path)

    else:
      arg_parser = vars(parser.parse_args())
      dir_path = arg_parser["dir"]
      metadata = arg_parser["metadata"]
      taxonomy = arg_parser["taxonomy"]
      taxonomic_data = du.open_tsv(taxonomy)
      global_data, sample_data, abundance_data = du.build_df(dir_path, metadata, abundance_path)
    
    Data = DataStorage(global_data, sample_data, taxonomic_data, norm_abundance_data, abundance_data)
    # Data.performance_benchmark()

    sh.run_shiny(Data)
