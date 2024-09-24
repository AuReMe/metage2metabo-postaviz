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
import os

import m2m_postaviz.data_utils as du
import m2m_postaviz.shiny_app as sh
from m2m_postaviz.data_struct import DataStorage

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="Directory containing the data")
parser.add_argument("-m", "--metadata", help="Tsv file containing metadata")
parser.add_argument("-t", "--taxonomy", help="Tsv file containing taxonomy data")
parser.add_argument("-a", "--abundance", help="Abundance data file as tsv.")
parser.add_argument("-o", "--output", help="Output path for saved plot of dataframe. If created if not valid. If not provided, save options disabled.")

parser.add_argument("--test", help="Run postaviz with test files only", action="store_true")
parser.add_argument("--dev", help="Run postaviz for dev only", action="store_true")

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(SRC_DIR))
TESTS_DIR = os.path.join(SRC_DIR, "postaviz_test_data/")

data_table_filepath = os.path.join(TESTS_DIR, "table_test_postaviz.tar.gz")


def main(args=None):
    arg_parser = parser.parse_args()

    if arg_parser.test:
        
        if not os.path.isdir(os.path.join(TESTS_DIR, "palleja/")):
            print("No data_test/ directory found. \nExtract test data tarfile...")
            du.extract_tarfile(data_table_filepath, TESTS_DIR)

        data_test_dir = os.path.join(TESTS_DIR, "palleja/")
        metadata_path = os.path.join(data_test_dir, "metadata_test_data.tsv")
        abundance_path = os.path.join(data_test_dir, "abundance_test_data.tsv")
        taxonomy_path = os.path.join(data_test_dir, "taxonomy_test_data.tsv")
        save_path = "test_path"
        file_format, hdf5_file_path, taxonomy_provided, abundance_provided = du.build_df(data_test_dir, metadata_path, abundance_path, taxonomy_path, save_path)

    elif arg_parser.dev:
        
        dir_path = "/home/lbrindel/output/western_diet_samples/res_smpl1/"
        # dir_path = "/home/lbrindel/output/western_diet_samples/all_samples/"
        
        metadata_path = "~/Downloads/western_diet_exp/metadata_drama.tsv"
        abundance_path = "~/Downloads/western_diet_exp/specI.mat"
        taxonomic_path = "~/Downloads/western_diet_exp/taxonomies.tsv"

        # save_path = "/home/lbrindel/output/TEST_BUILD/"
        save_path = "/home/lbrindel/output/test_res_smpl1/"
        # save_path = "/home/lbrindel/output/full_run_postaviz/"
        
        file_format, hdf5_file_path, taxonomy_provided, abundance_provided = du.build_df(dir_path, metadata_path, abundance_path, taxonomic_path, save_path)

    else:
        
        arg_parser = vars(parser.parse_args())
        dir_path = arg_parser["dir"]
        metadata_path = arg_parser["metadata"]
        taxonomy_path = arg_parser["taxonomy"]
        abundance_path = arg_parser["abundance"]
        save_path = arg_parser["output"]
        file_format, hdf5_file_path, taxonomy_provided, abundance_provided = du.build_df(dir_path, metadata_path, abundance_path, taxonomic_path, save_path)

    Data = DataStorage(save_path, file_format, hdf5_file_path, taxonomy_provided, abundance_provided)

    sh.run_shiny(Data)
