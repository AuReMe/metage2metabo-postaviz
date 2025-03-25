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
import tempfile

import m2m_postaviz.data_utils as du
import m2m_postaviz.shiny_app as sh
from pathlib import Path
from m2m_postaviz.data_struct import DataStorage

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="Directory containing the data")
parser.add_argument("-m", "--metadata", help="Tsv file containing metadata")
parser.add_argument("-t", "--taxonomy", help="Tsv file containing taxonomy data")
parser.add_argument("-a", "--abundance", help="Abundance data file as tsv.")
parser.add_argument("-o", "--output", help="Output path for saved plot of dataframe. If created if not valid. If not provided, save options disabled.")
parser.add_argument("-l", "--load", help="Run postaviz from save directory")
parser.add_argument("-c", "--metacyc", help="Run postaviz with the metacyc database as padmet file. This is usefull when the metabolite ID from the scopes use metacyc ID. Enable the research by category of metabolites.")

parser.add_argument("--test", help="Run postaviz with test files only", action="store_true")

SRC_DIR = Path(__file__).parent
PROJECT_DIR = Path(SRC_DIR).parent
TESTS_DIR = Path(SRC_DIR, "postaviz_test_data/")

data_table_filepath = Path(TESTS_DIR, "table_test_postaviz.tar.gz")

def main(args=None):
    arg_parser = parser.parse_args()

    if arg_parser.load:

        Data = DataStorage(arg_parser.load)
        sh.run_shiny(Data)

    elif arg_parser.test:

        if not Path(TESTS_DIR, "palleja/").is_dir:
            print("No data_test/ directory found. \nExtract test data tarfile...")
            du.extract_tarfile(data_table_filepath, TESTS_DIR)

        data_test_dir = Path(TESTS_DIR, "palleja/")
        metadata_path = Path(data_test_dir, "metadata_test_data.tsv")
        abundance_path = Path(data_test_dir, "abundance_test_data.tsv")
        taxonomy_path = Path(data_test_dir, "taxonomy_test_data.tsv")
        tempdir = tempfile.TemporaryDirectory()
        save_path = Path(tempdir.name)

        du.build_dataframes(data_test_dir, metadata_path, abundance_path, taxonomy_path, save_path)
        
        Data = DataStorage(save_path)

        sh.run_shiny(Data)

        tempdir.cleanup()

    else:

        arg_parser = vars(parser.parse_args())
        dir_path = Path(arg_parser["dir"]).resolve()
        metadata_path = Path(arg_parser["metadata"]).resolve()
        taxonomic_path = Path(arg_parser["taxonomy"]).resolve()
        abundance_path = Path(arg_parser["abundance"]).resolve()
        save_path = Path(arg_parser["output"]).resolve()
        metacyc = Path(arg_parser["metacyc"]).resolve()
        du.build_dataframes(dir_path, metadata_path, abundance_path, taxonomic_path, save_path, metacyc)

        Data = DataStorage(save_path)

        sh.run_shiny(Data)
