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

from m2m_postaviz import __version__ as VERSION
import m2m_postaviz.data_utils as du
import m2m_postaviz.shiny_app as sh
from pathlib import Path
from m2m_postaviz.data_struct import DataStorage

LICENSE = """Copyright (C) L. brindel and C. Frioux.\n
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.\n

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.\n

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>\n
"""
MESSAGE = """
M2M-postAViz is an application dedicated to the visualization of Metage2Metabo's metabolic complementarity results applied to multiple community compositions / samples. Type `m2m_postaviz --help` for help and see the documentation at https://metage2metabo-postaviz.readthedocs.io/en/latest/ for more information.
"""

parser = argparse.ArgumentParser("m2m_postaviz", description=MESSAGE + "\n\n" + LICENSE, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-d", "--dir", help="Directory containing the data")
parser.add_argument("-m", "--metadata", help="Tsv file containing metadata")
parser.add_argument("-t", "--taxonomy", help="Tsv file containing taxonomy data")
parser.add_argument("-a", "--abundance", help="Abundance data file as tsv.")
parser.add_argument("-o", "--output", help="Output path for saved plot of dataframe. If created if not valid. If not provided, save options disabled.")
parser.add_argument("-l", "--load", help="Run postaviz from save directory")
parser.add_argument("-c", "--metacyc", help="Run postaviz with the metacyc database as padmet file. This is usefull when the metabolite ID from the scopes use metacyc ID. Enable the research by category of metabolites.")

parser.add_argument("-v", "--version", action="version", 
                    help="Show the version of m2m-postaviz", version="%(prog)s " + VERSION + "\n" + LICENSE)



parser.add_argument("--test", help="Run postaviz with test files only", action="store_true")

SRC_DIR = Path(__file__).parent
PROJECT_DIR = Path(SRC_DIR).parent
TESTS_DIR = Path(SRC_DIR, "postaviz_test_data/")

data_table_filepath = Path(TESTS_DIR, "table_test_postaviz.tar.gz")

def main(args=None):
    arg_parser = parser.parse_args()

    if arg_parser.load:

        Data = DataStorage(Path(arg_parser.load))
        sh.run_shiny(Data)

    elif arg_parser.test:

        if not Path(TESTS_DIR, "palleja").is_dir():
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
        
        # Check if required arguments are provided
        if not arg_parser["dir"] or not arg_parser["metadata"]:
            print(MESSAGE)
            print("\nError: Required arguments missing.")
            print("Use --help for more information about required arguments.")
            return
        
        dir_path = Path(arg_parser["dir"]).resolve()
        metadata_path = Path(arg_parser["metadata"]).resolve()

        try:
          taxonomic_path = Path(arg_parser["taxonomy"]).resolve()
        except:
          taxonomic_path = None
        try:
          abundance_path = Path(arg_parser["abundance"]).resolve()
        except:
          abundance_path = None
          
        if not arg_parser["output"]:
            print("Warning: No output path provided. Save options will be disabled.")
            save_path = None
        else:
            save_path = Path(arg_parser["output"]).resolve()
        
        try:
          metacyc = Path(arg_parser["metacyc"]).resolve()
        except:
          metacyc = None

        du.build_dataframes(dir_path, metadata_path, abundance_path, taxonomic_path, save_path, metacyc)

        Data = DataStorage(save_path)

        sh.run_shiny(Data)
