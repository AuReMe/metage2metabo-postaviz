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

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="Directory containing the data")
parser.add_argument("-m", "--metadata", help="Tsv file containing metadata")


def main(args=None):
    arg_parser = vars(parser.parse_args())
    dir_path = arg_parser["dir"]
    metadata = arg_parser["metadata"]
    main_df = du.build_df(dir_path, metadata)
    data = du.open_tsv(metadata)
    sh.run_shiny(main_df, data)
