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
import sys
import tempfile
from typing import Optional

from m2m_postaviz import __version__ as VERSION
import m2m_postaviz.data_utils as du
from pathlib import Path
from m2m_postaviz.data_struct import DataStorage

try:
    from prompt_toolkit import prompt as pt_prompt
    from prompt_toolkit.completion import PathCompleter
except ModuleNotFoundError:
    pt_prompt = None
    PathCompleter = None

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
parser.add_argument(
    "-o",
    "--output",
    help="Output directory where precomputed PostAViz data are saved for reuse with --load, and where user-exported plots/data are written. Must be provided.",
)
parser.add_argument("-l", "--load", help="Run postaviz from save directory")
parser.add_argument("-c", "--metacyc", help="Run postaviz with the metacyc database as padmet file. This is usefull when the metabolite ID from the scopes use metacyc ID. Enable the research by category of metabolites.")

parser.add_argument("-v", "--version", action="version", 
                    help="Show the version of m2m-postaviz", version="%(prog)s " + VERSION + "\n" + LICENSE)



parser.add_argument("--test", help="Run postaviz with test files only", action="store_true")

SRC_DIR = Path(__file__).parent
PROJECT_DIR = Path(SRC_DIR).parent
TESTS_DIR = Path(SRC_DIR, "postaviz_test_data/")

data_table_filepath = Path(TESTS_DIR, "table_test_postaviz.tar.gz")


class InteractivePromptExit(Exception):
    """Raised when user exits interactive prompt."""


def _is_interactive_terminal() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def _check_exit_command(answer: str) -> None:
    """Raise InteractivePromptExit if user typed exit command."""
    if answer.lower() in {"quit", "exit", "q"}:
        raise InteractivePromptExit("User exited interactive prompt.")


def _read_path_input(label: str, only_directories: bool = False) -> str:
    # Use shell-like path completion when prompt_toolkit is available.
    hint = " (or 'quit'/'exit' to cancel)"
    full_label = label + hint
    if pt_prompt is not None and PathCompleter is not None:
        completer = PathCompleter(only_directories=only_directories, expanduser=True)
        return pt_prompt(f"{full_label}: ", completer=completer).strip()
    return input(f"{full_label}: ").strip()


def _prompt_yes_no(label: str, default: bool = False) -> bool:
    suffix = " [Y/n]" if default else " [y/N]"
    hint = " (or 'quit'/'exit' to cancel)"
    answer = input(f"{label}{suffix}{hint}: ").strip()
    _check_exit_command(answer)
    answer_lower = answer.lower()
    if not answer_lower:
        return default
    return answer_lower in {"y", "yes"}


def _prompt_required_path(label: str, must_exist: bool = True, is_dir: bool = False) -> Path:
    while True:
        answer = _read_path_input(label, only_directories=is_dir)
        _check_exit_command(answer)
        if not answer:
            print("Value cannot be empty.")
            continue

        path = Path(answer).resolve()

        if not must_exist:
            return path

        if not path.exists():
            print(f"Path does not exist: {path}")
            continue

        if is_dir and not path.is_dir():
            print(f"Expected a directory: {path}")
            continue

        if not is_dir and not path.is_file():
            print(f"Expected a file: {path}")
            continue

        return path


def _prompt_optional_path(label: str, expect_file: bool = True) -> Optional[Path]:
    while True:
        answer = _read_path_input(label, only_directories=not expect_file)
        _check_exit_command(answer)
        if not answer or answer.lower() in {"skip", "none"}:
            return None

        path = Path(answer).resolve()

        if not path.exists():
            print(f"Path does not exist: {path}. Press Enter to skip.")
            continue

        if expect_file and not path.is_file():
            print(f"Expected a file: {path}. Press Enter to skip.")
            continue

        if not expect_file and not path.is_dir():
            print(f"Expected a directory: {path}. Press Enter to skip.")
            continue

        return path


def _optional_path(value, expect_file: bool = True):
    if not value:
        return None
    path = Path(value).resolve()
    if expect_file and not path.is_file():
        print(f"Warning: file not found, ignoring: {path}")
        return None
    return path


def _run_shiny(data: DataStorage) -> None:
    # Import lazily to avoid importing heavy UI deps during CLI module import.
    import m2m_postaviz.shiny_app as sh

    sh.run_shiny(data)

def main(args=None):
    try:
        arg_parser = parser.parse_args(args)
    except KeyboardInterrupt:
        print("\nExiting...")
        return

    try:
        if arg_parser.load:
            Data = DataStorage(Path(arg_parser.load))
            _run_shiny(Data)

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

            _run_shiny(Data)

            tempdir.cleanup()

        else:
            # If required arguments are missing and we are in a terminal, ask interactively.
            if not arg_parser.dir or not arg_parser.metadata or not arg_parser.output:
                if _is_interactive_terminal():
                    print("Missing required arguments. Interactive mode enabled.")
                    print("Please provide paths for required inputs.")

                    if _prompt_yes_no("Do you want to load a precomputed instance instead (--load/-l)?", default=False):
                        load_path = _prompt_required_path("Precomputed save directory (--load/-l)", must_exist=True, is_dir=True)
                        Data = DataStorage(load_path)
                        _run_shiny(Data)
                        return

                    if not arg_parser.dir:
                        arg_parser.dir = str(_prompt_required_path("M2M data directory (--dir/-d)", must_exist=True, is_dir=True))
                    if not arg_parser.metadata:
                        arg_parser.metadata = str(_prompt_required_path("Metadata TSV file (--metadata/-m)", must_exist=True, is_dir=False))
                    if not arg_parser.output:
                        arg_parser.output = str(_prompt_required_path("Output directory (--output/-o)", must_exist=False, is_dir=True))

                    if not arg_parser.taxonomy:
                        taxonomy_prompt = "Taxonomy TSV file (--taxonomy/-t, optional; Enter to skip)"
                        taxonomy_path = _prompt_optional_path(taxonomy_prompt, expect_file=True)
                        if taxonomy_path is not None:
                            arg_parser.taxonomy = str(taxonomy_path)

                    if not arg_parser.abundance:
                        abundance_prompt = "Abundance TSV file (--abundance/-a, optional; Enter to skip)"
                        abundance_path = _prompt_optional_path(abundance_prompt, expect_file=True)
                        if abundance_path is not None:
                            arg_parser.abundance = str(abundance_path)

                    if not arg_parser.metacyc:
                        metacyc_prompt = "Metacyc PADMet file (--metacyc/-c, optional; Enter to skip)"
                        metacyc_path = _prompt_optional_path(metacyc_prompt, expect_file=True)
                        if metacyc_path is not None:
                            arg_parser.metacyc = str(metacyc_path)
                else:
                    print(MESSAGE)
                    print("\nError: Required arguments missing.")
                    print("Data directory (-d), metadata (-m) and output(-o) are required.")
                    print("Use --help for more information about required arguments.")
                    return

            dir_path = Path(arg_parser.dir).resolve()
            metadata_path = Path(arg_parser.metadata).resolve()

            taxonomic_path = _optional_path(arg_parser.taxonomy, expect_file=True)
            abundance_path = _optional_path(arg_parser.abundance, expect_file=True)

            save_path = Path(arg_parser.output).resolve()

            if not save_path.is_dir():
                try:
                    save_path.mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    print(e)
                    return

            metacyc = _optional_path(arg_parser.metacyc, expect_file=True)

            du.build_dataframes(dir_path, metadata_path, abundance_path, taxonomic_path, save_path, metacyc)

            Data = DataStorage(save_path)

            _run_shiny(Data)
    except KeyboardInterrupt:
        print("\nExiting...")
    except InteractivePromptExit:
        print("Exiting interactive mode.")
