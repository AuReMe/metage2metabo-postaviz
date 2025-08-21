import json
import sys
import tarfile
import tempfile
import time
from multiprocessing import Pool
from multiprocessing import cpu_count
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from padmet.classes import PadmetRef
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci
from scipy import stats
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
from statsmodels.stats.multitest import multipletests

try:
  import polars as pl

except Exception as e:
  print(e)
  print("CALLING subprocess to uninstall polars and polars-lts-cpu. Then reinstall polars-lts-cpu.")

  import subprocess
  import sys

  subprocess.check_call([sys.executable, "-m", "pip", "uninstall", "-y", "polars"])  # noqa: S603
  subprocess.check_call([sys.executable, "-m", "pip", "uninstall", "-y", "polars-lts-cpu"])  # noqa: S603
  subprocess.check_call([sys.executable, "-m", "pip", "install", "polars-lts-cpu"])  # noqa: S603

  import polars as pl

def file_exist(filename: str, directory_path: Path) -> bool:
    for file in directory_path.iterdir():
        if filename == file.name:
            return True
    return False


def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, "__dict__"):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def is_valid_dir(dirpath: Path):
    """Return True if directory exists or not

    Args:
        dirpath (str): path of directory
    Returns:
        bool: True if dir exists, False otherwise
    """
    if dirpath.is_dir:
        return True
    else:
        return False


def extract_tarfile(tar_file: Path, outdir):
    file = tarfile.open(tar_file, "r:gz")
    file.extractall(outdir, filter="data")


def has_only_unique_value(dataframe , input1, input2: str = "None"):
    """
    Return True if the dataframe's column(s) only has unique value, False otherwise.

    Args:
        dataframe (pd.DataFrame): _description_
        column_value (_type_): _description_
        input1 (_type_): _description_
        input2 (str, optional): _description_. Defaults to "None".
    """
    nb_row = len(dataframe)

    if isinstance(dataframe, pl.DataFrame):

        if input2 == "None":
            return True if nb_row == len(dataframe.get_column(input1).unique()) else False

        else:
            return True if nb_row == len(dataframe.get_column(input1).unique()) and nb_row == len(dataframe.get_column(input2).unique()) else False

    if isinstance(dataframe, pd.DataFrame):

        if input2 == "None":
            return True if nb_row == len(dataframe[input1].unique()) else False

        else:
            return True if nb_row == len(dataframe[input1].unique()) and nb_row == len(dataframe[input2].unique()) else False


def relative_abundance(abundance_path: Path, save_path: Path, cscope_dir: Path, scope: str):
    """Generate a second main_dataframe with the production based on weight from the abundance matrix.

    Args:
        abundance_matrix (Path): Pathlib Path of the abundance file.
        sample_cscope (Path): Pathlib Path of the cscope directory.
        save_path (Path): Pathlib Path of the output directory.

    Raises:
        RuntimeError: If more than one column of type other than numeric.

    Returns:
        Dataframe: production dataframe with sample in rows and compounds in column. Weighted by abundance.
    """

    if file_exist("normalised_abundance_dataframe_postaviz.tsv", save_path) and file_exist("abundance_file.tsv", save_path):
        print("Abundance dataframe already exist in save directory.")
        return
    # Read csv with pandas to avoid Polars schema lenght limit. Reset index and transform to polars.
    abundance_df = pd.read_csv(abundance_path, sep="\t").reset_index(names="pd_index")
    abundance_df = pl.DataFrame(abundance_df)

    # Check to see if reset index produced a default Int64 index. If yes remove it.
    if abundance_df.get_column("pd_index").dtype == pl.Int64:
        abundance_df = abundance_df.drop("pd_index")

    # The list of columns dtype String MAYBE USE NON NUMERIC ?
    non_numeric_col = abundance_df.select(pl.col(pl.String)).columns

    # Only one column should be Str dtype. (indentifier column)
    if len(non_numeric_col) == 1:
        abundance_df = abundance_df.with_columns(pl.col(non_numeric_col[0]).alias("smplID")).drop(non_numeric_col[0])
    else:
        raise RuntimeError(f"The numbers of non numeric columns ({len(non_numeric_col)}) isn't equal to one.\n Only one STR (non-numeric at least) column must be present in the abundance dataframe to be used as identifier index.")
    
    # Abundance_file.tsv needed in another function. Save with pandas index.
    abundance_matrix = abundance_df.to_pandas().set_index("smplID")
    abundance_matrix.to_csv(Path(save_path,"abundance_file.tsv"),sep="\t")

    # Loop replacing Pandas apply and to normalise dataframe.
    abundance_df = abundance_df.with_columns((pl.col(col)/pl.col(col).sum()).round(7) for col in abundance_df.select(pl.exclude("smplID")).columns)

    # Abundance_file_normalised needed in another function. Save it with pandas.
    abundance_matrix_normalised = abundance_df.to_pandas().set_index("smplID")
    abundance_matrix_normalised.to_csv(Path(save_path,"abundance_file_normalised.tsv"),sep="\t")

    # List filled with one rows dataframe from each samples.
    res = []

    for sample_path in cscope_dir.iterdir():

        sample_id = sample_path.name.split(".parquet")[0]
        sample_data = pl.read_parquet(sample_path)
        tmp_abundance_df = abundance_df.select(["smplID",sample_id])

        # Transpose allow to work on columns rather than rows which is more efficient.
        sample_df = sample_data.transpose(include_header=True, header_name="Compounds", column_names="smplID")
        sample_df = sample_df.with_columns(pl.col(col).mul(tmp_abundance_df.filter(pl.col("smplID") == col).get_column(sample_id).item()).name.keep() for col in sample_df.select(pl.exclude("Compounds")).columns)
        sample_df = sample_df.transpose(include_header=True, header_name="smplID", column_names="Compounds")
        sample_df = sample_df.sum().with_columns(pl.col("smplID").fill_null(sample_id))
        res.append(sample_df)

    final_df = pl.concat(res, how="diagonal").fill_null(0)
    final_df = final_df.to_pandas()
    final_df.set_index("smplID", inplace=True)
    if scope == "cscope":
        final_df.to_csv(Path(save_path,"normalised_abundance_dataframe_postaviz.tsv"),sep="\t")
    else:
        final_df.to_csv(Path(save_path,"normalised_iscope_abundance_dataframe_postaviz.tsv"),sep="\t")

    print("Abundance done.")


def open_tsv(file_name: str, convert_cpd_id: bool = False, rename_columns: bool = False, first_col: str = "smplID"):
    """Open tsv file as a pandas dataframe.

    Args:
        file_name (str): Path of the file
        rename_columns (bool, optional): Rename the first column and decode the metabolites names in sbml format into readable format. Defaults to False.
        first_col (str, optional): Label of the first col if rename_columns is True. Defaults to "smplID".

    Returns:
        Dataframe: Pandas dataframe
    """
    data = pd.read_csv(file_name, sep="\t")
    if rename_columns:
        data.columns.values[0] = first_col
    if convert_cpd_id:
        data.set_index(first_col,inplace=True,drop=True)
        data.columns = sbml_to_classic(data.columns.values)
    return data


def get_scopes(file_name, path: Path) -> pd.DataFrame:
    for root, _dirs, files in path.walk():
        if file_name in files:
            scope_matrix = open_tsv(Path(root, file_name), convert_cpd_id=True, rename_columns=True)
            return scope_matrix


def is_indexed_by_id(df: pd.DataFrame):
    if df.index.name == "smplID":
        return True
    else:
        return False


def sbml_to_classic(compounds_list):
    uncoded = []
    for coded in compounds_list:
        id, id_type, compart = cfci(coded)
        new_value = str(id)+"["+str(compart)+"]"
        uncoded.append(new_value)
    return uncoded


def retrieve_all_cscope(sample, dir_path: Path, cscope_directoy: Path, cscope_file_format):
    """Retrieve iscope, cscope, added_value and contribution_of_microbes files in the path given using os.listdir().

    Args:
        path (str): Directory path

    Returns:
        dict: Return a nested dict object where each key is a dictionnary of a sample. The key of those second layer dict [iscope, cscope, advalue, contribution] give acces to these files.
    """

    sample_directory_path = Path(dir_path, sample)

    cscope_dataframe = get_scopes("rev_cscope.tsv", sample_directory_path)

    if isinstance(cscope_dataframe,pd.DataFrame):

        bin_list = cscope_dataframe.index.tolist()
        cscope_dataframe.to_parquet(Path(cscope_directoy, sample + cscope_file_format), compression="gzip")
        return bin_list, sample

    else:

        return None, sample


def retrieve_all_iscope(sample, dir_path, iscope_directoy, iscope_file_format):
    """Retrieve iscope, cscope, added_value and contribution_of_microbes files in the path given using os.listdir().

    Args:
        path (str): Directory path

    Returns:
        dict: Return a nested dict object where each key is a dictionnary of a sample. The key of those second layer dict [iscope, cscope, advalue, contribution] give acces to these files.
    """

    sample_directory_path = Path(dir_path, sample)

    if sample_directory_path.is_dir():

        iscope_dataframe = get_scopes("rev_iscope.tsv", sample_directory_path)

        if isinstance(iscope_dataframe,pd.DataFrame):

            iscope_dataframe.to_parquet(Path(iscope_directoy, sample+iscope_file_format), compression="gzip")


def load_sample_cscope_data(dir_path: Path, cscope_directory: Path, cscope_file_format, save_path: Path): # Need rework
    """Open all directories given in -d path input. Get all cscopes, load and save them a dataframe in parquet.gzip format.
    No RAM used during process that way.

    Args:
        path (str): Path of directory

    Returns:
        dict: sample_data dictionnary
    """

    if Path(save_path,"sample_info.json").is_file():
        return

    print("Load/Save cscope data as parquet format...")

    nb_cpu = cpu_count() - 1
    if type(nb_cpu) is not int or nb_cpu < 1:
        nb_cpu = 1
    pool = Pool(nb_cpu)

    pool = Pool(nb_cpu)
    results_list = pool.starmap(retrieve_all_cscope,[(sample.name, dir_path, cscope_directory, cscope_file_format) for sample in dir_path.iterdir()])

    pool.close()
    pool.join()

    sample_info = {}
    sample_info["bins_list"] = []
    sample_info["bins_count"] = {}
    sample_info["bins_sample_list"] = {}

    for all_bins_in_sample, sample in results_list:

        if type(all_bins_in_sample) != list:  # noqa: E721
            continue

        sample_info["bins_list"] = sample_info["bins_list"] + all_bins_in_sample

        for bin in all_bins_in_sample:

            if bin not in sample_info["bins_count"]:
                sample_info["bins_count"][bin] = 0

            sample_info["bins_count"][bin] += 1

            if bin not in sample_info["bins_sample_list"]:
                sample_info["bins_sample_list"][bin] = []

            sample_info["bins_sample_list"][bin].append(str(sample))

    # Remove duplicate from list
    sample_info["bins_list"] = list(dict.fromkeys(sample_info["bins_list"]))

    with open(Path(save_path,"sample_info.json"), "w") as f:
            json.dump(sample_info, f)

    print("Done.")

    return sample_info


def load_sample_iscope_data(dir_path: Path, iscope_directory: Path, iscope_file_format):
    """Load and save iscope data as dataframe in parquet.gzip format.

    Args:
        dir_path (Path): Directory path given in cli (-d)
        iscope_directory (Path): Path of newly created save directory.
        iscope_file_format (bool): Format to save.
    """

    print("Load/Save iscope data as parquet format...")

    nb_cpu = cpu_count() - 1
    if type(nb_cpu) is not int or nb_cpu < 1:
        nb_cpu = 1
    pool = Pool(nb_cpu)

    pool = Pool(nb_cpu)
    pool.starmap(retrieve_all_iscope,[(sample.name, dir_path, iscope_directory, iscope_file_format) for sample in dir_path.iterdir()])

    pool.close()
    pool.join()


def build_main_dataframe(save_path: Path, cscope_directory: Path):
    """Create and save the main dataframe. Samples in rows and compounds in columns.
    It takes the compounds production in each samples cscope and return a pandas Series with 1 produced or 0 absent for each compounds.
    Merge all the series returned into a dataframe.

    Args:
        sample_data (dict): Samples's cscope.
        save_path (_type_): Save path given in CLI.
    """
    if file_exist("main_dataframe_postaviz.tsv", save_path):
        print("main dataframe already in save directory.")
        return

    print("Building main dataframe...")

    all_series = []

    for sample_filename in cscope_directory.iterdir():

        sample_id = sample_filename.name.split(".parquet")[0]

        current_sample_df = pd.read_parquet(Path(cscope_directory, sample_filename))

        # Get all the compounds produced in this cscope.
        serie_index = current_sample_df.columns.values

        serie_data = []

        # Since at least one bin must produce all of those compounds it is considered produced in the entire sample.
        for _i in range(len(serie_index)):

            serie_data.append(1)

        # Make the Pandas Series with compounds as index and 1 as value and append in global for concatenation.
        all_series.append(pd.Series(data=serie_data,index=serie_index,name=sample_id))

    # Concatenation of all the serie into one dataframe.
    results = pd.concat(all_series, axis=1).T
    results.fillna(0,inplace=True)
    results = results.astype(int)
    results.index.name = "smplID"

    results.to_csv(Path(save_path, "main_dataframe_postaviz.tsv"),sep="\t")

    print("Main dataframe done and saved.")


def build_dataframes(dir_path: Path, metadata_path: Path, abundance_path: Optional[Path] = None, taxonomic_path: Optional[Path] = None, save_path: Optional[Path] = None, metacyc: Optional[Path] = None):
    """Main function.
    dir_path, metadata_path and save_path are necessary.
    Generate most of the core dataframes to avoid calculation on the application side.

    Args:
        dir_path (Path): Directory path containing M2M output.
        metadata_path (Path): Metadata file path.
        abundance_path (Optional[Path], optional): Abundance file path. Defaults to None.
        taxonomic_path (Optional[Path], optional): Taxonomic file path. Defaults to None.
        save_path (Optional[Path], optional): Output path. Defaults to None.
        metacyc (Optional[Path], optional): Metacyc DB file path. Defaults to None.
    """
    if not is_valid_dir(dir_path):
        print(dir_path, "Sample directory path is not a valid directory")
        sys.exit(1)

    if save_path is None:
        sys.exit(1)

    if not is_valid_dir(save_path):
        save_path.mkdir(parents=True,exist_ok=True)

    metadata_processing(metadata_path, save_path)

    cscope_directory = Path(save_path,"sample_cscope_directory")

    if not cscope_directory.is_dir():

        cscope_directory.mkdir()

    cscope_file_format = ".parquet.gzip"

    load_sample_cscope_data(dir_path, cscope_directory, cscope_file_format, save_path)

    producers_dataframe(cscope_directory, save_path, "cscope")

    build_main_dataframe(save_path, cscope_directory)

    build_compounds_index(save_path)

    if abundance_path is not None:

        relative_abundance(abundance_path, save_path, cscope_directory, "cscope")

    if taxonomic_path is not None:

        taxonomy_processing(taxonomic_path, save_path)

    total_production_by_sample(save_path, abundance_path)

    build_pcoa_dataframe(save_path)

    iscope_directory = Path(save_path,"sample_iscope_directory")

    if not iscope_directory.is_dir():

        iscope_directory.mkdir()

        iscope_file_format = ".parquet.gzip"

        load_sample_iscope_data(dir_path, iscope_directory, iscope_file_format)

    relative_abundance(abundance_path, save_path, iscope_directory, "iscope")

    producers_dataframe(iscope_directory, save_path, "iscope")

    iscope_cscope_fill_difference(save_path)

    bin_dataframe_build(cscope_directory, "cscope", abundance_path, taxonomic_path, save_path)

    concat_bin_dataframe(save_path, "cscope")

    bin_dataframe_build(iscope_directory, "iscope", abundance_path, taxonomic_path, save_path)

    concat_bin_dataframe(save_path, "iscope")

    # Metacyc database TREE

    if metacyc is not None:

        padmet_to_tree(save_path, metacyc)


def metadata_processing(metadata_path: Path, save_path: Path):
    """Simple function to save the metadata as parquet file.
    allow for dtypes change of the file in application while creating a safe copy of the original file.

    Args:
        metadata_path (Path): Path to metadata file.
        save_path (Path): Saving path.

    Returns:
        None if file already exist in save_path.
    """
    if file_exist("metadata_dataframe_postaviz.parquet.gzip", save_path):
        return
    else:
        metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
        metadata = metadata.rename(columns={metadata.columns[0]: "smplID"})
        metadata.to_parquet(Path(save_path,"metadata_dataframe_postaviz.parquet.gzip"), index= True if is_indexed_by_id(metadata) else False)


def total_production_by_sample(save_path: Path, abundance_path: Optional[Path] = None):
    """Create and save the total production dataframe. This dataframe contain all samples in row and all compounds in columns.
    For each samples the compounds produced by each bins is sum up to get the estimated total production of compound by samples
    and the number of bins who produced those compounds.

    If the abundance is provided, each production (1) of bins is multiplied by their abundance in their sample which gives an
    estimated production of compounds weighted by the abundance of the bin producer.

    Args:
        save_path (_type_): Save path given in CLI
        abundance_path (str, optional): Abundance file path fiven in CLI. Defaults to None.
    """

    if file_exist("total_production_dataframe_postaviz.tsv", save_path):
        print("total_production_dataframe_postaviz already in save directory")
        return

    print("Building total production dataframe...")

    main_dataframe = open_tsv(Path(save_path, "main_dataframe_postaviz.tsv"))
    metadata_dataframe = pd.read_parquet(Path(save_path, "metadata_dataframe_postaviz.parquet.gzip"))

    if not is_indexed_by_id(main_dataframe):
        main_dataframe.set_index("smplID",inplace=True,drop=True)

    main_dataframe["Total_production"] = main_dataframe.apply(lambda row: row.to_numpy().sum(), axis=1)
    results = pd.DataFrame(main_dataframe["Total_production"])

    if abundance_path is not None:
        abundance_dataframe = open_tsv(Path(save_path, "normalised_abundance_dataframe_postaviz.tsv"))

        if not is_indexed_by_id(abundance_dataframe):
            abundance_dataframe.set_index("smplID",inplace=True,drop=True)

        abundance_dataframe["Total_abundance_weighted"] = abundance_dataframe.apply(lambda row: row.to_numpy().sum(), axis=1)
        abundance_dataframe = abundance_dataframe["Total_abundance_weighted"]

        results = pd.concat([results,abundance_dataframe], axis=1)

    results.reset_index(inplace=True)
    results = results.merge(metadata_dataframe,"inner","smplID")

    results.to_csv(Path(save_path, "total_production_dataframe_postaviz.tsv"),sep="\t", index=False)

    print("Total production dataframe done and saved.")

    return


def preprocessing_for_statistical_tests(dataframe: pd.DataFrame, y_value, input1, input2 = None, multipletests: bool = False, multipletests_method: str = "bonferroni"):
    """Create dataframe for each y_value in the list, to separate them and use wilcoxon_man_whitney function.
    Concat all results into one dataframe.

    Args:
        dataframe (pd.DataFrame): Dataframe to test.
        y_value (_type_): list of columns labels to separate into several dataframe. Must be at least of lenght 1.
        input1 (_type_): First user's input.
        input2 (_type_, optional): Second user's input. Defaults to None.

    Returns:
        Dataframe: Dataframe of statistical test.
    """
    all_results = []

    for y in y_value:

        if input2 is None:
            all_results.append(wilcoxon_man_whitney(dataframe[[y, input1]], y, input1, None, multipletests, multipletests_method))
        else:
            all_results.append(wilcoxon_man_whitney(dataframe[[y, input1, input2]], y , input1, input2, multipletests, multipletests_method))

    return pd.concat(all_results)


def wilcoxon_man_whitney(dataframe: pd.DataFrame, y, first_factor: str, second_factor: Optional[str] = None, multiple_correction: bool = False, correction_method: str = "hs"):
    """
    Takes one dataframe with only one value column y and return a dataframe of statistical tests.
    First all sub arrays by the first input then the second input are made and convert to numpy array.
    Then Wilcoxon or Mann Whitney test are run on each pair without doublon.
    If pairs array have the same lenght -> Wilcoxon, if not -> Mann Whitney

    Args:
    dataframe (pd.Dataframe): Pandas dataframe
    y (str): Column label containing the values to test.
    first_factor (str): Column label of the first user's input.
    second_factor (str): Column label of the second user's input. Default to None


    Returns:
        Dataframe: Dataframe of test's results.
    """
    # Array sorting by the first input and the second input if NOT None.
    sub_dataframes = {}

    for first_factor_array in dataframe[first_factor].unique():

        if second_factor is None:

            sub_dataframes[first_factor_array] = dataframe.loc[dataframe[first_factor] == first_factor_array][y].to_numpy()
            continue

        sub_dataframes[first_factor_array] = {}

        for second_factor_array in dataframe[second_factor].unique():

            sub_dataframes[first_factor_array][second_factor_array] = dataframe.loc[(dataframe[first_factor] == first_factor_array) & (dataframe[second_factor] == second_factor_array)][y].to_numpy()

    # Dataframe's structure declaration, Axis column added if second input is NOT None.
    if second_factor is None:

        results = pd.DataFrame(columns=["Compound", "Factor1", "Sample size1", "Factor2", "Sample size2", "Method", "Statistic", "Pvalue", "Significance"])

    else:

        results = pd.DataFrame(columns=["Compound", "Axis", "Factor1", "Sample size1", "Factor2", "Sample size2", "Method", "Statistic", "Pvalue", "Significance"])

    # Test each pairs avoiding useless duplicates and Array of lenght <= 1.
    for name in sub_dataframes.keys():  # noqa: PLC0206

        if second_factor is None: # One input selected

            for name2 in sub_dataframes.keys():  # noqa: PLC0206

                if name == name2:
                    continue

                if len(sub_dataframes[name]) < 1 or len(sub_dataframes[name2]) < 1:
                    continue

                if name2 in results["Factor1"].tolist():
                    continue

                if len(sub_dataframes[name]) == len(sub_dataframes[name2]):

                    test_value, pvalue = stats.wilcoxon(sub_dataframes[name], sub_dataframes[name2])
                    test_method = "Wilcoxon"

                else:

                    test_value, pvalue = stats.mannwhitneyu(sub_dataframes[name], sub_dataframes[name2])
                    test_method = "Mann-Whitney"

                if pvalue >= 0.05:
                    symbol = "ns"
                elif pvalue >= 0.01:
                    symbol = "*"
                elif pvalue >= 0.001:
                    symbol = "**"
                else:
                    symbol = "***"

                results.loc[len(results)] = {"Compound": y,
                                            "Factor1": name, "Sample size1": len(sub_dataframes[name]),
                                              "Factor2": name2, "Sample size2": len(sub_dataframes[name2]),
                                                "Statistic": test_value, "Pvalue": pvalue, "Significance": symbol, "Method": test_method}

        else: # Two inputs selected

            for name2 in sub_dataframes[name].keys():

                for name3 in sub_dataframes[name].keys():

                    if name2 == name3:
                        continue

                    if len(sub_dataframes[name][name2]) < 1 or len(sub_dataframes[name][name3]) < 1:
                        continue

                    if len(results.loc[(results["Factor1"] == str(second_factor+": "+str(name3))) & (results["Axis"] == name)]) > 0: # Avoid duplicate
                        continue

                    if len(sub_dataframes[name][name2]) == len(sub_dataframes[name][name3]):

                        test_value, pvalue = stats.wilcoxon(sub_dataframes[name][name2], sub_dataframes[name][name3])
                        test_method = "Wilcoxon"

                    else:

                        test_value, pvalue = stats.mannwhitneyu(sub_dataframes[name][name2], sub_dataframes[name][name3])
                        test_method = "Mann-Whitney"

                    if pvalue >= 0.05:
                        symbol = "ns"
                    elif pvalue >= 0.01:
                        symbol = "*"
                    elif pvalue >= 0.001:
                        symbol = "**"
                    else:
                        symbol = "***"

                    results.loc[len(results)] = {"Compound": y, "Axis": name,
                                                "Factor1": str(second_factor+": "+str(name2)), "Sample size1": len(sub_dataframes[name][name2]),
                                                "Factor2": str(second_factor+": "+str(name3)), "Sample size2": len(sub_dataframes[name][name3]),
                                                    "Statistic": test_value, "Pvalue": pvalue, "Significance": symbol, "Method": test_method}

    if multiple_correction:

        pvals_before_correction = results["Pvalue"].to_numpy()
        reject, pvals_after_correction, _, __ = multipletests(pvals_before_correction, method = correction_method)

        results["Pvalue corrected"] = pvals_after_correction
        results["Significance corrected"] = results["Pvalue corrected"].apply(lambda x:get_significance_symbol(x))
        results["Correction method"] = correction_method

    return results


def get_significance_symbol(pval: float) -> str:
    """Return Significance symbol depending on pvalue given.

    Args:
        pval (float): Pvalue of the test

    Returns:
        str: Significance's symbol
    """
    if pval >= 0.05:
        return "ns"
    elif pval >= 0.01:
        return "*"
    elif pval >= 0.001:
        return "**"
    else:
        return "***"


def build_pcoa_dataframe(save_path: Path) -> pd.DataFrame:
    """Compute Principal Coordinate Analysis from the main_dataframe given in input. Merge with metadata from the smplID column or index.

    Args:
        main_dataframe (pd.DataFrame): Dataframe from which the pcoa will be made.
        metadata (pd.DataFrame): Metadata dataframe (Must have smplID identifer column for the merge to work)

    Returns:
        pd.DataFrame: Pcoa dataframe with sample ID as index, PC1 and PC2 results and all metadata.
    """

    if file_exist("pcoa_dataframe_postaviz.tsv", save_path):
        print("Pcoa dataframe already in save directory.")
        return

    main_dataframe = open_tsv(Path(save_path, "main_dataframe_postaviz.tsv"))

    metadata = pd.read_parquet(Path(save_path, "metadata_dataframe_postaviz.parquet.gzip"))

    if not is_indexed_by_id(main_dataframe):
        main_dataframe = main_dataframe.set_index("smplID")

    if is_indexed_by_id(metadata):
        metadata = metadata.reset_index("smplID")

    dmatrix = main_dataframe.to_numpy()
    dist_m = pdist(dmatrix, "jaccard")
    square_df = squareform(dist_m)
    pcoa_result = pcoa(square_df,number_of_dimensions=2)
    coordinate = pcoa_result.samples

    df_pcoa = coordinate[["PC1","PC2"]]
    df_pcoa["smplID"] = main_dataframe.index.to_numpy()

    df_pcoa = df_pcoa.merge(metadata, "inner", "smplID")
    df_pcoa.set_index("smplID",inplace=True)

    df_pcoa.to_csv(Path(save_path,"pcoa_dataframe_postaviz.tsv"),sep="\t")

    return


def correlation_test(value_array, factor_array, factor_name, method:str = "pearson"):

    if method == "pearson":
        res = stats.pearsonr(value_array, factor_array)
    else:
        res = stats.spearmanr(value_array, factor_array)

    if res.pvalue >= 0.05:
        symbol = "ns"
    elif res.pvalue >= 0.01:
        symbol = "*"
    elif res.pvalue >= 0.001:
        symbol = "**"
    else:
        symbol = "***"

    return pd.DataFrame([[factor_name, len(value_array), res.statistic, res.pvalue, symbol, method]],columns=["Factor", "Sample size", "Statistic", "Pvalue", "Significance", "Method"])


def taxonomy_processing(taxonomy_filepath: Path, save_path: Path):
    """Open and save taxonomy file.

    Args:
        taxonomy_filepath (str): TSV or TXT format

    Raises:
        RuntimeError: Wrong file's format

    Returns:
        pd.DataFrame: Pandas dataframe
    """

    if file_exist("taxonomic_dataframe_postaviz.tsv", save_path):
        print("Taxonomic dataframe already exist in save directory.")
        return

    print("Building taxonomic dataframe...")

    if taxonomy_filepath.suffix == ".tsv":

        df = open_tsv(taxonomy_filepath)
        df = df.rename(columns={df.columns[0]: "mgs"})
        df.to_csv(Path(save_path,"taxonomic_dataframe_postaviz.tsv"), sep="\t", index=False)

        return

    if not taxonomy_filepath.suffix == ".txt":
        raise RuntimeError("Taxonomy file must be either a txt file or tsv file.")

    with open(taxonomy_filepath) as f:
        lines = f.readlines()

    df = pd.DataFrame(columns=["mgs","kingdom","phylum","class","order","family","genus"])

    del lines[0] # Delete header line

    for row in lines:

        mgs = row.split("\t")[0]
        genus = row.split("\t")[1:]

        k, p, c, o, f, g = genus[0].strip("\n").split(";")

        df.loc[len(df)] = [mgs,k,p,c,o,f,g]


    df.to_csv(Path(save_path,"taxonomic_dataframe_postaviz.tsv"), sep="\t", index=False)

    print("Taxonomic dataframe done and saved.")


def bin_dataframe_build(scope_directory: Path, scope_mode: str = "cscope", abundance_path: Optional[Path] = None, taxonomy_path: Optional[Path] = None, savepath: Optional[Path] = None):
    """Build a large dataframe with all the bins of the different samples as index, the dataframe contain the list of production, the abundance fot he bin in the sample
    and the count of production with or without abundance.

    Args:
        sample_info (dict): _description_
        sample_data (dict): _description_
        metadata (Dataframe): _description_
        abundance_file (Dataframe, optional): _description_. Defaults to None.
        taxonomy_path (Dataframe, optional): _description_. Defaults to None.

    Returns:
        pd.DataFrame: Pandas dataframe
    """

    if file_exist("bin_dataframe_chunk_1.parquet.gzip",savepath):
        print("Chunk of bin_dataframe already in save directory")
        return

    print("Building bin dataframe...")

    start = time.time()

    cpd_index = pd.read_csv(Path(savepath,"compounds_index.tsv"), sep="\t").squeeze()

    ##### Abundance normalisation, give percentage of abundance of bins in samples.
    if abundance_path is not None:

        abundance_file_normalised = pd.read_csv(Path(savepath,"abundance_file_normalised.tsv"),sep="\t",index_col=0)

    if taxonomy_path is not None:

        taxonomy_file = open_tsv(Path(savepath,"taxonomic_dataframe_postaviz.tsv"))

        ##### Checks if taxonomy file has default index which mean it is not indexed. TEMPORARY until found better way to deal with open/save from -t option OR load taxonomic_df option which return non indexed / indexed df

        if not pd.Index(np.arange(0, len(taxonomy_file))).equals(taxonomy_file.index):
                taxonomy_file = taxonomy_file.reset_index()

        mgs_col_taxonomy = taxonomy_file.columns[0]

    sample_unique_list = list(scope_directory.iterdir())

    ##### Create Generator to process data by chunks.
    print("Making chunk of sample list...")

    chunk_generator = chunks(sample_unique_list, 250)

    ##### Loop throught generator

    chunk_index = 0

    for current_chunk in chunk_generator:

        chunk_index += 1
        list_of_dataframe = []
        logs = []

        for sample in current_chunk:

            sample_id = sample.name.split(".parquet")[0]

            try:

                df = pd.read_parquet(Path(scope_directory, sample))
                df.insert(0 , "smplID", sample_id)
                df.index.name = "binID"
                list_of_dataframe.append(df)

            except Exception as e:

                logs.append(f"No dataframe named {sample} in scope_directory.\n{e}")

        results = pd.concat(list_of_dataframe)
        results.fillna(0,inplace=True)

        tmp = results.apply(lambda x: x.index[x == 1].tolist(), axis=1)

        results["Production"] = tmp.apply(lambda value: get_cpd_index(cpd_index, value))

        results["Count"] = results.apply(lambda x: len(x.index[x == 1].tolist()), axis=1)
        results = results[["smplID","Production","Count"]]

        del tmp

        if abundance_path is not None: # If abundance is provided, multiply each Count column with the relative abundance of the bins in their samples.

            abundance = results.apply(lambda row: abundance_file_normalised.at[row.name,row["smplID"]],axis=1)
            abundance.name = "Abundance"

            final_result = pd.concat([results,abundance], axis=1)
            del results

            final_result["Count_with_abundance"] = final_result["Count"] * final_result["Abundance"]

        else: # Seems quite useless but felt cute might deleted later.

            final_result = results
            del results

        final_result = final_result.reset_index() # .merge(metadata, "inner", "smplID")

        if taxonomy_path is not None: # If taxonomy is provided, merge the dataframe with the taxonomic_dataframe.

            final_result = final_result.merge(taxonomy_file, "inner", left_on="binID", right_on=mgs_col_taxonomy)

        ##### Save current chunks into parquet file

        filename = "bin_dataframe_chunk_"+str(chunk_index)+"-"+str(scope_mode)+".parquet.gzip"
        filepath = Path(savepath,filename)

        if len(final_result) == 0:
            print(f"Chunks {chunk_index} is empty !")

        final_result.to_parquet(filepath, compression="gzip")

        del final_result

    print("Took: ", time.time() - start, "Before saving")

    if len(logs) != 0:
        with open(Path(savepath,"bin_dataframe_logs.txt"),"w") as log_file:
            for line in logs:
                log_file.write(f"{line}\n")

    return


def chunks(lst, n):
    """Yield successive n-sized chunks from list."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def padmet_to_tree(save_path: Path, metacyc_file_path: Path):
    """Build a tree to be used in the Shiny application.
    Allow the user to select directly a compounds or a category of compounds and fill a list
    with all the compounds corresponding to that category.

    Use the function build_parent_child_dataframe to create a 2 columns (child_id/parent_id) dataframe.
    With the relation dataframe, build the tree using build_tree_from_root.

    Args:
        save_path (str): Path of the save directory.
    """

    if file_exist("padmet_compounds_category_tree.json",save_path):
        print("Padmet category tree already exist.")
        return

    print("Building compounds category tree...")

    padmet = PadmetRef(str(metacyc_file_path))

    cpd_id = [node.id for node in padmet.dicOfNode.values() if node.type == "compound"]

    # classes_id = [node.id for node in padmet.dicOfNode.values() if node.type == "class"]

    df = pd.DataFrame(columns=["child_id", "parent_id"])

    for mid in cpd_id:

        build_parent_child_dataframe(padmet, df, mid)

        if mid not in df["child_id"].values:

            print("ERROR for : ",mid)

    root = {}

    root["FRAMES"] = {}

    build_tree_from_root(root["FRAMES"], "FRAMES", df)
    root["All_metabolites"] = root.pop("FRAMES")

    with open(Path(save_path, "padmet_compounds_category_tree.json"), "w") as fp:
        json.dump(root, fp)

    print("Compounds category tree done.")


def build_parent_child_dataframe(padmet: PadmetRef, dataframe: pd.DataFrame, current_id, child_column = "child_id", parent_column = "parent_id"):
    """Build a child /parent relation dataframe between compounds category from a metacyc database in a padmet file format.

    Args:
        padmet (PadmetRef): PadmetRef object from padmet package.
        dataframe (pd.DataFrame): Transmission of the dataframe between recursive call.
        current_id (_type_): ID of the current category / compound.
        child_column (str, optional): Column label of the child column of the dataframe. Defaults to "child_id".
        parent_column (str, optional): Column label of the parent column of the dataframe. Defaults to "parent_id".
    """
    # current_id can be a Compound id or a class id children of the cpd_id / class_id.
    # Check of current_id is already in child column. STOP condition.
    if current_id in dataframe[child_column].values:

        return

    # Get list of relations of current_id
    try:
        rlt_classes = [rlt.id_out for rlt in padmet.dicOfRelationIn[current_id] if rlt.type == "is_a_class"]
    except KeyError as e:
        print(e)
        return

    if len(rlt_classes) == 0:

        return

    # Loop in classes id children of cpd_id OR one of its classes id.
    for rlt_c in rlt_classes:

        # If current_id is NOT in child columns. Write current_id in child column with its FIRST child as parent in parent column.
        if current_id not in dataframe[child_column].values:

            dataframe.loc[len(dataframe)] = {child_column : current_id, parent_column : rlt_c}

        # If the current child is already in child column pass to the next.
        if rlt_c in dataframe[child_column].values:

            continue

        # Then continue with the first child of current_id who is not in child column. Until no more child then pass to next CPD_ID.
        build_parent_child_dataframe(padmet, dataframe = dataframe, current_id = rlt_c)


def build_tree_from_root(node, id, df):
    """Build a tree from a dataframe and a dictionary with the first key as root.
    The first key is the first parent node from which the tree will be built starting with its first child.
    Any node that is not connected indirectly with the root node won't be in the tree.

    Example :   root = {}

                root["FRAMES"] = {}

                build_tree_from_root(root["FRAMES"], "FRAMES", dataframe)

    Args:
        node (dict): root node
        id (str): Root node key id, correspond to the string of the first node in the dataframe.
        df (pd.DataFrame): Dataframe with 2 columns: columns names must be child_id and parent_id. child_id column has only unique values.
    """
    children = df.loc[df["parent_id"] == id]["child_id"].tolist()

    if len(children) == 0:
        return

    for child in children:

        node[child] = {}

        build_tree_from_root(node[child], child, df)


def build_compounds_index(save_path: Path):

    main_dataframe = pd.read_csv(Path(save_path, "main_dataframe_postaviz.tsv"), sep="\t")

    cpd_list = main_dataframe.columns.tolist()
    if "smplID" in cpd_list:
        cpd_list.remove("smplID")
    cpd_index = pd.Series(cpd_list, range(len(cpd_list)))
    saving_path = Path(save_path, "compounds_index.tsv")

    cpd_index.to_csv(saving_path, sep="\t", index=False)


def get_cpd_label(cpd_index: pd.Series, cpd_list_index):

    return cpd_index[cpd_list_index]


def get_cpd_index(cpd_index: pd.Series, cpd_list_label):

    return cpd_index.index[cpd_index.isin(cpd_list_label)].tolist()


def concat_chunk(chunk_dir: Path, save_path: Path, scope_type: str):
    """Concatenation of all sub_dataframes produced.

    Args:
        chunk_dir (Path): Directory path where the chunk are.
        save_path (Path): Save result path
        scope_type (str): Cscope or Iscope
    """
    filename = "producers_"+scope_type+"_dataframe.parquet.gzip"
    df = None
    for chunk in chunk_dir.iterdir():
        if df is None:
            df = pl.read_parquet(chunk)
        else:
            tmp_df = pl.read_parquet(chunk)
            df = pl.concat([df,tmp_df], how="diagonal")

    df = df.fill_null(0)
    df.write_parquet(Path(save_path,filename), compression="gzip")


def sum_and_concat_by_chunk(directory_path: Path):
    """Produce dataframe from chunk of 250 samples, BETTER memory usage small performance price.

    Args:
        directory_path (Path): _description_
    """
    tmp_dir = tempfile.TemporaryDirectory(directory_path.name)

    tmp_dir_path = Path(tmp_dir.name)

    sample_unique_list = list(directory_path.iterdir())

    chunk_generator = chunks(sample_unique_list, 250)
    chunk_index = 0

    for current_chunk in chunk_generator:

        chunk_index += 1
        res = []

        for sample in current_chunk:
            sample_id = sample.name.split(".parquet")[0]
            sample_dataframe = pl.read_parquet(sample)
            sample_dataframe = sample_dataframe.sum()
            sample_dataframe._replace("smplID", pl.Series("smplID",[sample_id])) # replace_column(-1, pl.Series("smplID",[sample_id])) Could also replace smplID col by index since polar put index from parquet to the last column.
            res.append(sample_dataframe)

        final_dataframe = pl.concat(res, how="diagonal")
        final_dataframe = final_dataframe.drop("smplID").insert_column(0, final_dataframe.get_column("smplID")) # Not necessary at all but i'm used to smplID col index with pandas.
        final_dataframe = final_dataframe.fill_null(0)

        filename = "producers_dataframe_chunk_"+str(chunk_index)+".parquet.gzip"
        final_dataframe.write_parquet(Path(tmp_dir_path,filename), compression = "gzip")

    return tmp_dir_path, tmp_dir


def producers_dataframe(scope_directory: Path, save_path: Path, scope_type: str):

    tmp_chunks_dir_path, tmp_dir = sum_and_concat_by_chunk(scope_directory)

    concat_chunk(tmp_chunks_dir_path, save_path, scope_type)

    tmp_dir.cleanup()


def iscope_cscope_fill_difference(save_path: Path):

    iscope_df = pl.read_parquet(Path(save_path, "producers_iscope_dataframe.parquet.gzip"))

    cscope_df = pl.read_parquet(Path(save_path, "producers_cscope_dataframe.parquet.gzip"))

    ccol = cscope_df.columns
    icol = iscope_df.columns

    columns_difference = list(set(ccol) - set(icol))

    for col in columns_difference:

        iscope_df = iscope_df.with_columns(pl.lit(0).alias(col))

    iscope_df.write_parquet(Path(save_path, "producers_iscope_dataframe.parquet.gzip"), compression="gzip")


def concat_bin_dataframe(save_path: Path, scope_mode: str = "cscope"):

    df = None

    for file in save_path.iterdir():
        if file.stem.startswith("bin_dataframe_chunk") and file.stem.split("-")[1].split(".")[0] == scope_mode: # Check if file is a bin_dataframe_chunk AND check the scope of this file.

            if df is None:
                df = pl.read_parquet(file)
                file.unlink()
            else:
                tmp_df = pl.read_parquet(file)
                df = pl.concat([df, tmp_df], how="diagonal")
                file.unlink()

    print("Size of the full bin_dataframe: ",df.estimated_size("gb")," Gb")
    df.write_parquet(Path(save_path, f"bin_dataframe_{scope_mode}.parquet.gzip"), compression="gzip")
