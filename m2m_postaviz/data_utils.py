import json
import os
import os.path
import tarfile
import time
import sys

import numpy as np
import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci
from scipy import stats
from multiprocessing import Pool
from multiprocessing import cpu_count
from os import getpid
from functools import partial

def is_valid_dir(dirpath):
    """Return True if directory exists or not
    
    Args:
        dirpath (str): path of directory
    Returns:
        bool: True if dir exists, False otherwise
    """
    if os.path.isdir(dirpath) == True:
        return True
    else:
        return False


def is_valid_file(filepath):
    """Return True if filepath exists
    Args:
        filepath (str): path to file
    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        open(filepath, 'r').close()
        return True
    except OSError:
        return False


def get_bin_list(sample_data, mode: str = "cscope"):
    bin_list = {}
    for sample in sample_data.keys():
        if is_indexed_by_id(sample_data[sample][mode]):
            bin_list[sample] = sample_data[sample][mode].index.to_list()
        else:
            bin_list[sample] = sample_data[sample][mode]["smplID"].to_list()
    return bin_list


def extract_tarfile(tar_file, outdir):
    file = tarfile.open(tar_file, "r:gz")

    file.extractall(outdir, filter="data")
    # if sys.version_info >= (3, 12):
    # else:
    #     tar.extractall(outdir)


def benchmark_decorator(func):
    def wrapper(*args, **kwargs):
        results = list()
        n_repeats = 3
        for i in range(n_repeats):
            time_start = time.perf_counter()
            result = func(*args, **kwargs)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            results.append(time_duration)
            print(f">run {i+1} took {time_duration} seconds")
        avg_duration = sum(results) / n_repeats
        print(f"Took {avg_duration} seconds on average for {func.__name__} function.")
        return result

    return wrapper


def relative_abundance_calc(abundance_matrix: pd.DataFrame, sample_data: dict):
    """Use the abundance matrix given in input with the dataframe of the sample's production in community to create 2 abundance dataframe.
    1 with normalised bin quantity with x / x.sum(), the other without any transformation.

    Args:
        abundance_file_path (str): path to the abundance file.
        sample_data (dict): dictionnary of the cscope, iscope of each sample.

    Returns:
        Tuple: (Normalised abundance dataframe, Abundance dataframe)
    """
    smpl_norm_abundance = []
    smpl_norm_index = []

    str_filter = abundance_matrix.select_dtypes(include=["string","object","category"])

    if len(str_filter.columns) == 1:
        index_column = str_filter.columns.values[0]
        print(index_column, "column is str type, using it as index")
        abundance_matrix.set_index(index_column,drop=True,inplace=True)

    abundance_matrix_normalised = abundance_matrix.apply(lambda x: x / x.sum(), axis=0)

    for sample in sample_data.keys():

        sample_matrix = sample_data[sample]["cscope"].copy()

        if not is_indexed_by_id(sample_matrix):
            sample_matrix.set_index("smplID", inplace=True)
        
        sample_matrix = sample_matrix.apply(lambda row: row * abundance_matrix_normalised.at[row.name, sample], axis=1)
        sample_matrix = sample_matrix.apply(lambda col: col.to_numpy().sum(), axis=0)
        sample_matrix.name = sample

        smpl_norm_abundance.append(sample_matrix)
        smpl_norm_index.append(str(sample))

    normalized_abundance = pd.concat(smpl_norm_abundance,axis=1)
    normalized_abundance.fillna(0, inplace=True)
    normalized_abundance = normalized_abundance.T
    normalized_abundance.index.name = "smplID"

    return normalized_abundance


def sum_squash_table(abundance_table: pd.DataFrame, sample_id: str):
    """
    Return a dataframe with a unique row containing the sum of all metabolite produced
    by the different bin in a sample. In other word, transform a cpd production df by bin to a cpd production df by sample.
    Args:
        abundance_table (pd.DataFrame): Compound/Bin abundance table
        sample_id (str): Name of the sample

    Returns:
        Dataframe: Dataframe
    """
    # Prend la nouvelle matrice d'abondance du sample
    results = abundance_table.apply(lambda col: col.sum(), axis=0)
    results.name = sample_id
    return results


def get_metadata(sample_id: str, metadata: pd.DataFrame):
    return tuple(
        [metadata.loc[metadata["smplID"] == sample_id].values.tolist()[0], metadata.loc[metadata["smplID"] == sample_id].columns.to_list()]
    )


def taxonomy_groupby(
    metadata: pd.DataFrame,
    current_sample: str,
    bin_id_by_sample: dict,
    taxonomic_dataframe: pd.DataFrame,
    target_rank: str = "s",
    taxonomic_choice: list = [],
):
    """Generate a taxonomic count dataframe from a sample id, his dataframe and the rank choosen. Return only the selected taxonomic choice.

    Args:
        current_sample (str): Sample's id
        taxonomic_dataframe (pd.DataFrame): The taxonomic dataframe
        target_rank (str, optional): Selected rank in shiny's input. Defaults to "Genus".
        taxonomic_choice (list, optional): The list of taxonomic selection from shiny's input to keep in the returned Dataframe. Defaults to [].

    Returns:
        Dataframe: Pandas dataframe containing the selected taxonomic choice and rank.
    """
    df = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(bin_id_by_sample[current_sample])]
    for choice in taxonomic_choice:
        if choice not in df[target_rank].unique():
            taxonomic_choice.remove(choice)
            print(choice, " Removed.")

    df = df[["mgs", target_rank]]
    df = df.groupby([target_rank]).count()
    df = df.reset_index()
    df.columns = [target_rank, "Count"]
    df = df.loc[df[target_rank].isin(taxonomic_choice)]
    value, label = get_metadata(current_sample, metadata)
    for i in range(len(label)):
        df.insert(0, label[i], value[i])
    return df


def taxonomic_dataframe_from_input(
    taxonomic_rank: str, bin_id_by_sample: dict, taxonomic_choice: list, taxonomic_dataframe: pd.DataFrame, metadata: pd.DataFrame
):
    results = []
    taxonomic_choice = list(taxonomic_choice)
    if len(taxonomic_choice) == 0:
        print("The taxonomic choice list is empty")
        return
    for sample in bin_id_by_sample.keys():
        results.append(taxonomy_groupby(metadata, sample, bin_id_by_sample, taxonomic_dataframe, taxonomic_rank, taxonomic_choice))
    final_results = pd.concat(results, join="outer", ignore_index=True)
    return final_results


def get_taxonomy_size(sample_data: pd.DataFrame, taxonomic_dataframe: pd.DataFrame, only_metabolic_model_size: bool = False):
    """Return the numbers of different species in one sample.

    Args:
        sample_data (pd.DataFrame): dataframe of one sample.
        taxonomic_dataframe (pd.DataFrame): taxonomic dataframe of all samples.
        only_metabolic_model_size (bool, optional): Return only the number of individual in the sample. Defaults to False.

    Returns:
        int: Number of individual or number of different species.
    """
    if only_metabolic_model_size:
        taxonomy_size = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(sample_data["smplID"])]
        return len(taxonomy_size)
    taxonomy_size = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(sample_data["smplID"])]
    taxonomy_size = taxonomy_size[["mgs", "Genus"]]
    taxonomy_size = taxonomy_size.groupby(["Genus"]).count()
    taxonomy_size = taxonomy_size.reset_index()
    return len(taxonomy_size)


def add_row(df: pd.DataFrame, row: list):
    df.loc[len(df)] = row


def get_threshold_value(df: pd.DataFrame, threshold: int = 0, transpose: bool = False) -> dict:
    """Take a matrix with percentage as value instead of a count and return only value above threshold.
    value must be row indexed and sample in columns.

    Args:
        df (pd.DataFrame): Dataframe to extract values from.
        threshold (int, optional): The threshold used to extract value (by default all percentage above but not equal 0 are exctrated). Defaults to 0.
        transpose (bool, optional): True if the dataframe needs to be transpose (when value in columns and sample in rows). Defaults to False.

    Returns:
        dict: A dictionary, for which key correspond to a sample and value in a nested list of value above threshold and their index. return{sample1 : [[value1,value2],[index_value1,index_value2]]}
    """
    if transpose:
        df = df.T
    results = {}
    for sample in df:
        res = df.loc[df[sample] > 1]
        res = res[sample]
        results[sample] = res
    return results


def get_files(file_name: str, path: str, with_directory_name: bool = True):
    """Retrieve files from the name within the directory path given in CLI.

    Args:
        file_name (str): Name of the file
        path (str): directory path (-d from cli) from which os.walk begin the search.
        with_directory_name (bool, optional): If True, also return the directory name containing the file. Defaults to True.

    Returns:
        _type_: _description_
    """
    result = []
    for root, dirs, files in os.walk(path):
        if file_name in files:
            if with_directory_name:
                result.append(
                    [os.path.join(root, file_name), os.path.basename(os.path.dirname(os.path.dirname(os.path.join(root, file_name))))]
                )
            else:
                result.append(os.path.join(root, file_name))
    return result


def open_added_value(file_name, path):
    for root, dirs, files in os.walk(path):
        if file_name in files:
            added_file = sbml_to_classic(open_json(os.path.join(root, file_name))["addedvalue"])
            return added_file


def open_json(file_path):
    with open(file_path) as file:
        file_data = json.load(file)
    return file_data


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


def get_scopes(file_name, path) -> pd.DataFrame:
    for root, dirs, files in os.walk(path):
        if file_name in files:
            scope_matrix = open_tsv(os.path.join(root, file_name), convert_cpd_id=True, rename_columns=True)
            return scope_matrix


def is_indexed_by_id(df: pd.DataFrame):
    if df.index.name == "smplID":
        return True
    else:
        return False


def get_columns_index(df: pd.DataFrame, key_list):
    index_list = [0]
    for k in key_list:
        index_list.append(df.columns.get_loc(k))
    return index_list


def sbml_to_classic(compounds_list):
    uncoded = []
    for coded in compounds_list:
        id, id_type, compart = cfci(coded)
        new_value = str(id)+"["+str(compart)+"]"
        uncoded.append(new_value)
    return uncoded


def contribution_processing(file_opened: dict):
    for key in file_opened.keys():
        for second_key in file_opened[key].keys():
            file_opened[key][second_key] = sbml_to_classic(file_opened[key][second_key])
    return file_opened


def get_contributions(file_name, path):
    for root, dirs, files in os.walk(path):
        if file_name in files:
            contributions_file = open_json(os.path.join(root, file_name))
            contributions_file = contribution_processing(contributions_file)
            return contributions_file

def retrieve_all_sample_data(sample, path):
    """Retrieve iscope, cscope, added_value and contribution_of_microbes files in the path given using os.listdir().

    Args:
        path (str): Directory path

    Returns:
        dict: Return a nested dict object where each key is a dictionnary of a sample. The key of those second layer dict [iscope, cscope, advalue, contribution] give acces to these files.
    """
    sample_directory_path = os.path.join(path, sample)
    if os.path.isdir(sample_directory_path):

        cscope_dataframe = get_scopes("rev_cscope.tsv", sample_directory_path)
        cscope_total_production = cscope_dataframe.apply(lambda column: column.to_numpy().sum(),axis=0)
        cscope_total_production.name = sample

        # all_sample_data[sample]["iscope"] = get_scopes("rev_iscope.tsv", os.path.join(path, sample))
        # all_sample_data[sample]["advalue"] = open_added_value("addedvalue.json", os.path.join(path, sample))
        # all_sample_data[sample]["contribution"] = get_contributions("contributions_of_microbes.json", os.path.join(path, sample))
    else:
        return None, None, None

    return cscope_dataframe, sample, cscope_total_production


def multiprocess_retrieve_data(path):
    retrieve_data = partial(retrieve_all_sample_data, path=path)
    nb_cpu = cpu_count() - 1
    if not type(nb_cpu) == int or nb_cpu < 1:
        nb_cpu = 1
    pool = Pool(nb_cpu)
    results_list = pool.map(retrieve_data,[sample for sample in os.listdir(path)])
    pool.close()
    pool.join
    cpd_producers = []
    all_data = {}
    for df, smpl, cpd_prod in results_list:
        if not df is None: 
            cpd_producers.append(cpd_prod)
            all_data[smpl] = {}
            all_data[smpl]["cscope"] = df

    cpd_producers = pd.concat(cpd_producers,axis=1).T
    cpd_producers.fillna(0,inplace=True)
    cpd_producers = cpd_producers.astype(int)
    cpd_producers.index.name = "smplID"

    return all_data, cpd_producers


def merge_metadata_with_df(main_dataframe, metadata):
    return pd.merge(metadata, main_dataframe, how="left")

def merge_df(left_df, right_df, how: str = "left"):
    data = left_df.iloc[:, 0]
    filter = right_df.loc[right_df["mgs"].isin(data)]
    return filter

def multiply_production_abundance(row: pd.Series, abundance_matrix: pd.DataFrame, sample_id):
    row = row.astype(float)
    abundance_value = abundance_matrix.at[row.name, sample_id]
    row *= abundance_value
    return row

def build_main_dataframe(sample_data: dict):
    all_series = []
    for sample in sample_data.keys():
        current_sample_df = sample_data[sample]["cscope"]
        serie_index = current_sample_df.columns.values
        serie_data = []
        for i in range(len(serie_index)):
            serie_data.append(1)
        all_series.append(pd.Series(data=serie_data,index=serie_index,name=sample))

    results = pd.concat(all_series, axis=1).T
    results.fillna(0,inplace=True)
    results = results.astype(int)
    results.index.name = "smplID"
    
    return results



# @benchmark_decorator
def build_df(dir_path, metadata, abundance_path: str = None, taxonomic_path: str = None):
    """
    Extract community scopes present in directory from CLI then build a single dataframe from the metabolites produced by each comm_scopes.

    Args:
        dir_path (str): Directory path containing comm scopes
        metadata (tsv file): tsv file containing the metadata of the scopes. The number of row must be equal to the number of comm_scopes given in dir_path.

    Returns:
        global_data: dict
        sample_data: dict
        abundance_data: pandas dataframe
    """
    if not is_valid_dir(dir_path):
        print(dir_path, "Not valid directory")
        sys.exit(1)
    
    # if not is_valid_file(metadata):
    #     print(metadata, "Not valid file metadata")
    #     sys.exit(1)

    all_data = {}

    all_data["sample_data"], cpd_producers = multiprocess_retrieve_data(dir_path)

    main_df = build_main_dataframe(all_data["sample_data"])

    all_data["metadata"] = open_tsv(metadata)

    all_data["main_dataframe"] = main_df

    all_data["producers_long_format"] = producer_long_format(all_data["main_dataframe"],all_data["metadata"], cpd_producers, all_data["metadata"].columns[1:])
    # all_data["producers_long_format_new_method"] = producer_long_format_new_method(all_data["main_dataframe"],all_data["metadata"], cpd_producers, all_data["metadata"].columns[1:])
    if abundance_path is not None:
        try:
            raw_abundance_file = open_tsv(abundance_path)
            normalised_abundance_dataframe = relative_abundance_calc(raw_abundance_file, all_data["sample_data"])
        except Exception as e:
            print("Abundance process went wrong.",e)
            normalised_abundance_dataframe = None
    else:
        normalised_abundance_dataframe = None

    if taxonomic_path is not None:
        try:
            raw_taxonomic_data = open_tsv(taxonomic_path)
            long_taxonomic_data = taxonomic_data_long_format(
                    raw_taxonomic_data, get_bin_list(all_data["sample_data"]), all_data["metadata"].columns[1:],all_data["metadata"]
                )
        except Exception as e:
            print("Taxonomy process went wrong.", e)
            long_taxonomic_data = None
    else:
        long_taxonomic_data = None


    total_production_dataframe = total_production_by_sample(all_data["main_dataframe"], all_data["sample_data"], all_data["metadata"], normalised_abundance_dataframe)

    return all_data, normalised_abundance_dataframe, long_taxonomic_data, total_production_dataframe


def list_to_boolean_serie(model_list: list, with_quantity: bool = True):
    results = {}
    value = 1
    for model in model_list:
        if model not in results.keys():
            results[model] = value
        else:
            if with_quantity:
                results[model] += value
    return pd.Series(results)

def taxonomy_matrix_build(taxonomic_df: pd.DataFrame, binlist_by_id: dict):
    rank = "s"
    id_col = "mgs"
    all_series = {}

    for sample in binlist_by_id.keys():
        res = taxonomic_df.loc[taxonomic_df[id_col].isin(binlist_by_id[sample])][rank]
        all_series[sample] = list_to_boolean_serie(res)

    matrix = pd.DataFrame(all_series)
    matrix.fillna(0, inplace=True)
    matrix = matrix.T
    matrix.index.name = "smplID"
    matrix.reset_index(inplace=True)
    # !!! NAME OF ID isnt SMPLID !!! CAN LEAD TO DRAMA
    return matrix

def taxonomic_data_long_format(taxonomic_df: pd.DataFrame, binlist_by_id: dict, metadata_label: list, metadata: pd.DataFrame):
    """
    Produce long format taxonomic dataframe for plot purpose.
    Returns:
        Dataframe: Dataframe in long format
    """

    df = taxonomy_matrix_build(taxonomic_df, binlist_by_id)

    if is_indexed_by_id(df):
        df.reset_index()

    df = df.melt("smplID", var_name="Taxa", value_name="Quantity")
    brand_new_df = {}

    # Assign all metadata a new columns in df.
    for factor in metadata_label:
        brand_new_df[factor] = add_factor_column(metadata, df["smplID"], factor)
    df = df.assign(**brand_new_df)
    df = df.astype(str)
    df["smplID"] = df["smplID"].astype("category")
    df["Quantity"] = df["Quantity"].astype(float)
    df["Nb_taxon"] = df["smplID"].apply(lambda row: search_long_format(row, df))

    return df

def search_long_format(id_value, df):
    value = df.loc[(df["smplID"] == id_value) & (df["Quantity"] != 0)]["Taxa"].unique()
    return len(value)

def get_cpd_quantity(sample_data: dict, sample_id: str):
    results = {}
    for cpd in sample_data[sample_id]["cscope"].columns[1:]:
        results[cpd] = sample_data[sample_id]["cscope"][cpd].sum()
    return pd.Series(results, name=sample_id)

def total_production_by_sample_and_compound(main_df: pd.DataFrame, sample_data: dict):
    prod_df = []
    if not is_indexed_by_id(main_df):
        main_df = main_df.set_index("smplID",drop=True)

    for sample in main_df.index:
        prod_df.append(get_cpd_quantity(sample_data, sample))
    
    final_df = pd.concat(prod_df, axis=1)
    final_df = final_df.fillna(0)
    final_df = final_df.T
    return final_df

def add_factor_column(metadata, serie_id, factor_id):
    if not is_indexed_by_id(metadata):
        metadata = metadata.set_index("smplID", drop=True)
    new_col = []
    for value in serie_id:
        new_col.append(str(metadata.at[value, factor_id]))
    return new_col

def total_production_by_sample(main_dataframe: pd.DataFrame, sample_data: dict, metadata_dataframe: pd.DataFrame, abundance_matrix: pd.DataFrame = None):
    boolean_production_df = main_dataframe.copy()
    if not is_indexed_by_id(boolean_production_df):
        boolean_production_df.set_index("smplID",inplace=True,drop=True)
    boolean_production_df["Total_production"] = boolean_production_df.apply(lambda row: row.to_numpy().sum(), axis=1)
    results = pd.DataFrame(boolean_production_df["Total_production"])

    if abundance_matrix is not None:
        abundance_production_df = abundance_matrix.copy()
        if not is_indexed_by_id(boolean_production_df):
            abundance_production_df.set_index("smplID",inplace=True,drop=True)
        abundance_production_df["Total_abundance_weighted"] = abundance_production_df.apply(lambda row: row.to_numpy().sum(), axis=1)
        abundance_production_df = abundance_production_df["Total_abundance_weighted"]

        results = pd.concat([results,abundance_production_df], axis=1)
    metadata_dict = {}
    # Makes all metadata columns
    for factor in metadata_dataframe.columns[1:]:
        metadata_dict[factor] = add_factor_column(metadata_dataframe, boolean_production_df.index, factor)

    results = results.assign(**metadata_dict)
    results.reset_index(inplace=True)
    results["smplID"] = results["smplID"].astype("category")
    return results

def processing_number_producers(df,cpd):
    df["Nb_producers"] = df.apply(lambda row: cpd[row["smplID"]][row["Compound"]],axis=1)
    return df

def producer_long_format_new_method(main_dataframe: pd.DataFrame, metadata: pd.DataFrame, cpd_producers: pd.DataFrame, metadata_label: list):
    main_dataframe = main_dataframe.reset_index()
    main_dataframe = main_dataframe.melt("smplID",var_name="Compound",value_name="Value")
    main_dataframe = main_dataframe.loc[main_dataframe["Value"] != 0]

    nb_cpu = cpu_count() - 1
    print(nb_cpu, "CPU available")
    pool = Pool(nb_cpu)
    all_sub_dataframes = []
    nb_prod_count = partial(processing_number_producers,cpd=cpd_producers)

    for sample in main_dataframe["smplID"].unique():
        tmp_df = main_dataframe.loc[main_dataframe["smplID"] == sample]
        all_sub_dataframes.append(tmp_df)

    res_sub_dataframe = pool.map(nb_prod_count,all_sub_dataframes)
    pool.close()
    pool.join()

    main_dataframe = pd.concat(res_sub_dataframe)

    metadata_dict = {}
    # Makes all metadata columns
    for factor in metadata_label:
        metadata_dict[factor] = add_factor_column(metadata, main_dataframe["smplID"], factor)

    main_dataframe = main_dataframe.assign(**metadata_dict)
    main_dataframe["smplID"] = main_dataframe["smplID"].astype("category")
    
    return main_dataframe

def producer_long_format(main_dataframe: pd.DataFrame, metadata: pd.DataFrame, cpd_producers: dict, metadata_label: list):
    main_dataframe = main_dataframe.reset_index()
    main_dataframe = main_dataframe.melt("smplID",var_name="Compound",value_name="Value")
    main_dataframe = main_dataframe.loc[main_dataframe["Value"] != 0]

    main_dataframe["Nb_producers"] = main_dataframe.apply(lambda row: cpd_producers.at[row["smplID"],row["Compound"]],axis=1)

    metadata_dict = {}
    # Makes all metadata columns
    for factor in metadata_label:
        metadata_dict[factor] = add_factor_column(metadata, main_dataframe["smplID"], factor)

    main_dataframe = main_dataframe.assign(**metadata_dict)
    main_dataframe["smplID"] = main_dataframe["smplID"].astype("category")
    
    return main_dataframe


def stat_on_plot(data: dict, layer: int):
    """Apply Wilcoxon or Mann-Whitney test on each pair of a dataframe.

    Args:
        data (dict): Dictionnary of each pair to test.
        layer (int): Number of layer to the dict. 2 for nested dict.

    Returns:
        Dataframe: Dataframe of test's results.
    """
    
    if layer == 1:
        error_log = []
        res = pd.DataFrame(columns=["group_1","n1","group_2","n2","test_value","p_value","Significance","Test"])
        for pair_1 in data.keys():
            for pair_2 in data.keys():
                if pair_1 != pair_2:
                    if not pair_2 in res["group_1"].tolist():
                        pair1_data = data[pair_1]["data"]
                        pair2_data = data[pair_2]["data"]
                        n1, n2 = data[pair_1]["n_data"], data[pair_2]["n_data"]
                        if len(pair1_data) != 0 and len(pair2_data) != 0:
                            if len(pair1_data) == len(pair2_data):
                                try:
                                    test_value, p_value = stats.wilcoxon(pair1_data, pair2_data)
                                    test_type = "Wilcoxon"
                                except:
                                    error_log.append([pair_1,pair1_data,pair_2,pair2_data])
                                    continue
                            else:
                                try:
                                    test_value, p_value = stats.mannwhitneyu(pair1_data, pair2_data)
                                    test_type = "Mann-Whitney"
                                except:
                                    error_log.append([pair_1,pair1_data,pair_2,pair2_data])
                                    continue
                            if p_value >= 0.05:
                                symbol = "ns"
                            elif p_value >= 0.01:
                                symbol = "*"
                            elif p_value >= 0.001:
                                symbol = "**"
                            else:
                                symbol = "***"
                            new_row = {"group_1": pair_1,"n1": n1, "group_2":pair_2,"n2": n2, "test_value": test_value, "p_value": p_value,"Significance": symbol,"Test": test_type}
                            res.loc[len(res)] = new_row

    if layer == 2:
        error_log = []
        res = pd.DataFrame(columns=["Axis","group_1","n1","group_2","n2","test_value","p_value","Significance","Test"])
        for current_layer in data.keys():
            for pair_1 in data[current_layer].keys():
                for pair_2 in data[current_layer].keys():
                    # Don't test same pair
                    if pair_1 != pair_2:
                        pair1_data = data[current_layer][pair_1]["data"]
                        pair2_data = data[current_layer][pair_2]["data"]
                        n1, n2 = data[current_layer][pair_1]["n_data"], data[current_layer][pair_2]["n_data"]
                        # Don't test if no value
                        if len(pair1_data) != 0 and len(pair2_data) != 0:
                            # Don't test pair already tested.
                            if len(res.loc[(res["group_1"] == pair_2) & (res["Axis"] == current_layer)] ) == 0:
                                # If len of both pair is same value then Wilcoxon, else Mann Whitney
                                if len(pair1_data) == len(pair2_data):
                                    try:
                                        test_value, p_value = stats.wilcoxon(pair1_data, pair2_data)
                                        test_type = "Wilcoxon"
                                    except:
                                        error_log.append([current_layer,pair_1,pair1_data,current_layer,pair_2,pair2_data])
                                        continue
                                else:
                                    try:
                                        test_value, p_value = stats.mannwhitneyu(pair1_data, pair2_data)
                                        test_type = "Mann-Whitney"
                                    except:
                                        error_log.append([current_layer,pair_1,pair1_data,current_layer,pair_2,pair2_data])
                                        continue

                                if p_value >= 0.05:
                                    symbol = "ns"
                                elif p_value >= 0.01:
                                    symbol = "*"
                                elif p_value >= 0.001:
                                    symbol = "**"
                                else:
                                    symbol = "***"
                                new_row = {"Axis":current_layer,
                                            "group_1": pair_1,
                                            "n1": n1,
                                            "group_2": pair_2,
                                            "n2": n2,
                                            "test_value": test_value,
                                            "p_value": p_value,
                                            "Significance": symbol,
                                            "Test": test_type
                                            }
                                res.loc[len(res)] = new_row
    print("At least: ",len(error_log)," errors occured.")
    print(error_log)
    return res
