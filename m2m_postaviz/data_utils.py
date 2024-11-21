import json
import os
import os.path
import tarfile
import time
import sys
import pickle

import numpy as np
import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci
from scipy import stats
from multiprocessing import Pool
from multiprocessing import cpu_count
from functools import partial
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import squareform,pdist
from statsmodels.stats.multitest import multipletests


def pickle_write(path, obj):
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def pickle_read(path):
    with open(path, "rb") as f:
        obj = pickle.load(f)
    return obj


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


def has_only_unique_value(dataframe: pd.DataFrame, input1, input2: str = "None"):
    """
    Return True if df column value contain only one value per input, False otherwise.

    Args:
        dataframe (pd.DataFrame): _description_
        column_value (_type_): _description_
        input1 (_type_): _description_
        input2 (str, optional): _description_. Defaults to "None".
    """
    nb_row = len(dataframe)

    if input2 == "None":
        return True if nb_row == len(dataframe[input1].unique()) else False
        
    else:
        return True if nb_row == len(dataframe[input1].unique()) and nb_row == len(dataframe[input2].unique()) else False


def relative_abundance_calc(abundance_matrix: pd.DataFrame, sample_data: dict) -> pd.DataFrame:
    """Generate a second main_dataframe with the production based on weight from the abundance matrix.

    Args:
        abundance_matrix (pd.DataFrame): abundance matrix given in input.
        sample_data (dict): Dictionnary of sample's cscopes.

    Raises:
        RuntimeError: If more than one column of type other than INT.

    Returns:
        Dataframe: production dataframe with sample in rows and compounds in column. Weighted by abundance.
    """
    smpl_norm_abundance = []
    smpl_norm_index = []

    # Checking if all column are INT type, if one is not its used as index, if more than 1 raise RuntimeERROR.
    str_filter = abundance_matrix.select_dtypes(include=["string","object","category"])

    if len(str_filter.columns) == 1:
        index_column = str_filter.columns.values[0]
        print(index_column, "column used as index")
        abundance_matrix.set_index(index_column,drop=True,inplace=True)
        
    elif len(str_filter.columns) > 1:
        raise RuntimeError("More than one non-numeric columns in abundance dataframe.")
    
    # Normalisation
    abundance_matrix_normalised = abundance_matrix.apply(lambda x: x / x.sum(), axis=0)

    # For all sample's cscopes, multiply each row (bin's production) by the normalised abundance matrix. 
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
        if cscope_dataframe is None:
            return None, sample

    else:
        return None, sample

    return cscope_dataframe, sample


def producers_by_compounds_and_samples_multi(all_data: dict, metadata: pd.DataFrame):

    if not bool(all_data):
        raise Exception("Sample data empty.")

    cpu_available = cpu_count() - 1
    if not type(cpu_available) == int or cpu_available < 1:
        cpu_available = 1
    pool = Pool(cpu_available)
    all_producers = pool.starmap(individual_producers_processing,[(all_data[sample]["cscope"], sample) for sample in all_data.keys()])
    pool.close()
    pool.join()

    res = pd.concat(all_producers,axis=1).T
    res.fillna(0,inplace=True)
    res.index.name = "smplID"
    res.reset_index(inplace=True)
    res = res.merge(metadata,'inner',"smplID")

    return res


def individual_producers_processing(sample_cscope: pd.DataFrame , sample: str):
    """sums the count of metabolites produced in the sample's cscope.

    Args:
        sample_cscope (pd.DataFrame): Sample cscope dataframe
        sample (str): Sample ID

    Returns:
        pd.Series: Pandas serie with all metabolites columns sum
    """
    serie_value = []
    serie_index = []

    for i in range(len(sample_cscope.columns)):
        serie_index.append(sample_cscope.columns[i])
        serie_value.append(sample_cscope[sample_cscope.columns[i]].to_numpy().sum())
    return pd.Series(serie_value,index=serie_index,name=sample)


def multiprocess_retrieve_data(path):
    """Open all directories given in -d path input. Get all cscopes tsv and load them in memory as pandas
    dataframe. 

    Args:
        path (str): Path of directory

    Returns:
        dict: sample_data dictionnary
    """
    retrieve_data = partial(retrieve_all_sample_data, path=path)

    nb_cpu = cpu_count() - 1
    if not type(nb_cpu) == int or nb_cpu < 1:
        nb_cpu = 1
    pool = Pool(nb_cpu)
    results_list = pool.map(retrieve_data,[sample for sample in os.listdir(path)])
    
    pool.close()
    pool.join()

    all_data = {}
    for df, smpl in results_list:
        if not df is None: 
            all_data[smpl] = {}
            all_data[smpl]["cscope"] = df

    sample_info = {}
    sample_info["bins_list"] = []
    sample_info["bins_count"] = {}
    sample_info["bins_sample_list"] = {}

    for sample in all_data.keys():

        dataframe = all_data[sample]["cscope"]
        all_bins_in_sample = dataframe.index.tolist()

        sample_info["bins_list"] = sample_info["bins_list"] + all_bins_in_sample

        for bin in all_bins_in_sample:
            
            if not bin in sample_info["bins_count"]:
                sample_info["bins_count"][bin] = 0
            
            sample_info["bins_count"][bin] += 1

            if not bin in sample_info["bins_sample_list"]:
                sample_info["bins_sample_list"][bin] = []
            
            sample_info["bins_sample_list"][bin].append(str(sample))

    # Remove duplicate from list
    sample_info["bins_list"] = list(dict.fromkeys(sample_info["bins_list"]))

    return sample_info, all_data


def melt_df_multi(dataframe: pd.DataFrame) -> pd.DataFrame:
    dataframe.reset_index(inplace=True)
    return dataframe.melt("smplID",var_name="Compound",value_name="Value")


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


def build_df(dir_path, metadata_path: str, abundance_path: str = None, taxonomic_path: str = None, save_path: str = None):
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
        print(dir_path, "Sample directory path is not a valid directory")
        sys.exit(1)

    sample_info, sample_data, metadata, main_dataframe, normalised_abundance_dataframe, abundance_file, taxonomic_dataframe, producers_dataframe, total_production_dataframe, pcoa_dataframe, bin_dataframe_loaded, bin_dataframe = check_for_files(save_path)
    # main_dataframe, metadata, producers_dataframe, normalised_abundance_dataframe, taxonomy, total_production_dataframe, pcoa_dataframe = load_hdf5_datafames(save_path)

    if metadata is None: 
        print("Open metadata dataframe.")
        metadata = open_tsv(metadata_path)

    if not bool(sample_data): 
        print("Fetching sample's data...")
        sample_info, sample_data = multiprocess_retrieve_data(dir_path)

    sample_info = iscope_production(dir_path=dir_path, sample_info_dict=sample_info)

    if producers_dataframe is None:
        print("Building metabolite production dataframe...")
        producers_dataframe = producers_by_compounds_and_samples_multi(sample_data,metadata) 

    if main_dataframe is None:
        print("Building main dataframe...")
        main_dataframe = build_main_dataframe(sample_data)

    if normalised_abundance_dataframe is None:
        if abundance_path is not None:
            try:
                abundance_file = open_tsv(abundance_path)
                normalised_abundance_dataframe = relative_abundance_calc(abundance_file, sample_data)
            except Exception as e:
                print("Abundance process went wrong.",e)
                abundance_file = None
                normalised_abundance_dataframe = None
        else:
            abundance_file = None
            normalised_abundance_dataframe = None

    if taxonomic_dataframe is None:
        if taxonomic_path is not None:
            try:
                taxonomic_dataframe = taxonomy_processing(taxonomic_path)

            except Exception as e:
                taxonomic_dataframe = None
        else:
            taxonomic_dataframe = None

    if total_production_dataframe is None:
        print("Building global production dataframe...")
        total_production_dataframe = total_production_by_sample(main_dataframe, metadata, normalised_abundance_dataframe)

    if pcoa_dataframe is None:
        print("Running pcoa with main dataframe...")
        pcoa_dataframe = pcoa_alternative_method(main_dataframe, metadata)

    if not bin_dataframe_loaded:
        print("Building bin dataframe...")
        bin_dataframe_build(sample_info, sample_data, metadata, abundance_file, taxonomic_dataframe, save_path)

    save_all_dataframe(sample_info, sample_data, metadata, main_dataframe, normalised_abundance_dataframe, taxonomic_dataframe, producers_dataframe, total_production_dataframe, pcoa_dataframe, save_path, abundance_file)
    file_format = "tsv"

    taxonomy_provided = False if taxonomic_dataframe is None else True
    abundance_provided = False if normalised_abundance_dataframe is None else True

    return file_format, taxonomy_provided, abundance_provided


def save_dataframe_hdf_format(metadata, main_dataframe, norm_abundance_df: pd.DataFrame = None, long_taxo_df: pd.DataFrame = None, producers_dataframe: pd.DataFrame = None, total_production_df: pd.DataFrame = None, pcoa_dataframe: pd.DataFrame = None, savepath: str = None):
    """Save every dataframe to save_path input.

    Args:
        all_data (dict): all_sample dataframe, metadata, main_dataframe and producer_dataframe
        norm_abundance_df (Dataframe): abundance dataframe normalised
        long_taxo_df (Dataframe): taxonomic dataframe in long format
        total_production_df (Dataframe): total production dataframe
        savepath (str): path to save all files.
    """
    if savepath is None:
        print("save_path is None, data can't be saved into HDF5 format.")
        sys.exit(1)
    
    filename = os.path.join(savepath,"postaviz_dataframes.h5")

    if not os.path.exists(savepath):
        try:
            os.makedirs(savepath)
        except FileExistsError:
            pass
    
    if os.path.isfile(filename):
        print(f"Removing existing files at {filename}")
        try:
            os.remove(filename)
        except OSError:
            pass

    pd.set_option('io.hdf.default.format', 'table')
    
    store = pd.HDFStore(filename)
    
    if norm_abundance_df is not None:
        store['normalised_abundance_dataframe'] = norm_abundance_df
    else:
        print("Unable to save normalised_abundance_dataframe")
    
    if total_production_df is not None:
        store['global_production_dataframe'] = total_production_df
    else:
        print("Unable to save global_production_dataframe")
    
    if producers_dataframe is not None:
        store['metabolite_production_dataframe'] = producers_dataframe
    else:
        print("Unable to save metabolite_production_dataframe")
    
    if pcoa_dataframe is not None:
        store['pcoa_dataframe'] = pcoa_dataframe
    else:
        print("Unable to save pcoa_dataframe")
    
    if long_taxo_df is not None:
        store['taxonomic_dataframe'] = long_taxo_df
    else:
        print("Unable to save taxonomic_dataframe")
    
    if main_dataframe is not None:
        store['main_dataframe'] = main_dataframe
    else:
        print("Unable to save main_dataframe")
        
    if metadata is not None:
        store['metadata'] = metadata
    else:
        print("Unable to save metadata")
    
    store.close()

    #### Usage of sample data disabled for now. ####

    # store_sample = pd.HDFStore(os.path.join(savepath,"postaviz_samples.h5"))

    # import warnings
    # warnings.filterwarnings("ignore")

    # for sample in global_data["sample_data"].keys():
    #     store_sample[sample] = global_data["sample_data"][sample]["cscope"]

    # store_sample.close()

    return


def load_hdf5_datafames(path: str):
    keys = ["main_dataframe", "metadata", "metabolite_production_dataframe", "normalised_abundance_dataframe", "long_taxonomic_data", "global_production_dataframe", "pcoa_dataframe"]
    results = []

    with pd.HDFStore(path=os.path.join(path,"postaviz_dataframes.h5"),mode='r') as storage:
        print("Loading following dataframes... ", storage.keys())

        for k in keys:

            try:
                results.append(storage[k])

            except KeyError as e:
                print(e)
                results.append(None)

    return results


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


def add_factor_column(metadata, serie_id, factor_id):
    if not is_indexed_by_id(metadata):
        metadata = metadata.set_index("smplID", drop=True)
    new_col = []
    for value in serie_id:
        new_col.append(str(metadata.at[value, factor_id]))
    return new_col


def total_production_by_sample(main_dataframe: pd.DataFrame, metadata_dataframe: pd.DataFrame, abundance_matrix: pd.DataFrame = None):
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

    results.reset_index(inplace=True)
    results = results.merge(metadata_dataframe,'inner','smplID')
    # results["smplID"] = results["smplID"].astype("category")

    return results


def preprocessing_for_statistical_tests(dataframe: pd.DataFrame, y_value, input1, input2 = None, multipletests: bool = False, multipletests_method: str = "bonferroni"):
    """Create dataframe for each y_value in the list, to separate them and use wilcoxon_man_whitney function.
    Concat all results into one dataframe.

    Args:
        dataframe (pd.DataFrame): Dataframe to test.
        y_value (_type_): list of columns labels to separata into several dataframe. Must be of lenght 1 at least. 
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


def wilcoxon_man_whitney(dataframe: pd.DataFrame, y, first_factor: str, second_factor: str = None, multiple_correction: bool = False, correction_method: str = "hs"):
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
    for name in sub_dataframes.keys():

        if second_factor is None: # One input selected

            for name2 in sub_dataframes.keys():
                
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


def check_for_files(save_path = None):
    """Open dataframe produced from previous postaviz session in savepath directory.

    Args:
        save_path (str): save_path directory path

    Returns:
        Tuple: Tuple of dataframe (dict, dataframe, dataframe, dataframe)
    """
    
    metadata = None
    main_dataframe = None
    producers_dataframe = None
    sample_info = {}
    sample_data = {}
    total_production_dataframe = None
    taxonomic_dataframe = None
    normalised_abundance_dataframe = None
    abundance_file = None
    pcoa_dataframe = None
    bin_dataframe = None
    bin_dataframe_loaded = False

    if save_path is None or not is_valid_dir(save_path):
        print("save_path is none, ignoring load data function.")
        return sample_info, sample_data, metadata, main_dataframe, normalised_abundance_dataframe, abundance_file, taxonomic_dataframe, producers_dataframe, total_production_dataframe, pcoa_dataframe, bin_dataframe_loaded, bin_dataframe

    for root, dirname, filename in os.walk(save_path):

        if "all_samples_dataframe_postaviz" in dirname:
            dirpath = os.path.join(root,"all_samples_dataframe_postaviz")
            for sample in os.listdir(dirpath):
                if sample.endswith("_cscope.pkl"):
                    try:
                        sample_name = sample.split("_", 1)[0]
                        sample_data[sample_name] = {}
                        sample_data[sample_name]["cscope"] = pd.read_pickle(os.path.join(dirpath, sample))
                    except Exception as e:
                        print("Read pickle file: ",sample,"in", os.path.join(dirpath,sample)," failed.", "\n",e)
                        continue

        if "sample_info.json" in filename:
            sample_info_path = os.path.join(root,"sample_info.json")
            with open(sample_info_path) as f:
                sample_info = json.load(f)

        if "metadata_dataframe_postaviz.tsv" in filename:
            metadata = open_tsv(os.path.join(root, "metadata_dataframe_postaviz.tsv"))

        if "main_dataframe_postaviz.tsv" in filename:
            main_dataframe = open_tsv(os.path.join(root, "main_dataframe_postaviz.tsv"))

        if "producers_dataframe_postaviz.tsv" in filename:
            producers_dataframe = open_tsv(os.path.join(root, "producers_dataframe_postaviz.tsv"))

        if "total_production_dataframe_postaviz.tsv" in filename:
            total_production_dataframe = open_tsv(os.path.join(root, "total_production_dataframe_postaviz.tsv"))

        if "taxonomic_dataframe_postaviz.tsv" in filename:
            taxonomic_dataframe = pd.read_csv(os.path.join(root, "taxonomic_dataframe_postaviz.tsv"), sep="\t", index_col=0)

        if "normalised_abundance_dataframe_postaviz.tsv" in filename:
            normalised_abundance_dataframe = open_tsv(os.path.join(root, "normalised_abundance_dataframe_postaviz.tsv"))
        
        if "abundance_file.tsv" in filename:
            abundance_file = pd.read_csv(os.path.join(root, "abundance_file.tsv"), sep="\t", index_col=0)
        
        if "pcoa_dataframe_postaviz.tsv" in filename:
            pcoa_dataframe = open_tsv(os.path.join(root, "pcoa_dataframe_postaviz.tsv"))

        if "bin_dataframe_chunk_1.parquet.gzip" in filename:
            bin_dataframe_loaded = True

    return sample_info, sample_data, metadata, main_dataframe, normalised_abundance_dataframe, abundance_file, taxonomic_dataframe, producers_dataframe, total_production_dataframe, pcoa_dataframe, bin_dataframe_loaded, bin_dataframe


def save_all_dataframe(sample_info,
                        sample_data,
                        metadata,
                        main_dataframe,
                        normalised_abundance_dataframe,
                        taxonomic_dataframe, producers_dataframe,
                        total_production_dataframe,
                        pcoa_dataframe,
                        savepath,
                        raw_abundance_file):
    """Save every dataframe to save_path input.

    Args:
        all_data (dict): all_sample dataframe, metadata, main_dataframe and producer_dataframe
        norm_abundance_df (Dataframe): abundance dataframe normalised
        long_taxo_df (Dataframe): taxonomic dataframe in long format
        total_production_df (Dataframe): total production dataframe
        savepath (str): path to save all files.
    """

    if savepath is None:
        print("save_path is None, data will not be saved.")
        return

    # Sample cscopes PICKLE
    data_dir = os.path.join(savepath,"all_samples_dataframe_postaviz")

    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)
    
    for sample in sample_data.keys():
        try:
            filename = str(sample+"_cscope.pkl")
            file_path = os.path.join(data_dir,filename)
            sample_data[sample]["cscope"].to_pickle(file_path)
        except Exception as e:
            print("Saving ERROR",e)

    # Sample_info JSON
    if os.path.isfile(os.path.join(savepath,"sample_info.json")):
        print(os.path.join(savepath,"sample_info.json"), "directory already exist in save_path.")
        pass
    else:
        with open(os.path.join(savepath,"sample_info.json"), 'w') as f:
            json.dump(sample_info, f)
    
    # Main_dataframe
    if os.path.isfile(os.path.join(savepath,"main_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"main_dataframe_postaviz.tsv"), " file already exist.")
    else:
        if main_dataframe is not None: main_dataframe.to_csv(os.path.join(savepath,"main_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(main_dataframe) else False)

    # Metadata_dataframe
    if os.path.isfile(os.path.join(savepath,"metadata_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"metadata_dataframe_postaviz.tsv"), " file already exist.")
    else:
        if metadata is not None: metadata.to_csv(os.path.join(savepath,"metadata_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(metadata) else False)
 
    # Producers_dataframe
    if os.path.isfile(os.path.join(savepath,"producers_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"producers_dataframe_postaviz.tsv"), " file already exist.")
    else:
        if producers_dataframe is not None: producers_dataframe.to_csv(os.path.join(savepath,"producers_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(producers_dataframe) else False)

    # Normalised abundance dataframe
    if os.path.isfile(os.path.join(savepath,"normalised_abundance_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"normalised_abundance_dataframe_postaviz.tsv"), " file already exist")
    else:
        if normalised_abundance_dataframe is not None: normalised_abundance_dataframe.to_csv(os.path.join(savepath,"normalised_abundance_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(normalised_abundance_dataframe) else False)
 
    # Taxonomic dataframe
    if os.path.isfile(os.path.join(savepath,"taxonomic_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"taxonomic_dataframe_postaviz.tsv"), " file already exist")
    else:
        if taxonomic_dataframe is not None:
            taxonomic_dataframe.index.name = "binID"
            taxonomic_dataframe.to_csv(os.path.join(savepath,"taxonomic_dataframe_postaviz.tsv"), sep="\t", index=False)
 
    # Total production dataframe
    if os.path.isfile(os.path.join(savepath,"total_production_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"total_production_dataframe_postaviz.tsv"), " file already exist")
    else:
        if total_production_dataframe is not None: total_production_dataframe.to_csv(os.path.join(savepath,"total_production_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(total_production_dataframe) else False)

    # Pcoa dataframe
    if os.path.isfile(os.path.join(savepath,"pcoa_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"pcoa_dataframe_postaviz.tsv"), " file already exist")
    else:
        if pcoa_dataframe is not None: pcoa_dataframe.to_csv(os.path.join(savepath,"pcoa_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(pcoa_dataframe) else False)

    # Abundance dataframe
    if os.path.isfile(os.path.join(savepath,"abundance_file.tsv")):
        print(os.path.join(savepath,"abundance_file.tsv"), " file already exist")
    else:
        if raw_abundance_file is not None:
            raw_abundance_file.index.name = "binID"
            raw_abundance_file.to_csv(os.path.join(savepath,"abundance_file.tsv"),sep="\t", index= True)

    # # Bin_dataframe ##### No longer use since bin_dataframe function save chunks of dataframe
    # if bin_dataframe_loaded:
    #     print(os.path.join(savepath,"bin_dataframe.parquet.gzip"), " file already exist")
    # else:
    #     bin_dataframe.to_parquet(os.path.join(savepath,"bin_dataframe.parquet.gzip"), compression='gzip')

    return


def pcoa_alternative_method(main_dataframe: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Comptute Principal Coordinate Analysis from the main_dataframe given in input. Merge with metadata from the smplID column or index.

    Args:
        main_dataframe (pd.DataFrame): Dataframe from which the pcoa will be made.
        metadata (pd.DataFrame): Metadata dataframe (Must have smplID identifer column for the merge to work)

    Returns:
        pd.DataFrame: Pcoa dataframe with sample ID as index, PC1 and PC2 results and all metadata.
    """
    if not is_indexed_by_id(main_dataframe):
        main_dataframe = main_dataframe.set_index("smplID")

    if  is_indexed_by_id(metadata):
        metadata = metadata.reset_index("smplID")

    dmatrix = main_dataframe.to_numpy()
    dist_m = pdist(dmatrix, "jaccard")
    square_df = squareform(dist_m)
    pcoa_result = pcoa(square_df,number_of_dimensions=2)
    coordinate = pcoa_result.samples

    df_pcoa = coordinate[['PC1','PC2']]
    df_pcoa['smplID'] = main_dataframe.index.to_numpy()

    df_pcoa = df_pcoa.merge(metadata, "inner", "smplID")
    df_pcoa.set_index("smplID",inplace=True)

    return df_pcoa


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


def serie_is_float(ser: pd.Series):

    if np.issubdtype(ser.dtype, np.integer) or np.issubdtype(ser.dtype, np.floating):
        return True
    
    return False


def taxonomy_processing(taxonomy_filepath):
    """Open taxonomy file en process it if in txt format.

    Args:
        taxonomy_filepath (str): TSV or TXT format

    Raises:
        RuntimeError: Wrong file's format

    Returns:
        pd.DataFrame: Pandas dataframe
    """

    if taxonomy_filepath.endswith(".tsv"):
        df = open_tsv(taxonomy_filepath)
        df = df.rename(columns={df.columns[0]: 'mgs'})
        return df

    if not taxonomy_filepath.endswith(".txt"):
        raise RuntimeError("Taxonomy file must be either a txt file or tsv file.")
    
    with open(taxonomy_filepath, 'r') as f:
        lines = f.readlines()

    df = pd.DataFrame(columns=["mgs","kingdom","phylum","class","order","family","genus"])

    del lines[0] # Delete header line

    for row in lines:

        mgs = row.split("\t")[0]
        genus = row.split("\t")[1:]

        k, p, c, o, f, g = genus[0].strip("\n").split(";")

        df.loc[len(df)] = [mgs,k,p,c,o,f,g]

    return df


def iscope_production(dir_path: str, sample_info_dict: dict):

    indiv_scope_path = "indiv_scopes/rev_iscope.tsv"
    sample_info_dict["iscope"] = {}
    # Takes first of like of sample where bin is present then get iscope production via file rev_iscope.tsv 
    for bin in sample_info_dict["bins_sample_list"].keys():

        if bin in sample_info_dict["iscope"]:
            continue

        sample_used = sample_info_dict["bins_sample_list"][bin][0]
        file_path = os.path.join(os.path.join(dir_path, sample_used), indiv_scope_path)
        df = open_tsv(file_path,True,True)
        bin_row = df.loc[bin]

        sample_info_dict["iscope"][bin] = bin_row.index.tolist()

    return sample_info_dict


def bin_dataframe_build(sample_info: dict, sample_data: dict, metadata, abundance_file = None, taxonomy_file = None, savepath = None):
    """Build a large dataframe with all the bins of the different samples as index, the dataframe contain the list of production, abundance, count,
    the metadata and the taxonomic rank associated.

    Args:
        sample_info (dict): _description_
        sample_data (dict): _description_
        metadata (Dataframe): _description_
        abundance_file (Dataframe, optional): _description_. Defaults to None.
        taxonomy_file (Dataframe, optional): _description_. Defaults to None.

    Returns:
        pd.DataFrame: Pandas dataframe
    """
    
    start = time.time()

    ##### Abundance normalisation, give percentage of abundance of bins in samples.
    if abundance_file is not None:

        abundance_matrix_normalised = abundance_file.apply(lambda x: x / x.sum(), axis=0)

    
    ##### Iterate thought the bins key in sample_info_dict to get list of sample where they are present. 

    sample_list = []
    bin_list = []
    for bin in sample_info["bins_sample_list"].keys():
        sample_list += sample_info["bins_sample_list"][bin]
        bin_list.append(bin)

    #   Delete replicate in list
    sample_unique_list = list(dict.fromkeys(sample_list))

    ##### Create Generator to process data by chunks.
    print("Making chunk of sample list...")

    chunk_generator = chunks(sample_unique_list, 250)

    ##### Loop throught generator

    chunk_index = 0
    for current_chunk in chunk_generator:
        start_chunk = time.time()
        chunk_index += 1
        list_of_dataframe = []
        logs = []

    ##### Loop in sample unique list, get their index (where bins are listed) then select row with isin(bins)

        for sample in current_chunk:

            try:
                df = sample_data[sample]["cscope"]
                rows = df.loc[df.index.isin(bin_list)]
                rows.insert(0 , "smplID", sample)
                rows.index.name = "binID"
                list_of_dataframe.append(rows)

            except:
                logs.append(f'No dataframe named {sample} in sample_data dictionnary')

        results = pd.concat(list_of_dataframe)
        results.fillna(0,inplace=True)
        print(f'Chunk {chunk_index} first concat with {sys.getsizeof(results)/1000000000} Gb memory size')

        s = results.apply(lambda row: get_production_list_from_bin_dataframe(row), axis=1)
        s.name = "Production"
        print(f'Chunk {chunk_index} production serie produced with {sys.getsizeof(s)/1000000000} Gb memory size of production serie')

        count = results.drop("smplID", axis=1).apply(np.sum,axis=1,raw=True)
        count.name = "Count"
        print(f'Chunk {chunk_index} count serie produced with {sys.getsizeof(count)/1000000000} Gb memory size')
    # del sample_list
    # del sample_unique_list
    # del bin_list
        results = pd.concat([results["smplID"],count,s],axis=1)
        del count
        del s

        if abundance_file is not None:

            abundance = results.apply(lambda row: abundance_matrix_normalised.at[row.name,row["smplID"]],axis=1)
            abundance.name = "Abundance"

            final_result = pd.concat([results,abundance], axis=1)
            del results

            final_result["Count_with_abundance"] = final_result["Count"] * final_result["Abundance"]

        else:
            
            final_result = results
            del results
            
        final_result = final_result.reset_index().merge(metadata, "inner", "smplID")


        if taxonomy_file is not None:

            ##### Checks if taxonomy file has default index which mean it is not indexed. TEMPORARY until found better way to deal with open/save from -t option OR load taxonomic_df option which return non indexed / indexed df
            if not pd.Index(np.arange(0, len(taxonomy_file))).equals(taxonomy_file.index):
                print("INDEXED !")
                taxonomy_file = taxonomy_file.reset_index()

            mgs_col_taxonomy = taxonomy_file.columns[0]

            final_result = final_result.merge(taxonomy_file, "inner", left_on="binID", right_on=mgs_col_taxonomy)

        ##### Save current chunks into parquet file

        filename = 'bin_dataframe_chunk_'+str(chunk_index)+'.parquet.gzip'
        filepath = os.path.join(savepath,filename)
        if len(final_result) == 0:
            print(f'Chunks {chunk_index} is empty !')

        final_result.to_parquet(filepath, compression='gzip')

        print(f'Chunk {chunk_index} done in {time.time() - start_chunk} with {sys.getsizeof(final_result) / 1000000000} Gb memory size.')
        del final_result


    print("Took: ", time.time() - start, "Before saving")

    if len(logs) != 0:
        with open(os.path.join(savepath,'bin_dataframe_logs.txt'),'w') as log_file:
            for line in logs:
                log_file.write("%s\n" % line)

    return


def get_production_list_from_bin_dataframe(serie: pd.Series) -> list:

    list_of_cpd_produced = []

    for label, value in serie.items():
        if label == "smplID":
            continue
        if value > 0:
            list_of_cpd_produced.append(label)

    return list_of_cpd_produced


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]