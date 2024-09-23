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
from functools import partial
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import squareform,pdist

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
        raise RuntimeError("More than one non-numeric columns.")
    
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

    # Correct melted version sys.getsizeof() == 10053680472 !
    # res = res.melt("smplID",var_name="Compound",value_name="Value")
    # metadata_dict = {}
    # for factor in metadata.columns[1:]:
    #     metadata_dict[factor] = add_factor_column(metadata, res["smplID"], factor)

    # res = res.assign(**metadata_dict)

    return res

def individual_producers_processing(sample_cscope: pd.DataFrame , sample: str):
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

    return all_data

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

    sample_data, metadata, main_dataframe, normalised_abundance_dataframe, taxonomy, producers_dataframe, total_production_dataframe, pcoa_dataframe = check_for_files(save_path)
    # main_dataframe, metadata, producers_dataframe, normalised_abundance_dataframe, taxonomy, total_production_dataframe, pcoa_dataframe = load_hdf5_datafames(save_path)

    if metadata is None: 
        print("Open metadata dataframe.")
        metadata = open_tsv(metadata_path)

    if not bool(sample_data): 
        print("Fetching sample's cscope...")
        sample_data = multiprocess_retrieve_data(dir_path) 

    if producers_dataframe is None:
        print("Building metabolite production dataframe...")
        producers_dataframe = producers_by_compounds_and_samples_multi(sample_data,metadata) 

    if main_dataframe is None:
        print("Building main dataframe...")
        main_dataframe = build_main_dataframe(sample_data)

    if normalised_abundance_dataframe is None and abundance_path is not None:
        try:
            raw_abundance_file = open_tsv(abundance_path)
            normalised_abundance_dataframe = relative_abundance_calc(raw_abundance_file, sample_data)
        except Exception as e:
            print("Abundance process went wrong.",e)
            normalised_abundance_dataframe = None

    if taxonomy is None and taxonomic_path is not None:
        try:
            raw_taxonomic_data = open_tsv(taxonomic_path)
            taxonomy = taxonomic_data_long_format(
                    raw_taxonomic_data, get_bin_list(sample_data), metadata.columns[1:],metadata
                )
        except Exception as e:
            print("Taxonomy process went wrong.", e)
            taxonomy = None
    else:
        taxonomy = None

    if total_production_dataframe is None:
        print("Building global production dataframe...")
        total_production_dataframe = total_production_by_sample(main_dataframe, sample_data, metadata, normalised_abundance_dataframe)

    if pcoa_dataframe is None:
        print("Running pcoa with main dataframe...")
        pcoa_dataframe = pcoa_alternative_method(main_dataframe, metadata)

    try:
        hdf5_file_path = save_dataframe_hdf_format(metadata, main_dataframe, normalised_abundance_dataframe, taxonomy, producers_dataframe, total_production_dataframe, pcoa_dataframe, save_path)

    except Exception as e:
        print(e)
        print("Saving as TSV format instead")
        save_all_dataframe(sample_data , metadata, main_dataframe, normalised_abundance_dataframe, taxonomy, producers_dataframe, total_production_dataframe, pcoa_dataframe, save_path)

    taxonomy_provided = False if taxonomy is None else True
    abundance_provided = False if normalised_abundance_dataframe is None else True

    return hdf5_file_path, taxonomy_provided, abundance_provided


def save_dataframe_hdf_format(metadata, main_dataframe, norm_abundance_df: pd.DataFrame = None, long_taxo_df: pd.DataFrame = None, producers_dataframe: pd.DataFrame = None, total_production_df: pd.DataFrame = None, pcoa_dataframe: pd.DataFrame = None, savepath: str = None):
    """Save every dataframe to save_path input.

    Args:
        all_data (dict): all_sample dataframe, metadata, main_dataframe and producer_dataframe
        norm_abundance_df (Dataframe): abundance dataframe normalised
        long_taxo_df (Dataframe): taxonomic dataframe in long format
        total_production_df (Dataframe): total production dataframe
        savepath (str): path to save all files.
    """
    filename = os.path.join(savepath,"postaviz_dataframes.h5")

    if savepath is None:
        print("save_path is None, data can't be saved into HDF5 format.")
        sys.exit(1)
    
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

    if os.path.isfile(filename):
        return filename
    
    else:
        print('Error when creating HDF5 dataframe storage.')
        sys.exit(1)

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

    print(results)
    results.reset_index(inplace=True)
    results = results.merge(metadata_dataframe,'inner','smplID')
    print(results)
    # results["smplID"] = results["smplID"].astype("category")

    return results


def stat_on_plot(data: dict, layer: int):
    """Apply Wilcoxon or Mann-Whitney test on each pair of a dataframe.

    Args:
        data (dict): Dictionnary of each pair to test.
        layer (int): Number of layer to the dict. 2 for nested dict.

    Returns:
        Dataframe: Dataframe of test's results.
    """
    start_timer = time.time()
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
                                except Exception as e:
                                    error_log.append([pair_1,pair1_data,pair_2,pair2_data,e])
                                    continue
                            else:
                                try:
                                    test_value, p_value = stats.mannwhitneyu(pair1_data, pair2_data)
                                    test_type = "Mann-Whitney"
                                except Exception as e:
                                    error_log.append([pair_1,pair1_data,pair_2,pair2_data,e])
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
    print(time.time() - start_timer)  
    return res


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
    sample_data = {}
    
    total_production_dataframe = None
    long_taxonomic_dataframe = None
    normalised_abundance_dataframe = None
    pcoa_dataframe = None

    if save_path is None or not is_valid_dir(save_path):
        print("save_path is none, ignoring load data function.")
        return sample_data, metadata, main_dataframe, normalised_abundance_dataframe, long_taxonomic_dataframe, producers_dataframe, total_production_dataframe, pcoa_dataframe
    
    for root, dirname, filename in os.walk(save_path):

        if "all_samples_dataframe_postaviz" in dirname:
            dirpath = os.path.join(root,"all_samples_dataframe_postaviz")
            for sample in os.listdir(dirpath):
                if sample.endswith("_cscope.tsv"):
                    try:
                        sample_name = sample.split("_", 1)[0]
                        sample_data[sample_name] = {}
                        sample_data[sample_name]["cscope"] = open_tsv(os.path.join(dirpath, sample),True,True)
                    except Exception as e:
                        print("Read tsv file: ",sample,"in", os.path.join(dirpath,sample)," failed.")
                        continue

        if "metadata_dataframe_postaviz.tsv" in filename:
            metadata = open_tsv(os.path.join(root, "metadata_dataframe_postaviz.tsv"))

        if "main_dataframe_postaviz.tsv" in filename:
            main_dataframe = open_tsv(os.path.join(root, "main_dataframe_postaviz.tsv"))

        if "producers_dataframe_postaviz.tsv" in filename:
            producers_dataframe = open_tsv(os.path.join(root, "producers_dataframe_postaviz.tsv"))

        if "total_production_dataframe_postaviz.tsv" in filename:
            total_production_dataframe = open_tsv(os.path.join(root, "total_production_dataframe_postaviz.tsv"))

        if "taxonomic_dataframe_postaviz.tsv" in filename:
            long_taxonomic_dataframe = open_tsv(os.path.join(root, "taxonomic_dataframe_postaviz.tsv"))

        if "normalised_abundance_dataframe_postaviz.tsv" in filename:
            normalised_abundance_dataframe = open_tsv(os.path.join(root, "normalised_abundance_dataframe_postaviz.tsv"))
        
        if "pcoa_dataframe_postaviz.tsv" in filename:
            pcoa_dataframe = open_tsv(os.path.join(root, "pcoa_dataframe_postaviz.tsv"))

    return sample_data, metadata, main_dataframe, normalised_abundance_dataframe, long_taxonomic_dataframe, producers_dataframe, total_production_dataframe, pcoa_dataframe

def save_all_dataframe(sample_data, metadata, main_dataframe, normalised_abundance_dataframe, long_taxonomic_dataframe, producers_dataframe, total_production_dataframe, pcoa_dataframe, savepath):
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

    data_dir = os.path.join(savepath,"all_samples_dataframe_postaviz")

    if os.path.isdir(data_dir):
        print(data_dir, "directory already exist in save_path.")
        pass
    else:
        os.makedirs(data_dir)
        print("No existing directory for path: ", data_dir)
 
        for sample in sample_data.keys():
            try:
                filename = str(sample+"_cscope.tsv")
                filename_path = os.path.join(data_dir,filename)
                sample_data[sample]["cscope"].to_csv(filename_path,sep="\t")
            except Exception as e:
                print("Saving ERROR",e)
    
    if os.path.isfile(os.path.join(savepath,"main_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"main_dataframe_postaviz.tsv"), " file already exist.")
    else:
        if main_dataframe is not None: main_dataframe.to_csv(os.path.join(savepath,"main_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(main_dataframe) else False)

    if os.path.isfile(os.path.join(savepath,"metadata_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"metadata_dataframe_postaviz.tsv"), " file already exist.")
    else:
        if metadata is not None: metadata.to_csv(os.path.join(savepath,"metadata_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(metadata) else False)
 
    if os.path.isfile(os.path.join(savepath,"producers_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"producers_dataframe_postaviz.tsv"), " file already exist.")
    else:
        if producers_dataframe is not None: producers_dataframe.to_csv(os.path.join(savepath,"producers_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(producers_dataframe) else False)

    if os.path.isfile(os.path.join(savepath,"normalised_abundance_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"normalised_abundance_dataframe_postaviz.tsv"), " file already exist")
    else:
        if normalised_abundance_dataframe is not None: normalised_abundance_dataframe.to_csv(os.path.join(savepath,"normalised_abundance_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(normalised_abundance_dataframe) else False)
 
    if os.path.isfile(os.path.join(savepath,"taxonomic_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"taxonomic_dataframe_postaviz.tsv"), " file already exist")
    else:
        if long_taxonomic_dataframe is not None: long_taxonomic_dataframe.to_csv(os.path.join(savepath,"taxonomic_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(long_taxonomic_dataframe) else False)
 
    if os.path.isfile(os.path.join(savepath,"total_production_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"total_production_dataframe_postaviz.tsv"), " file already exist")
    else:
        if total_production_dataframe is not None: total_production_dataframe.to_csv(os.path.join(savepath,"total_production_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(total_production_dataframe) else False)

    if os.path.isfile(os.path.join(savepath,"pcoa_dataframe_postaviz.tsv")):
        print(os.path.join(savepath,"pcoa_dataframe_postaviz.tsv"), " file already exist")
    else:
        if pcoa_dataframe is not None: pcoa_dataframe.to_csv(os.path.join(savepath,"pcoa_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(pcoa_dataframe) else False)

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