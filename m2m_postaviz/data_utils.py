import json
import os
import os.path
import tarfile
import time
import numpy as np
from scipy import stats

import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci


def extract_tarfile(tar_file, outdir):
    file = tarfile.open(tar_file, "r:gz")

    file.extractall(outdir, filter="data")
    # if sys.version_info >= (3, 12):
    # else:
    #     tar.extractall(outdir)


def multiply_production_abundance(df_row: pd.Series, abundance_matrix: pd.DataFrame, sample_id):
    df_row = df_row.astype(float)
    abundance_value = abundance_matrix.at[df_row.name, sample_id]
    df_row *= abundance_value
    return df_row


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
            print(f'>run {i+1} took {time_duration} seconds')
        avg_duration = sum(results) / n_repeats
        print(f'Took {avg_duration} seconds on average for {func.__name__} function.')
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

    smpl_abundance = []
    smpl_index = []
    
    for sample in sample_data.keys():

        sample_matrix = sample_data[sample]["cscope"].copy()

        if not is_indexed_by_id(sample_matrix):
            sample_matrix.set_index("smplID", inplace=True)

        sample_matrix = sample_matrix.apply(lambda row: row.astype(float)*abundance_matrix.at[row.name,sample], axis=1)
        sample_matrix = sum_squash_table(sample_matrix, sample)

        smpl_abundance.append(sample_matrix)
        smpl_index.append(str(sample))

    for sample in sample_data.keys():

        abundance_matrix_normalised = abundance_matrix.apply(lambda x: x / x.sum(),axis=0)

        sample_matrix = sample_data[sample]["cscope"].copy()

        if not is_indexed_by_id(sample_matrix):
            sample_matrix.set_index("smplID", inplace=True)

        sample_matrix = sample_matrix.apply(lambda row: row.astype(float)*abundance_matrix_normalised.at[row.name,sample], axis=1)
        sample_matrix = sum_squash_table(sample_matrix, sample)

        smpl_norm_abundance.append(sample_matrix)
        smpl_norm_index.append(str(sample))

    normalized_abundance = pd.concat(smpl_norm_abundance)
    normalized_abundance.fillna(0, inplace=True)
    normalized_abundance.insert(0, "smplID", smpl_norm_index)
    normalized_abundance.set_index("smplID", inplace=True, drop=True)

    vanilla_abundance = pd.concat(smpl_abundance)
    vanilla_abundance.fillna(0, inplace=True)
    vanilla_abundance.insert(0, "smplID", smpl_index)
    vanilla_abundance.set_index("smplID", inplace=True, drop=True)

    return normalized_abundance, vanilla_abundance


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
    new_dataframe = {}
    # Flip avec métabolites en index et bin en column
    abundance_table = abundance_table.T
    for index, row in abundance_table.iterrows():
        # Pour chaque métabolites, calcul le total crée par l'ensemble des bin de l'échantillon
        new_dataframe[index] = row.values.sum()
    return pd.DataFrame(new_dataframe, index=[sample_id])


def taxonomic_overview(list_bin_id, taxonomic_dataframe, metadata, mode: str = "cscope", with_count: bool = True):
    current_selection = []
    current_rank = taxonomic_dataframe.columns[-1]

    for sample in list_bin_id.keys():
        tmp_df = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(list_bin_id[sample])]
        tmp_df = tmp_df[["mgs", current_rank]]
        if with_count:
            tmp_df = tmp_df.groupby([current_rank]).count()
            tmp_df = tmp_df.reset_index()
            tmp_df.columns = [current_rank, "Count"]
        value, label = get_metadata(sample, metadata)
        for i in range(len(label)):
            tmp_df.insert(0, label[i], value[i])
        tmp_df.insert(0, "nb_taxon", tmp_df[current_rank].nunique())

        current_selection.append(tmp_df)

    return pd.concat(current_selection, join="outer", ignore_index=True)


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


def get_added_value_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["advalue"])


def get_individual_production_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["iscope"].columns)


def get_total_production_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["cscope"].columns)


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


def get_metabolic_info(sample_data: dict, metadataframe: pd.DataFrame, taxonomic_dataframe: pd.DataFrame):
    """Really ugly and probably useless need rework.
    Return a dataframe with global information of all the sample in input.

    Args:
        sample_data (dict): _description_
        metadataframe (pd.DataFrame): _description_
        taxonomic_dataframe (pd.DataFrame): _description_

    Returns:
        _type_: _description_
    """
    tot_size = []
    ind_size = []
    ad_size = []
    taxo_size = []
    model_size = []
    for sample in metadataframe["smplID"]:
        tot_size.append(get_total_production_size(sample, sample_data))
        ind_size.append(get_individual_production_size(sample, sample_data))
        ad_size.append(get_added_value_size(sample, sample_data))
        taxo_size.append(get_taxonomy_size(sample_data[sample]["iscope"], taxonomic_dataframe))
        model_size.append(get_taxonomy_size(sample_data[sample]["iscope"], taxonomic_dataframe, only_metabolic_model_size=True))
    new_df = metadataframe
    new_df["prod_community"] = tot_size
    new_df["prod_individual"] = ind_size
    new_df["added_value_total"] = ad_size
    new_df["Numbers of models"] = model_size
    new_df["Numbers of species"] = taxo_size
    return new_df


def get_metabolic_model_prod_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["cscope"].columns)


def add_row(df: pd.DataFrame, row: list):
    df.loc[len(df)] = row


def deal_with_duplicated_row(df: pd.DataFrame, dropcol=None):
    if not dropcol == None:
        df = df.drop(dropcol, axis=1)
    df = pd.get_dummies(df, prefix="", prefix_sep="", columns=["category"])
    df = df.groupby(["compound_id"]).sum()
    return df


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


def get_percent_value(df: pd.DataFrame, transpose: bool = False) -> pd.DataFrame:
    """Calculate the total of each columns and replace value by their percentage.

    Args:
        df (pd.DataFrame): Count matrix as dataframe value must be row indexed and sample in columns.
        transpose (bool, optional): Transpose the matrix. Defaults to False.

    Returns:
        pd.DataFrame: ERROR 404
    """
    if transpose:
        df = df.T
    total = df.sum()
    percentage = (df / total) * 100
    return percentage


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


def open_tsv(file_name: str, rename_columns: bool = False, first_col: str = "smplID"):
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
        data.columns.values[1:] = sbml_to_classic(data.columns.values[1:])
    return data


def get_scopes(file_name, path):
    for root, dirs, files in os.walk(path):
        if file_name in files:
            iscope_dataframe = open_tsv(os.path.join(root, file_name), rename_columns=True)
            return iscope_dataframe


def convert_to_dict(list_of_produced_cpd):
    decoded_cpd_list = sbml_to_classic(list_of_produced_cpd)
    data = {}
    for i in decoded_cpd_list:
        value = [1]
        data[i] = value
    return data


def is_indexed_by_id(df: pd.DataFrame):
    if df.index.name == "smplID":
        return True
    else:
        return False


def get_column_size(df: pd.DataFrame):
    size = len(df.columns)
    return size


def get_row_size(df: pd.DataFrame):
    size = len(df)
    return size


def get_columns_index(df: pd.DataFrame, key_list):
    index_list = [0]
    for k in key_list:
        index_list.append(df.columns.get_loc(k))
    return index_list


def convert_to_dataframe(all_file):
    all_dataframe = []
    for file in all_file:
        current_file = open_json(file)["com_scope"]
        all_dataframe.append(pd.DataFrame.from_dict(convert_to_dict(current_file)))
    return all_dataframe


def sbml_to_classic(list_of_metabolites):
    uncoded = []
    for coded in list_of_metabolites:
        id, id_type, compart = cfci(coded)
        uncoded.append(id)
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


def retrieve_all_sample_data(path):
    """Retrieve iscope, cscope, added_value and contribution_of_microbes files in the path given using os.listdir().

    Args:
        path (str): Directory path given in CLI.

    Returns:
        dict: Return a nested dict object where each key is a dictionnary of a sample. The key of those second layer dict [iscope, cscope, advalue, contribution] give acces to these files.
    """
    data_by_sample = {}
    for sample in os.listdir(path):
        if os.path.isdir(os.path.join(path, sample)):
            data_by_sample[sample] = {}
            # data_by_sample[sample]["iscope"] = get_scopes("rev_iscope.tsv", os.path.join(path, sample))
            data_by_sample[sample]["cscope"] = get_scopes("rev_cscope.tsv", os.path.join(path, sample))
            # data_by_sample[sample]["advalue"] = open_added_value("addedvalue.json", os.path.join(path, sample))
            # data_by_sample[sample]["contribution"] = get_contributions("contributions_of_microbes.json", os.path.join(path, sample))
    return data_by_sample


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


def build_df(dir_path, metadata, abundance_path):
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
    data_by_sample = {}
    sample_data = retrieve_all_sample_data(dir_path)

    ### Multi-Threading for building data, slower than vanilla (for NOW).

    # all_task = []
    # data_by_sample = {}
    # for sample_directory in os.listdir(dir_path):
    #     if os.path.isdir(os.path.join(dir_path, sample_directory)):
    #         task = threading.Thread(target=retrieve_sample_data,args=(os.path.join(dir_path, sample_directory), sample_directory, data_by_sample))
    #         task.start()
    #         all_task.append(task)

    # for task in all_task:
    #     task.join()

    ### Main dataframe and metadata dataframe. ###

    global_data = {}
    dir_files = get_files("comm_scopes.json", dir_path)
    file_list = []
    dir_list = []
    for file, dir in dir_files:
        file_list.append(file)
        dir_list.append(dir)

    all_df = convert_to_dataframe(file_list)

    main_df = pd.concat(all_df, join="outer", ignore_index=True)
    main_df = main_df.fillna(0)
    main_df = main_df.astype(int)
    main_df.insert(0, "smplID", dir_list)

    global_data["metadata"] = open_tsv(metadata)

    global_data["main_dataframe"] = main_df

    abundance_file = open_tsv(abundance_path)

    norm_abundance_data, abundance_data = relative_abundance_calc(abundance_file, sample_data)

    return global_data, sample_data, norm_abundance_data, abundance_data


def build_test_data(test_dir_path):
    """Load, open and return the test data to be used by the DataStorage object.

    Args:
        test_dir_path (_type_): Path to the archive

    Returns:
        dict, dict, pd.dataframe: 2 dict, 1 pandas dataframe
    """
    global_data = {}
    sample_data = {}
    for file in os.listdir(test_dir_path):
        filename = os.fsdecode(file)
        if not filename.endswith(".tsv"):
            continue
        if not filename.startswith("sample"):
            if filename == "main_table.tsv":
                global_data["main_dataframe"] = open_tsv(os.path.join(test_dir_path, filename))
            if filename == "metadata_table.tsv":
                global_data["metadata"] = open_tsv(os.path.join(test_dir_path, filename))
            if filename == "abundance_table.tsv":
                abundance_table = open_tsv(os.path.join(test_dir_path, filename))
        else:
            current_file = os.path.splitext(filename)[0]
            current_file = current_file.split("_")
            sample_data[current_file[2]] = {}

    for sample_file in os.listdir(test_dir_path):
        filename = os.fsdecode(sample_file)
        if filename.startswith("sample"):
            current_file = os.path.splitext(filename)[0]
            current_file = current_file.split("_")
            sample_data[current_file[2]][current_file[1]] = open_tsv(os.path.join(test_dir_path, filename))

    return global_data, sample_data, abundance_table


def produce_test_data(global_data, sample_data, abundance_data):
    global_data["main_dataframe"].to_csv("/home/lbrindel/output/postaviz/data_test/main_table.tsv", sep="\t", index=False)
    global_data["metadata"].to_csv("/home/lbrindel/output/postaviz/data_test/metadata_table.tsv", sep="\t", index=False)
    abundance_data.to_csv("/home/lbrindel/output/postaviz/data_test/abundance_table.tsv", sep="\t", index=False)
    for sample in sample_data.keys():
        sample_data[sample]["cscope"].to_csv(
            "/home/lbrindel/output/postaviz/data_test/sample_cscope_" + sample + ".tsv", sep="\t", index=False
        )
        sample_data[sample]["iscope"].to_csv(
            "/home/lbrindel/output/postaviz/data_test/sample_iscope_" + sample + ".tsv", sep="\t", index=False
        )
    quit()


def unit_test_abundance():
    mock_cscope1 = pd.DataFrame(
        data=[
            [1, 0, 0, 0, 1, 1, 0],
            [0, 0, 1, 0, 1, 1, 0],
            [1, 0, 0, 1, 1, 1, 1],
            [0, 1, 0, 0, 0, 0, 0],
        ],              
        index=["bin1", "bin3", "spec437", "bin1030"],
        columns=["CPD1", "CPD2", "CPD3", "CPD4", "CPD5", "CPD6", "CPD7"]
    )
    mock_cscope1.index.name = "smplID"
    mock_cscope2 = pd.DataFrame(
        data=[
            [1, 0, 0, 0, 1, 1, 0],
            [0, 0, 1, 0, 1, 1, 0],
            [1, 0, 0, 1, 1, 1, 1],
            [0, 1, 0, 0, 0, 0, 0],
        ],
        index=["bin1", "spec88", "bin1030", "bin502"],
        columns=["CPD1", "CPD2", "CPD3", "CPD4", "CPD5", "CPD6", "CPD7"]
    )
    mock_cscope2.index.name = "smplID"

    mock_abundance_df = pd.DataFrame(
        data=[
            [15, 25, 2, 5],
            [1, 30, 12, 2],
            [8, 0, 0, 1],
            [2, 1, 18, 0],
            [8, 2, 2, 5],
            [0, 10, 12, 2],
            [7, 2, 6, 1],
            [0, 1, 18, 0],
        ],
        index=["bin1", "spec88", "bin1030", "bin502", "bin3", "spec437","bin999","spec999"],
        columns=["mock1", "mock2", "mock3", "mock4"]
    )

    sample_mock = dict()
    sample_mock["mock1"] = {}
    sample_mock["mock1"]["cscope"] = mock_cscope1
    sample_mock["mock2"] = {}
    sample_mock["mock2"]["cscope"] = mock_cscope2

    normalised_mock_ab = mock_abundance_df.apply(lambda x: x / x.sum(), axis=0)
    expected_results = sample_mock["mock1"]["cscope"].T.dot(normalised_mock_ab.loc[normalised_mock_ab.index.isin(sample_mock["mock1"]["cscope"].index)]["mock1"])
    expected_results2 = sample_mock["mock2"]["cscope"].T.dot(normalised_mock_ab.loc[normalised_mock_ab.index.isin(sample_mock["mock2"]["cscope"].index)]["mock2"])

    expected_results.rename("mock1", inplace=True)
    expected_results2.rename("mock2", inplace=True)
    expected_df = pd.DataFrame([expected_results,expected_results2])
    # print(expected_df)
    
    observed_results = relative_abundance_calc(mock_abundance_df,sample_mock)
    # print(observed_results[0])
    # print("--------------")
    # print(observed_results[1])

    assert expected_df.equals(observed_results[0]), "Expected abundance dataframe from unit_test_abundance() and abundance dataframe from tested function are not equals."
    return True



def add_p_value_annotation(fig, array_columns, subplot=None, _format=dict(interline=0.07, text_height=1.07, color='black')):
    ''' Adds notations giving the p-value between two box plot data (t-test two-sided comparison)
    
    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare 
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    '''
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        y_range[i] = [1.01+i*_format['interline'], 1.02+i*_format['interline']]

    # Get values from figure
    fig_dict = fig.to_dict()
    # Get indices if working with subplots
    if subplot:
        if subplot == 1:
            subplot_str = ''
        else:
            subplot_str =str(subplot)
            print('Subplot str is : ',subplot_str)
        indices = [] #Change the box index to the indices of the data for that subplot
        for index, data in enumerate(fig_dict['data']):
            if data['xaxis'] == 'x' + subplot_str:
                indices = np.append(indices, index)
        indices = [int(i) for i in indices]
        print((indices))
    else:
        subplot_str = ''

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        print(column_pair)
        if subplot:
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair

        # Mare sure it is selecting the data and subplot you want
        # print(data_pair)
        # print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['xaxis'])
        # print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['xaxis'])

        # Get the p-value
        pvalue = stats.ttest_ind(
            fig_dict['data'][data_pair[0]]['y'],
            fig_dict['data'][data_pair[1]]['y'],
            equal_var=False,
        )[1]
        if pvalue >= 0.05:
            symbol = 'ns'
        elif pvalue >= 0.01: 
            symbol = '*'
        elif pvalue >= 0.001:
            symbol = '**'
        else:
            symbol = '***'
        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][0], 
            x1=column_pair[0], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        # Horizontal line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][1], 
            x1=column_pair[1], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[1], y0=y_range[index][0], 
            x1=column_pair[1], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(dict(font=dict(color=_format['color'],size=14),
            x=(column_pair[0] + column_pair[1])/2,
            y=y_range[index][1]*_format['text_height'],
            showarrow=False,
            text=symbol,
            textangle=0,
            xref="x"+subplot_str,
            yref="y"+subplot_str+" domain"
        ))
    return fig