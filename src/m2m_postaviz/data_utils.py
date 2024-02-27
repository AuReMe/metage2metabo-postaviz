import json
import os
import os.path
from m2m_postaviz.time_decorator import timeit
import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci
import cProfile
import time
import threading
from scipy import stats

def taxonomic_overview(list_bin_id, taxonomic_dataframe, metadata, mode: str = "cscope"):
    current_selection = []
    current_rank = taxonomic_dataframe.columns[-1]
    print(type(current_rank))

    for sample in list_bin_id.keys():
        tmp_df = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(list_bin_id[sample])]
        tmp_df = tmp_df[['mgs',current_rank]]
        tmp_df = tmp_df.groupby([current_rank]).count()
        tmp_df = tmp_df.reset_index()
        tmp_df.columns = [current_rank,'Count']
        value, label = get_metadata(sample, metadata)
        for i in range(len(label)):
            tmp_df.insert(0,label[i],value[i])

        current_selection.append(tmp_df)
        
    return pd.concat(current_selection, join='outer', ignore_index=True)


def get_metadata(sample_id: str, metadata: pd.DataFrame):
    return tuple([metadata.loc[metadata["Name"] == sample_id].values.tolist()[0],metadata.loc[metadata["Name"] == sample_id].columns.to_list()])


def taxonomy_groupby(metadata: pd.DataFrame, current_sample: str, bin_id_by_sample: dict, taxonomic_dataframe: pd.DataFrame, target_rank:str = "Genus", taxonomic_choice: list = []):
    """Generate a taxonomic count dataframe from a sample id, his dataframe and the rank choosen. Return only the selected taxonomic choice. 

    Args:
        current_sample (str): Sample's id
        taxonomic_dataframe (pd.DataFrame): The taxonomic dataframe
        target_rank (str, optional): Selected rank in shiny's input. Defaults to "Genus".
        taxonomic_choice (list, optional): The list of taxonomic selection from shiny's input to keep in the returned Dataframe. Defaults to [].

    Returns:
        Dataframe: Pandas dataframe containing the selected taxonomic choice and rank.
    """
    taxonomic_choice = list(taxonomic_choice)
    if len(taxonomic_choice) == 0:
        print("The taxonomic choice list is empty")
        return
    
    df = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(bin_id_by_sample[current_sample])]
    for choice in taxonomic_choice:
        if not choice in df[target_rank].unique():
            taxonomic_choice.remove(choice)
            print(choice, " Removed.")
    
    df = df[["mgs",target_rank]]
    df = df.groupby([target_rank]).count()
    df = df.reset_index()
    df.columns = [target_rank,"Count"]
    df = df.loc[df[target_rank].isin(taxonomic_choice)]
    value,label = get_metadata(current_sample, metadata)
    for i in range(len(label)):
        df.insert(0, label[i], value[i])
    return df


def taxonomic_dataframe_from_input(taxonomic_rank: str, bin_id_by_sample: dict, taxonomic_choice: list, taxonomic_dataframe: pd.DataFrame, metadata: pd.DataFrame):
    results = []
    for sample in bin_id_by_sample.keys():
        results.append(taxonomy_groupby(metadata, sample, bin_id_by_sample, taxonomic_dataframe, taxonomic_rank, taxonomic_choice))
    
    final_results = pd.concat(results, join="outer", ignore_index=True)
    return final_results


def wilcoxon_test(group1,group2):
    try:
         res = stats.wilcoxon(group1,group2)
    except Exception as e:
        res = e
    return res


def student_test(group1,group2):
    try:
         res = stats.ttest_ind(group1,group2)
    except Exception as e:
        res = e
    return res


def get_added_value_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["advalue"])


def get_individual_production_size(sample:str, metadataframe: dict):
    return len(metadataframe[sample]["iscope"].columns)


def get_total_production_size(sample:str, metadataframe: dict):
    return len(metadataframe[sample]["cscope"].columns)


def get_taxonomy_size(sample_data: pd.DataFrame, taxonomic_dataframe: pd.DataFrame, only_metabolic_model_size: bool = False):
    if only_metabolic_model_size:
       taxonomy_size = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(sample_data["Name"])]
       return len(taxonomy_size) 
    taxonomy_size = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(sample_data["Name"])]
    taxonomy_size = taxonomy_size[["mgs","Genus"]]
    taxonomy_size = taxonomy_size.groupby(["Genus"]).count()
    taxonomy_size = taxonomy_size.reset_index()
    return len(taxonomy_size)


def get_metabolic_info(sample_data: dict, metadataframe: pd.DataFrame, taxonomic_dataframe: pd.DataFrame):
    tot_size = []
    ind_size = []
    ad_size = []
    taxo_size = []
    model_size = []
    for sample in metadataframe["Name"]:
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


def remove_metadata(df: pd.DataFrame) -> pd.DataFrame:
    new_df = df.drop("Test", axis=1)
    new_df = new_df.drop("Days", axis=1)
    # df = df.drop("Name", axis=1)
    new_df = new_df.drop("Antibiotics", axis=1)
    new_df = new_df.set_index("Name")
    return new_df


def add_row(df: pd.DataFrame, row: list):
    df.loc[len(df)] = row


def deal_with_duplicated_row(df: pd.DataFrame, dropcol = None):
    if not dropcol == None:
        df = df.drop(dropcol,axis=1)
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


def get_percent_value(df: pd.DataFrame, transpose: bool = False)  -> pd.DataFrame:
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


def get_files(file_name, path, with_directory_name : bool = True):
    result = []
    for root, dirs, files in os.walk(path):
        if file_name in files:
            if with_directory_name:
                result.append([os.path.join(root, file_name), os.path.basename(os.path.dirname(os.path.dirname(os.path.join(root, file_name))))])
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


def open_tsv(file_name, rename_columns : bool = False, first_col : str = "Name"):
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


def convert_to_dict(file_as_list):
    decoded_list = sbml_to_classic(file_as_list)
    data = {}
    for i in decoded_list:
        value = [1]
        data[i] = value
    return data


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
    for root, dirs , files in os.walk(path):
        if file_name in files:
            contributions_file = open_json(os.path.join(root, file_name))
            contributions_file = contribution_processing(contributions_file)
            return contributions_file


def retrieve_all_sample_data(path):
    data_by_sample = {}
    for sample in os.listdir(path):
        if os.path.isdir(os.path.join(path, sample)):
            data_by_sample[sample] = {}
            data_by_sample[sample]["iscope"] = get_scopes("rev_iscope.tsv", os.path.join(path, sample))
            data_by_sample[sample]["cscope"] = get_scopes("rev_cscope.tsv", os.path.join(path, sample))
            data_by_sample[sample]["advalue"] = open_added_value("addedvalue.json", os.path.join(path, sample))
            data_by_sample[sample]["contribution"] = get_contributions("contributions_of_microbes.json", os.path.join(path, sample))
    return data_by_sample


def retrieve_sample_data(sample_directory_path, sample_name, sample_dictionnary: dict):
    sample_dictionnary[sample_name] = {}
    sample_dictionnary[sample_name]["iscope"] = get_scopes("rev_iscope.tsv", sample_directory_path)
    sample_dictionnary[sample_name]["cscope"] = get_scopes("rev_cscope.tsv", sample_directory_path)
    sample_dictionnary[sample_name]["advalue"] = open_added_value("addedvalue.json", sample_directory_path)
    sample_dictionnary[sample_name]["contribution"] = get_contributions("contributions_of_microbes.json", sample_directory_path)
    return 

def merge_metadata_with_df(main_dataframe, metadata):
    return pd.merge(metadata, main_dataframe, how="left")


def merge_df(left_df, right_df, how : str = "left"):
    data = left_df.iloc[:,0]
    filter = right_df.loc[right_df["mgs"].isin(data)]
    return filter


def multiply_production_abundance(row: pd.Series, abundance_matrix: pd.DataFrame,sample_id):
    row = row.astype(float)
    abundance_value = abundance_matrix.at[row.name,sample_id]
    row *= abundance_value
    return row


def generate_stoichiometric_matrix(binary_matrix: pd.DataFrame, abundance_matrix: pd.DataFrame, sample_id: str):
    """Produce a stoichiometric matrix from a binary matrix (presence or absence) and an abundance matrix.

    Args:
        binary_matrix (pandas.DataFrame): Binary dataframe containing the presence or absence of compounds.
        abundance_matrix (pandas.DataFrame): Dataframe containing the abundance of each bin in sample.
        sample_id (str): Name of the sample.

    Returns:
        pandas Dataframe: Dataframe with the theorical quantity of compounds produced by each bin. 
    """
    binary_matrix.set_index("Name",inplace=True)
    binary_matrix = binary_matrix.apply(lambda row: multiply_production_abundance(row, abundance_matrix,sample_id),axis=1) 
    return binary_matrix


# @timeit(repeat=3,number=10)
def build_df(dir_path, metadata):
    """
    Extract community scopes present in directory then build the a single dataframe from the metabolites produced by each comm_scopes.

    Args:
        dir_path (str): Directory path containing comm scopes
        metadata (tsv file): tsv file containing the metadata of the scopes. The number of row must be equal to the number of comm_scopes given in dir_path.

    Returns:
        pandas_DataFrame:
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
    main_df.insert(0, "Name", dir_list)

    global_data["metadata"] = open_tsv(metadata)

    global_data["main_dataframe"] = main_df
    
    return global_data, sample_data

def performance_test(dir,meta):
    start_time = time.perf_counter()

    # pr = cProfile.Profile()
    # pr.enable()

    build_df(dir,meta)


    end_time = time.perf_counter()

    Total_time = end_time - start_time
    print("TEST TIME: ", Total_time)
    # pr.disable()
    # pr.print_stats(sort='time')

    # cProfile.runctx('build_df(dir,meta)', globals(), locals(), sort=1)