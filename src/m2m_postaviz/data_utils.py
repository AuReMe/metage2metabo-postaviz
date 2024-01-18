import json
import os
import os.path

import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci

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


def get_added_value(file_name, directory_path):
    addedvalue = {}
    addedvalue_files = get_files(file_name, directory_path)
    for content, sample in addedvalue_files:
        addedvalue[sample] = sbml_to_classic(open_json(content)["addedvalue"])
    return addedvalue


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


def get_scopes(files_name, path):
    all_scopes = {}
    data = get_files(files_name, path, with_directory_name=True)
    for files, dir in data:
        all_scopes[dir] = open_tsv(files, rename_columns=True)
    return all_scopes


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


def merge_metadata_with_df(main_dataframe, metadata):
    return pd.merge(metadata, main_dataframe, how="left")


def merge_df(left_df, right_df, how : str = "left"):
    data = left_df.iloc[:,0]
    filter = right_df.loc[right_df["mgs"].isin(data)]
    return filter


def build_df(dir_path, metadata):
    """
    Extract community scopes present in directory then build the a single dataframe from the metabolites produced by each comm_scopes.

    Args:
        dir_path (str): Directory path containing comm scopes
        metadata (tsv file): tsv file containing the metadata of the scopes. The number of row must be equal to the number of comm_scopes given in dir_path.

    Returns:
        pandas_DataFrame:
    """

    all_data = {}

    all_data["cscope"] = get_scopes("rev_cscope.tsv", dir_path)
    all_data["iscope"] = get_scopes("rev_iscope.tsv", dir_path)
    all_data["added_value"] = get_added_value("addedvalue.json", dir_path)
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

    metadata = open_tsv(metadata)

    all_data["main_dataframe"] = main_df

    return all_data
