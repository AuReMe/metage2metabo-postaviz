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


def get_scopes_dirname(file_name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if file_name in files:
            result.append(
                [os.path.join(root, file_name), os.path.basename(os.path.dirname(os.path.dirname(os.path.join(root, file_name))))]
            )
    return result


def json_to_df(list_of_filepath):
    foi = []
    for file in list_of_filepath:
        result = open_json(file)["com_scope"]
        foi.append(result)
    return foi


def open_json(file_path):
    with open(file_path) as file:
        file_data = json.load(file)
    return file_data


def open_tsv(file_name):
    data = pd.read_csv(file_name, sep="\t")
    return data


def get_all_iscope(file_name, path):
    all_iscopes = {}
    for root, dirs, files in os.walk(path):
        if file_name in files:
            dir_name = os.path.basename(os.path.dirname(os.path.dirname(os.path.join(root, file_name))))
            iscope_df = open_tsv(os.path.join(root, file_name))
            column_name_series = {}
            for col in iscope_df.columns.values:
                id, id_type, compart = cfci(col)
                column_name_series[col] = id
            iscope_df.rename(columns=column_name_series,inplace=True)
            all_iscopes[dir_name] = iscope_df
    return all_iscopes


def convert_to_dict(file_as_list):
    data = {}
    for i in file_as_list:
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
    if not len(all_file) < 1:
        for file in all_file:
            all_dataframe.append(pd.DataFrame.from_dict(convert_to_dict(file)))
    return all_dataframe


def sbml_to_classic_all(list_of_cscope_path: list):
    uncoded_files = []
    for file in list_of_cscope_path:
        uncoded_files.append(sbml_to_classic(file))
    return uncoded_files


def sbml_to_classic(cscope_file):
    c_file = open_json(cscope_file)["com_scope"]
    uncoded = []
    for coded in c_file:
        id, id_type, compart = cfci(coded)
        uncoded.append(id)
    return uncoded


def merge_df(left_df, right_df, how : str = "left"):
    # print("-----------")
    # print("left_df:" ,left_df)
    # print("-----------")
    # print("rigth_df:" ,rigth_df)
    data = left_df.iloc[:,0]
    total = left_df.set_index("Unnamed: 0")
    print("TESTOU\n", total)
    total = left_df.sum()
    filter = right_df.loc[right_df["mgs"].isin(data)]
    print("TESTOU\n", total)
    
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
    iscopes_files = get_all_iscope("rev_iscope.tsv", dir_path)
    dir_files = get_scopes_dirname("comm_scopes.json", dir_path)
    file_list = []
    dir_list = []
    for file in dir_files:
        file_list.append(file[0])
        dir_list.append(file[1])

    all_uncoded_cscope = sbml_to_classic_all(file_list)
    all_df = convert_to_dataframe(all_uncoded_cscope)

    main_df = pd.concat(all_df, join="outer", ignore_index=True)
    main_df = main_df.fillna(0)
    main_df = main_df.astype(int)
    main_df.insert(0, "Name", dir_list)

    metadata = open_tsv(metadata)

    main_df = pd.merge(metadata, main_df, how="left")

    return main_df , iscopes_files
