import json
import os
import os.path

import pandas as pd


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


def get_columns_index(df: pd.DataFrame,key_list):
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


def build_df(dir_path, metadata):
    """
    Extract community scopes present in directory then build the a single dataframe from the metabolites produced by each comm_scopes.

    Args:
        dir_path (str): Directory path containing comm scopes
        metadata (tsv file): tsv file containing the metadata of the scopes. The number of row must be equal to the number of comm_scopes given in dir_path.

    Returns:
        pandas_DataFrame:
    """
    dir_files = get_scopes_dirname("comm_scopes.json", dir_path)
    file_list = []
    dir_list = []
    for file in dir_files:
        file_list.append(file[0])
        dir_list.append(file[1])

    all_file = json_to_df(file_list)
    print("Community scopes extracted", len(all_file))
    all_df = convert_to_dataframe(all_file)

    main_df = pd.concat(all_df, join="outer", ignore_index=True)
    main_df = main_df.fillna(0)
    main_df = main_df.astype(int)
    main_df.insert(0, "Name", dir_list)
    # main_df.set_index("Name", inplace=True)

    metadata = open_tsv(metadata)
    print("Metadata tsv file extracted", len(metadata.columns))
    # metadata.set_index("Name", inplace=True)

    main_df = pd.merge(metadata, main_df, how="left")

    return main_df
