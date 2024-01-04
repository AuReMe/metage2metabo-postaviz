import csv
import json

import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci


def open_txt(filename, with_pandas: bool = False):
    if not with_pandas:
        with open(filename) as file:
            content = file.readlines()
            return content
    else:
        content = pd.read_csv(filename, sep="\t")
        return content


def write_tsv(filename: str, header: list, rows: list):
    """Write tsv metadata

    Args:
        filename (str): name of file
        header (list): list of columns name
        rows (list): nested list of rows value
    """
    with open(filename, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            pass


def extract_genus(str_container, full: bool = True):
    for i in range(len(str_container)):
        if str_container[i] == "\t":
            res = ""
            for y in str_container[i + 1 : -1]:
                res += y
            if full:
                return res
            else:
                res = res.split(";")
                return res[-1]


def decoding_sbml(cscope_file):
    c_file = open_json(cscope_file)["com_scope"]
    uncoded = []
    for coded in c_file:
        id, id_type, compart = cfci(coded)
        uncoded.append(tuple((id, id_type, compart)))
    return uncoded


def open_json(file_path):
    with open(file_path) as file:
        file_data = json.load(file)
    return file_data


def open_tsv(file_name):
    data = pd.read_csv(file_name, sep="\t")
    return data


def naive_search_db(uncoded_compound: list, db_path):
    hit = []
    db = open_tsv(db_path)
    for metabolite_tuple in uncoded_compound:
        metabolite = metabolite_tuple[0]
        for db_metabolite in db["compound_id"]:
            if metabolite in db_metabolite:
                hit.append(metabolite)
    print(len(hit), "match between", db.shape[0], "metabolites in db")
    return hit


uncoded_list = decoding_sbml("/home/lbrindel/inria_m2m/output_dir/metadata_ouput/E1/community_analysis/comm_scopes.json")
# print(uncoded_list[1])

# all_match = naive_search_db(uncoded_list, "/home/lbrindel/Downloads/compounds_26_5_level1.tsv")

# all_matchv2 = naive_search_db(uncoded_list, "/home/lbrindel/Downloads/compounds_26_5_level2.tsv")


sample_file = pd.read_csv("/home/lbrindel/Downloads/ERAS1d0.txt", header=None)
taxo_file = open_txt("/home/lbrindel/Downloads/mapping_mgs_genus.txt", True)

sample_file.columns = ["mgs"]

sample_file["test"] = [i for i in range(len(sample_file))]

merged_results = sample_file.merge(taxo_file,how='left')

print(sample_file)
print(merged_results)