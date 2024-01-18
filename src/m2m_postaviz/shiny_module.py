from shiny import ui
import os.path
import pandas as pd


def generate_input_select(id: str,label: str,choice: list):
    return ui.input_select(id=id,label=label,choices=choice)


def generate_card(id: list, label: list,choice_list: list):
    all_input = []
    for i in range(len(id)):
        all_input.append(generate_input_select(id[i],label[i],choice_list[i]))
    card = ui.card(
                ui.row(
                        inp for inp in all_input
                    )
                )
    return card


def read_rev_iscope(file_name, dirpath, list_name: str):
    result = {}
    for root, dirs, files in os.walk(dirpath):
        if file_name in files:
            current_directory = os.path.basename(os.path.dirname(os.path.dirname(os.path.join(root, file_name))))
            current_file = pd.read_csv(os.path.join(root, file_name), sep="\t")
            result[current_directory] = list(current_file[list_name])
    return result


def read_padmet_pathway(file_name, dirpath, pathway_namelist: str):
    result = {}
    for root, dirs, files in os.walk(dirpath):
        if file_name in files:
            current_directory = os.path.basename(os.path.dirname(os.path.join(root, file_name))).rstrip(".padmet_report")
            current_file = pd.read_csv(os.path.join(root, file_name), sep="\t")
            result[current_directory] = list(current_file[pathway_namelist])
    return result


def aggregate_and_count(pathway_list: list, database: pd.DataFrame):
    tmp_df = database.loc[database["pathway_id"].isin(pathway_list)]
    tmp_df = tmp_df.groupby(['category']).agg({'category': ['count']})
    tmp_df.reset_index(inplace=True)
    tmp_df.columns = ['category', 'count']
    return tmp_df


def pathway_data_processing():
    """Open both references tsv file and all pathway tsv file from the padmet report_network.
    WARNING PERFORMANCE WORK IN PROGRESS. 

    Returns:
        list: A list containing two dict, one for pathway lvl1 and the other for lvl2, each key is a bin name and each value is a pandas dataframe.
    """
    pathway_data = read_padmet_pathway("all_pathways.tsv", "/home/lbrindel/all_padmet_report/", "dbRef_id")
    database1 = pd.read_csv("/home/lbrindel/Downloads/pathways_26_5_level1.tsv", sep="\t")
    database2 = pd.read_csv("/home/lbrindel/Downloads/pathways_26_5_level2.tsv", sep="\t")

    category_count_lvl1 = {}
    category_count_lvl2 = {}
    for key in pathway_data.keys():
        category_count_lvl1[key] = aggregate_and_count(pathway_data[key], database1)
        category_count_lvl2[key] = aggregate_and_count(pathway_data[key], database2)
        
    all_iscope_category_count = [category_count_lvl1, category_count_lvl2]
    return all_iscope_category_count


def taxonomic_processing(iscope_bins: list, taxonomic_dataframe: pd.DataFrame):
    iscope_dataframe = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(iscope_bins)]
    all_taxo_data = {}
    for rank in iscope_dataframe.columns.values[1:]:
        tmp_df = iscope_dataframe.groupby([rank]).agg({rank: 'count'})
        all_taxo_data[rank] = tmp_df
    return all_taxo_data