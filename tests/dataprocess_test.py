import pandas as pd
import os

from m2m_postaviz import data_utils as du

def test_data_processing():

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    MODULE_DIR = os.path.dirname(BASE_DIR)
    SOURCE_DIR = os.path.join(MODULE_DIR, "m2m_postaviz")
    TEST_DIR = os.path.join(SOURCE_DIR, "postaviz_test_data")
    if not os.path.isdir(os.path.join(TEST_DIR, "palleja/")):
        print("Palleja test directory does not exist.")
        file_to_extract = os.path.join(TEST_DIR, "table_test_postaviz.tar.gz")
        print("No Palleja/ directory found. \nExtract test data tarfile...")
        print("tar file path: ", file_to_extract)
        du.extract_tarfile(file_to_extract, TEST_DIR)
        TEST_DATA_CONTAINER = os.path.join(TEST_DIR, "palleja/")

    else:
        print("Palleja test directory exist.")
        TEST_DATA_CONTAINER = os.path.join(TEST_DIR, "palleja/")

    metadata_file = os.path.join(TEST_DATA_CONTAINER, "metadata_test_data.tsv")
    taxonomy_file = os.path.join(TEST_DATA_CONTAINER, "taxonomy_test_data.tsv")
    abundance_file = os.path.join(TEST_DATA_CONTAINER, "abundance_test_data.tsv")
    data_dictionnary, norm_abundance_data, taxonomic_data, production_data = du.build_df(TEST_DATA_CONTAINER, metadata_file, abundance_file, taxonomy_file)

    # Check if resulting all_data dictionnary containt all the keys.
    keys = ["main_dataframe", "metadata", "sample_data", "producers_long_format"]
    assert all(k in data_dictionnary for k in keys), "Missing keys in all_data dictionnary."

    # Check if all sample_data are dataframe and if the numbers of metabolites produced is the same by checking len(columns in df) and len(metabolites produced list in json)
    for sample in data_dictionnary["sample_data"].keys():
        assert isinstance(data_dictionnary["sample_data"][sample]["cscope"], pd.DataFrame), f"sample {sample} in sample_data is not a pandas dataframe object."
        sample_json = du.open_json(os.path.join(TEST_DATA_CONTAINER, sample+"/"+"community_analysis/"+"comm_scopes.json"))["com_scope"]
        assert len(data_dictionnary["sample_data"][sample]["cscope"].columns) == len(sample_json),"Length of sample columns and length of json metabolites production list are not the same."
        
        # .loc by sample index the main dataframe into a serie. then extract the values and sum() the get the of cpd amount produced. Then compare with json file.
        assert data_dictionnary["main_dataframe"].loc[data_dictionnary["main_dataframe"].index == sample].values.sum() == len(sample_json),f"The amount of metabolites produced for {sample} in main_dataframe is not the same as com_scopes json file."


