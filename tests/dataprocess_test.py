import pandas as pd
import os
import warnings

from m2m_postaviz import data_utils as du
from m2m_postaviz.data_struct import DataStorage

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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdf5_file_path, taxonomy_provided, abundance_provided = du.build_df(dir_path=TEST_DATA_CONTAINER, metadata=metadata_file, abundance_path=abundance_file, taxonomic_path=taxonomy_file)

    Data = DataStorage(None, hdf5_file_path, taxonomy_provided, abundance_provided)

    main_dataframe = Data.get_main_dataframe()
    metadata = Data.get_metadata()
    global_production_dataframe = Data.get_global_production_dataframe()
    metabolite_production_dataframe = Data.get_metabolite_production_dataframe()

    assert any(global_production_dataframe["Total_abundance_weighted"].isna),"Nan value in global_production dataframe[Total_abundance_w]"
    assert any(global_production_dataframe["Total_production"].isna),"Nan value in global_production dataframe[Total_production]"

    metadata_against_global_dataframe = metadata.loc[metadata["smplID"].isin(global_production_dataframe["smplID"])]
    assert len(metadata_against_global_dataframe) == len(global_production_dataframe), "len() of metadata and global dataframe not the same after .loc"
    assert global_production_dataframe[:,2:].equals(metadata_against_global_dataframe), "Metadata in global_dataframe and metadata aren't equals."




""" Disable since sample_data is no longer use after data processing.
   
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

    # Producers_long_format tests
    assert len(data_dictionnary["sample_data"].keys()) == len(data_dictionnary["producers_long_format"]["smplID"].unique()), "length of sample_data keys (numbers of samples) and length of producers dataframe smplID.unique() are not the same."

    # Production data tests
    assert all(k in production_data.columns.tolist() for k in data_dictionnary["metadata"].columns.tolist()), "Not all metadata have been inserted in production dataframe."

"""

