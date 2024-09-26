import pandas as pd
import os
import warnings

from m2m_postaviz import data_utils as du
from m2m_postaviz.data_struct import DataStorage

def test_data_processing():

    # From this test file, get to the test directory then the POSTAVIZ dir.
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    MODULE_DIR = os.path.dirname(BASE_DIR)

    # From Postaviz directory go the m2m_postaviz then postaviz test_data.
    SOURCE_DIR = os.path.join(MODULE_DIR, "m2m_postaviz")
    TEST_DIR = os.path.join(SOURCE_DIR, "postaviz_test_data")

    # If test_data dir doest not contain 'palleja' directory who contain the test data then extract tarfile.
    if not os.path.isdir(os.path.join(TEST_DIR, "palleja/")):
        # Extracting
        du.extract_tarfile(os.path.join(TEST_DIR, "table_test_postaviz.tar.gz"), TEST_DIR)
        # Load path of newly extracted dir into variable
        TEST_DATA_CONTAINER = os.path.join(TEST_DIR, "palleja/")

    else:
        print("Palleja test directory exist.")
        TEST_DATA_CONTAINER = os.path.join(TEST_DIR, "palleja/")

    metadata_file = os.path.join(TEST_DATA_CONTAINER, "metadata_test_data.tsv")
    taxonomy_file = os.path.join(TEST_DATA_CONTAINER, "taxonomy_test_data.tsv")
    abundance_file = os.path.join(TEST_DATA_CONTAINER, "abundance_test_data.tsv")

    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")

    # Loading dataframes into variables

    sample_data = du.multiprocess_retrieve_data(TEST_DATA_CONTAINER)

    metadata = du.open_tsv(metadata_file)

    main_dataframe = du.build_main_dataframe(sample_data)

    metabolite_production_dataframe = du.producers_by_compounds_and_samples_multi(sample_data, metadata)

    global_production_dataframe = du.total_production_by_sample(sample_data, metadata, du.relative_abundance_calc(du.open_tsv(abundance_file),sample_data))


    # Global_production_dataframe are not supposed to contain NaN values in their processed columns.
    assert any(global_production_dataframe["Total_abundance_weighted"].isna),"Nan value in global_production dataframe[Total_abundance_w]"

    assert any(global_production_dataframe["Total_production"].isna),"Nan value in global_production dataframe[Total_production]"

    # Select rows of metadata matching with global_dataframe, just in case they weren't a perfect match (metadata is bigger i think)
    metadata_against_global_dataframe = metadata.loc[metadata["smplID"].isin(global_production_dataframe["smplID"])]

    assert len(metadata_against_global_dataframe) == len(global_production_dataframe), "len() of metadata and global dataframe not the same after .loc"

    # After removing Total_production and Total_abundance_weighted column from global, both dataframe are supposed to be the same. This is a check for the merge which showed weird behavior.
    assert global_production_dataframe[:,2:].equals(metadata_against_global_dataframe), "Metadata in global_dataframe and metadata aren't equals."

    # Check if all sample_data are dataframe and if the numbers of metabolites produced is the same by checking len(columns in df) and len(metabolites produced list in json)
    for sample in sample_data.keys():
        # JSON equivalent of TSV format cscope. usefull to check if any values has been lost in the processing.
        sample_json = du.open_json(os.path.join(TEST_DATA_CONTAINER, sample+"/"+"community_analysis/"+"comm_scopes.json"))["com_scope"]

        assert isinstance(sample_data[sample]["cscope"], pd.DataFrame), f"sample {sample} in sample_data is not a pandas dataframe object."
        assert len(sample_data[sample]["cscope"].columns) == len(sample_json),"Length of sample columns and length of json metabolites production list are not the same."

        # Select 1 row of main_dataframe by matching index ID with the sample name. Then extract the values and sum() the get the amount of cpd produced. Compare first with TSV cscope then JSON file.
        if not du.is_indexed_by_id(main_dataframe):
            main_dataframe.set_index("smplID",inplace=True)
        assert main_dataframe.loc[main_dataframe.index == sample].values.sum() == len(sample_json),f"The amount of metabolites produced for {sample} in main_dataframe is not the same as com_scopes json file."


""" Disable since sample_data is no longer use after data processing.


    # Producers_long_format tests
    assert len(sample_data.keys()) == len(data_dictionnary["producers_long_format"]["smplID"].unique()), "length of sample_data keys (numbers of samples) and length of producers dataframe smplID.unique() are not the same."

    # Production data tests
    assert all(k in production_data.columns.tolist() for k in data_dictionnary["metadata"].columns.tolist()), "Not all metadata have been inserted in production dataframe."

"""

