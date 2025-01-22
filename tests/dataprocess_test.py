import pandas as pd
import os
import warnings
import plotly.express as px
import plotly
import shutil

from m2m_postaviz import data_utils as du
from m2m_postaviz.data_struct import DataStorage
from m2m_postaviz import shiny_module as sm

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

TMP_DIR = os.path.join(TEST_DATA_CONTAINER,"test_save_dir")

# If the directory already exist, directory must be removed with all content. This is usefull for local test since most of data-processing is ignored due to previous local test that already created files.

if os.path.isdir(TMP_DIR):
    shutil.rmtree(TMP_DIR)
    os.makedirs(TMP_DIR)

metadata_file = os.path.join(TEST_DATA_CONTAINER, "metadata_test_data.tsv")
taxonomy_file = os.path.join(TEST_DATA_CONTAINER, "taxonomy_test_data.tsv")
abundance_file = os.path.join(TEST_DATA_CONTAINER, "abundance_test_data.tsv")

CSCOPE_DIR = os.path.join(TEST_DATA_CONTAINER, "cscope_directory")

os.makedirs(CSCOPE_DIR)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    du.build_dataframes(TEST_DATA_CONTAINER,metadata_file,abundance_file,taxonomy_file,TMP_DIR)



def test_data_processing():

    data = DataStorage(TMP_DIR)

    sample_info = du.load_sample_cscope_data(TEST_DATA_CONTAINER, CSCOPE_DIR, ".parquet.gz")

    metadata = data.get_metadata()

    main_dataframe = data.get_main_dataframe()

    metabolite_production_dataframe = data.get_metabolite_production_dataframe()

    global_production_dataframe = data.get_global_production_dataframe()

    # Global_production_dataframe are not supposed to contain NaN values in their processed columns.
    assert any(global_production_dataframe["Total_abundance_weighted"].notna()),"Nan value in global_production dataframe[Total_abundance_w]"

    assert any(global_production_dataframe["Total_production"].notna()),"Nan value in global_production dataframe[Total_production]"

    # Select rows of metadata matching with global_dataframe, just in case they weren't a perfect match (metadata is bigger i think)
    metadata_against_global_dataframe = metadata.loc[metadata["smplID"].isin(global_production_dataframe["smplID"])]

    assert len(metadata_against_global_dataframe) == len(global_production_dataframe), "len() of metadata and global dataframe not the same after .loc"

    # Set sample id as index
    global_production_dataframe.set_index("smplID",inplace=True)
    metadata_against_global_dataframe.set_index("smplID",inplace=True)
    # Sort index to align both dataframe
    global_production_dataframe.sort_index(inplace=True)
    metadata_against_global_dataframe.sort_index(inplace=True)

    # After removing Total_production and Total_abundance_weighted column from global, both dataframe are supposed to be the same. This is a check for the merge which showed weird behavior.
    print(global_production_dataframe.drop(columns=["Total_production", "Total_abundance_weighted"]).compare(metadata_against_global_dataframe,align_axis=0))
    assert global_production_dataframe.drop(columns=["Total_production", "Total_abundance_weighted"]).equals(metadata_against_global_dataframe), "Metadata in global_dataframe and metadata aren't equals."

    # Check if all sample_data are dataframe and if the numbers of metabolites produced is the same by checking len(columns in df) and len(metabolites produced list in json)
    # for sample in sample_data.keys():
    #     # JSON equivalent of TSV format cscope. usefull to check if any values has been lost in the processing.
    #     sample_json = du.open_json(os.path.join(TEST_DATA_CONTAINER, sample+"/"+"community_analysis/"+"comm_scopes.json"))["com_scope"]

    #     assert isinstance(sample_data[sample]["cscope"], pd.DataFrame), f"sample {sample} in sample_data is not a pandas dataframe object."
    #     assert len(sample_data[sample]["cscope"].columns) == len(sample_json),"Length of sample columns and length of json metabolites production list are not the same."

        # Select 1 row of main_dataframe by matching index ID with the sample name. Then extract the values and sum() the get the amount of cpd produced. Compare first with TSV cscope then JSON file.
        # if not du.is_indexed_by_id(main_dataframe):
        #     main_dataframe.set_index("smplID",inplace=True)
        # assert main_dataframe.loc[main_dataframe.index == sample].values.sum() == len(sample_json),f"The amount of metabolites produced for {sample} in main_dataframe is not the same as com_scopes json file."


""" Disable since sample_data is no longer use after data processing.


    # Producers_long_format tests
    assert len(sample_data.keys()) == len(data_dictionnary["producers_long_format"]["smplID"].unique()), "length of sample_data keys (numbers of samples) and length of producers dataframe smplID.unique() are not the same."

    # Production data tests
    assert all(k in production_data.columns.tolist() for k in data_dictionnary["metadata"].columns.tolist()), "Not all metadata have been inserted in production dataframe."

"""

def test_query_parquet():

    data = DataStorage(TMP_DIR)

    taxonomic_rank_input = "c"

    taxonomic_rank_unique_input = "Fusobacteriia"

    list_of_bin_in_rank = data.get_bin_list_from_taxonomic_rank(taxonomic_rank_input, taxonomic_rank_unique_input)
    
    query = [("binID", "in", list_of_bin_in_rank)]

    df = data.get_bin_dataframe(condition=query)

    # Testing dataframe. 5 bins in query

    assert isinstance(df,pd.DataFrame), "bin_dataframe is not an instance of pandas dataframe."

    assert df.empty == False, "bin_dataframe is empty."

    assert all(df['c'].to_numpy() == "Fusobacteriia"), f"Bin_dataframe with rank choice of {taxonomic_rank_input} should only contain {taxonomic_rank_unique_input}."

    
def test_shiny_module():

    data = DataStorage(TMP_DIR)

    # Test for bins exploration tab

    production_histplot, production_boxplot, df, timing, abundance_plot = sm.bin_exploration_processing(data,
                                                                                                        "Group",
                                                                                                        ["Control","Treatment"],
                                                                                                        "c",
                                                                                                        "Clostridia",
                                                                                                        True, "Days")
    
    # Object type check.

    assert isinstance(production_histplot, plotly.graph_objs._figure.Figure), "Production histogram is not a plotly express histplot"

    assert isinstance(production_boxplot, plotly.graph_objs._figure.Figure), "Production boxplot is not a plotly express boxplot"

    assert isinstance(df, pd.DataFrame), "Dataframe returned by bin_exploration_processing is not a pandas dataframe."

    assert isinstance(abundance_plot, plotly.graph_objs._figure.Figure), "Abundance barplot is not a plotly express barplot"

    # Object is empty check.

    assert production_histplot != tuple(), "Production histogram is empty."

    assert production_boxplot != tuple(), "Production boxplot is empty."

    assert abundance_plot != tuple(), "Abundance barplot is empty."


    # Test for custom PCOA

    with warnings.catch_warnings():

        warnings.simplefilter("ignore")

        custom_pcoa_category_factor = sm.make_pcoa(data, "Group", ["Control","Treatment"], True, "Days")

        custom_pcoa_integer_factor = sm.make_pcoa(data, "Days", [0, 180], False, "Group")

    assert custom_pcoa_category_factor != tuple(), "Custom pcoa function returned empty plot"

    assert custom_pcoa_integer_factor != tuple(), "Custom pcoa function returned empty plot"

    # Test for total production reactive plot.

    reactive_total_production_plot = sm.render_reactive_total_production_plot(data, "Group", "Days", True)
    reactive_total_production_plot_abundance = sm.render_reactive_total_production_plot(data, "Group", "Days", True)

    assert isinstance(reactive_total_production_plot, plotly.graph_objs._figure.Figure), "reactive_total_production_plot is supposed to be a plotly graph object."

    assert isinstance(reactive_total_production_plot_abundance, plotly.graph_objs._figure.Figure), "reactive_total_production_plot is supposed to be a plotly graph object."

    # Test for metabolites production reactive plot.

    reactive_metabolites_production_plot = sm.render_reactive_metabolites_production_plot(data, ["CPD-15709[c]", "CPD-372[c]"], "Group", "Days")
    reactive_metabolites_production_plot_abundance = sm.render_reactive_metabolites_production_plot(data, ["CPD-15709[c]", "CPD-372[c]"],"Group", "Days")

    assert isinstance(reactive_metabolites_production_plot, plotly.graph_objs._figure.Figure), "reactive_metabolites_production_plot is supposed to be a plotly graph object."

    assert isinstance(reactive_metabolites_production_plot_abundance, plotly.graph_objs._figure.Figure), "reactive_metabolites_production_plot is supposed to be a plotly graph object."


def test_statistic_method():
    
    data = DataStorage(TMP_DIR)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Return None
        total_production_test_dataframe = sm.global_production_statistical_dataframe(data, "None", "Days", True, "simes-hochberg", True)
        assert total_production_test_dataframe == None, "global_production_statistical_dataframe function should've return None with user_input1 == None."

        # Return Wilcoxon/Man-whitney dataframe
        total_production_test_dataframe = sm.global_production_statistical_dataframe(data, "Group", "Days", True, "simes-hochberg", True)
        assert total_production_test_dataframe["Method"].unique() == "Wilcoxon", "global_production_statistical_dataframe function should've return dataframe with Wilcoxon test method"

        # Return Correlation dataframe
        total_production_test_dataframe = sm.global_production_statistical_dataframe(data, "Days", "Group", True, "simes-hochberg", True)
        assert total_production_test_dataframe["Method"].unique() == "pearson", "global_production_statistical_dataframe function should've return dataframe with Pearson test method"

        # Return Wilcoxon/Man-whitney dataframe
        metabolites_production_test_dataframe = sm.metabolites_production_statistical_dataframe(data, ["CPD-15709[c]", "CPD-372[c]"], "Group", "None", True, "simes-hochberg", True)
        assert metabolites_production_test_dataframe["Method"].unique() == "Mann-Whitney", "metabolites_production_statistical_dataframe function should've return dataframe with Mann-Whitney test method"
        
        # Return Correlation dataframe        
        metabolites_production_test_dataframe = sm.metabolites_production_statistical_dataframe(data, ["CPD-15709[c]", "CPD-372[c]"], "Days", "None", True, "simes-hochberg", True )
        assert metabolites_production_test_dataframe["Method"].unique() == "pearson", "metabolites_production_statistical_dataframe function should've return dataframe with Pearson test method"

        # Return Wilcoxon/Man-whitney dataframe
        metabolites_production_test_dataframe = sm.metabolites_production_statistical_dataframe(data, ["CPD-15709[c]", "CPD-372[c]"], "Group", "Days", True, "simes-hochberg", True )
        assert metabolites_production_test_dataframe["Method"].unique() == "Wilcoxon", "metabolites_production_statistical_dataframe function should've return dataframe with Wilcoxon test method"

        # Return Correlation dataframe
        metabolites_production_test_dataframe = sm.metabolites_production_statistical_dataframe(data, ["CPD-15709[c]", "CPD-372[c]"], "Days", "Group", True, "simes-hochberg", True )
        assert metabolites_production_test_dataframe["Method"].unique() == "pearson", "metabolites_production_statistical_dataframe function should've return dataframe with Pearson test method"

    # assert isinstance(total_production_test_dataframe, pd.DataFrame) or total_production_test_dataframe == None, "Total production dataframe statistical test is not None or a pandas dataframe."

    # assert isinstance(metabolites_production_test_dataframe, pd.DataFrame) or metabolites_production_test_dataframe == None, "Metabolites production dataframe statistical test is not None or a pandas dataframe."
