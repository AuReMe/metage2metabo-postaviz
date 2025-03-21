import pandas as pd
# import os
from pathlib import Path
import warnings
import plotly.express as px
import plotly
import shutil
import seaborn as sns

from m2m_postaviz import data_utils as du
from m2m_postaviz.data_struct import DataStorage
from m2m_postaviz import shiny_module as sm

# From this test file, get to the test directory then the POSTAVIZ dir.
BASE_DIR = Path(__file__).parent
MODULE_DIR = BASE_DIR.parent

# From Postaviz directory go the m2m_postaviz then postaviz test_data.
SOURCE_DIR = Path(MODULE_DIR, "m2m_postaviz")
TEST_DIR = Path(SOURCE_DIR, "postaviz_test_data")

# If test_data dir doest not contain 'palleja' directory who contain the test data then extract tarfile.
if not Path(TEST_DIR, "palleja/").is_dir():
    # Extracting
    du.extract_tarfile(Path(TEST_DIR, "table_test_postaviz.tar.gz"), TEST_DIR)
    # Load path of newly extracted dir into variable
    TEST_DATA_CONTAINER = Path(TEST_DIR, "palleja/")

else:
    print("Palleja test directory exist.")
    TEST_DATA_CONTAINER = Path(TEST_DIR, "palleja/")

TMP_DIR = Path(TEST_DATA_CONTAINER,"test_save_dir")

# If the directory already exist, directory must be removed with all content. This is usefull for local test since most of data-processing is ignored due to previous local test that already created files.

if TMP_DIR.is_dir():
    shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir()
else:
    TMP_DIR.mkdir(parents=True)
    
metadata_file = Path(TEST_DATA_CONTAINER, "metadata_test_data.tsv")
taxonomy_file = Path(TEST_DATA_CONTAINER, "taxonomy_test_data.tsv")
abundance_file = Path(TEST_DATA_CONTAINER, "abundance_test_data.tsv")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    du.build_dataframes(TEST_DATA_CONTAINER,metadata_file,abundance_file,taxonomy_file,TMP_DIR)



def test_data_processing():

    data = DataStorage(TMP_DIR)

    # sample_info = du.load_sample_cscope_data(TEST_DATA_CONTAINER, CSCOPE_DIR, ".parquet.gz")

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
    
    production_histplot2, production_boxplot2, df2, timing2, abundance_plot2 = sm.bin_exploration_processing(data,
                                                                                                    "Group",
                                                                                                    ["Control","Treatment"],
                                                                                                    "all",
                                                                                                    "Clostridia",
                                                                                                    True, "Days")
    
    production_histplot3, production_boxplot3, df3, timing3, abundance_plot3 = sm.bin_exploration_processing(data,
                                                                                                    "Group",
                                                                                                    ["Control","Treatment"],
                                                                                                    "mgs",
                                                                                                    "MB2bin10 taxonomy",
                                                                                                    True, "Days")
    
    # Object type check.

    assert isinstance(production_histplot, plotly.graph_objs._figure.Figure), "production_histplot is not a plotly express histplot"

    assert isinstance(production_boxplot, plotly.graph_objs._figure.Figure), "Production boxplot is not a plotly express boxplot"

    assert isinstance(df, pd.DataFrame), "Dataframe returned by bin_exploration_processing is not a pandas dataframe."

    assert isinstance(abundance_plot, plotly.graph_objs._figure.Figure), "Abundance barplot is not a plotly express barplot"

    assert isinstance(production_histplot2, plotly.graph_objs._figure.Figure), "production_histplot2 is not a plotly express histplot"

    assert isinstance(production_boxplot2, plotly.graph_objs._figure.Figure), "Production boxplot2 is not a plotly express boxplot"

    assert isinstance(df2, pd.DataFrame), "Dataframe returned by bin_exploration_processing is not a pandas dataframe."

    assert isinstance(abundance_plot2, plotly.graph_objs._figure.Figure), "Abundance2 barplot is not a plotly express barplot"

    assert isinstance(production_histplot3, plotly.graph_objs._figure.Figure), "production_histplot3 is not a plotly express histplot"

    assert isinstance(production_boxplot3, plotly.graph_objs._figure.Figure), "Production boxplot3 is not a plotly express boxplot"

    assert isinstance(df3, pd.DataFrame), "Dataframe returned by bin_exploration_processing is not a pandas dataframe."

    assert abundance_plot3 is None, "Abundance3 barplot should be Nonetype."

    # Object is empty check.

    assert production_histplot.data != tuple(), "Production histogram is empty."

    assert production_boxplot.data != tuple(), "Production boxplot is empty."

    assert abundance_plot.data != tuple(), "Abundance barplot is empty."

    assert production_histplot2.data != tuple(), "Production histogram is empty."

    assert production_boxplot2.data != tuple(), "Production boxplot is empty."

    assert abundance_plot.data != tuple(), "Abundance barplot is empty."

    # Getter test of DataStorage object.

    list_of_bins = data.get_bins_list()

    assert list_of_bins, "List of bins is empty."

    taxonomic_rank = data.get_taxonomy_rank()

    assert taxonomic_rank, "List of taxonomic rank is empty."

    converted_bin_list = data.associate_bin_taxonomy(list_of_bins)

    assert converted_bin_list, "List of bins associated with their taxonomy is empty."

    list_of_factors = data.get_factors()

    assert list_of_factors, "List of metadata factors is empty."

    # Test for custom PCOA

    with warnings.catch_warnings():

        warnings.simplefilter("ignore")

        custom_pcoa_category_factor = sm.make_pcoa(data, "Group", ["Control","Treatment"], True, "Days")

        custom_pcoa_integer_factor = sm.make_pcoa(data, "Days", [0, 180], False, "Group")

    assert custom_pcoa_category_factor.data != tuple(), "Custom pcoa function returned empty plot"

    assert custom_pcoa_integer_factor.data != tuple(), "Custom pcoa function returned empty plot"

    # Test for total production reactive plot.

    reactive_total_production_plot = sm.render_reactive_total_production_plot(data, "Group", "Days", False)
    reactive_total_production_plot_abundance = sm.render_reactive_total_production_plot(data, "Group", "Days", True)
    reactive_total_production_plot_2 = sm.render_reactive_total_production_plot(data, "Group", "Group", False)

    assert isinstance(reactive_total_production_plot, plotly.graph_objs._figure.Figure), "reactive_total_production_plot is supposed to be a plotly graph object."

    assert isinstance(reactive_total_production_plot_abundance, plotly.graph_objs._figure.Figure), "reactive_total_production_plot is supposed to be a plotly graph object."

    assert isinstance(reactive_total_production_plot_2, plotly.graph_objs._figure.Figure), "reactive_total_production_plot2 is supposed to be a plotly graph object."

    assert reactive_total_production_plot.data != tuple(), "reactive_total_production_plot data empty."

    assert reactive_total_production_plot_abundance.data != tuple(), "reactive_total_production_plot_abundance data empty."

    assert reactive_total_production_plot_2.data != tuple(), "reactive_total_production_plot_2 data empty."

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
        assert total_production_test_dataframe == None, "global_production_statistical_dataframe function should have return None with user_input1 == None."

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

    # Test Save function with DF.

    data.save_dataframe(metabolites_production_test_dataframe, "test_save_producers_stat_test")

    saved_path = Path(data.output_path, "test_save_producers_stat_test.tsv")

    if not saved_path.is_file():
        raise AssertionError("Save dataframe function -- File do not exist in: %s" % saved_path)
    # assert isinstance(total_production_test_dataframe, pd.DataFrame) or total_production_test_dataframe == None, "Total production dataframe statistical test is not None or a pandas dataframe."

    # assert isinstance(metabolites_production_test_dataframe, pd.DataFrame) or metabolites_production_test_dataframe == None, "Metabolites production dataframe statistical test is not None or a pandas dataframe."

def test_recursive_tree_padmet():


    from collections import Counter

    data = DataStorage(TMP_DIR)

    if data.USE_METACYC_PADMET:

        dataframe = pd.read_csv(Path(TEST_DIR,"padmet_dataframe_unit_test.tsv"), sep="\t")

        tree = {}

        tree["FRAMES"] = {}

        du.build_tree_from_root(tree["FRAMES"], "FRAMES", dataframe)

        res = data.get_all_tree_keys(tree)

        cpd_list = []

        data.get_compounds_from_category(tree, cpd_list)

        expected_cpd_list = ['CPD-C', 'CPD-E', 'CPD-D', 'CPD-A', 'CPD-B', 'Proton-pump', 'Proton-pump-inhibitor', 'Amylase']

        expected_keys_list = ['FRAMES', 'Sucre', 'hexose', 'pentose', 'Lipide', 'acide-gras', 'acide-gras-non-saturé', 'Proteins', 'Pumps', 'Enzymes']

        assert Counter(res) == Counter(expected_keys_list), "all keys found in padmet tree are not equal to the expected keys list."

        assert Counter(cpd_list) == Counter(expected_cpd_list), "all compounds found in padmet tree are not equal to the expected compounds list."

        list_of_category = data.get_metacyc_category_list()

        assert list_of_category, "List of category is empty."

def test_compounds_exploration_module():

    data = DataStorage(TMP_DIR)

    cpd_short_list = data.get_compound_list()[:50]

    ### Heatmap

    cscope_hm, icscope_hm, added_value_hm = sm.sns_clustermap(data, cpd_short_list, "Days", False, False, "All", None)

    assert isinstance(cscope_hm, sns.matrix.ClusterGrid), "cscope heatmap is not a plotly graph object."

    # assert cscope_hm.data != tuple(), "cscope heatmap is empty."

    assert isinstance(icscope_hm, sns.matrix.ClusterGrid), "iscope heatmap is not a plotly graph object."

    # assert icscope_hm.data != tuple(), "iscope heatmap is empty."

    assert isinstance(added_value_hm, sns.matrix.ClusterGrid), "added_value heatmap is not a plotly graph object."

    # assert added_value_hm.data != tuple(), "added_value heatmap is empty."

    ### Sample producers graph.

    cscope_producers, iscope_producers = sm.percentage_smpl_producing_cpd(data, cpd_short_list, "Days")

    assert isinstance(cscope_producers, plotly.graph_objs._figure.Figure), "Percent of sample producers cscope graph is not an plotly object"

    assert cscope_producers.data != tuple(), "smpl_prod cscope graph empty."

    assert isinstance(iscope_producers, plotly.graph_objs._figure.Figure), "Percent of sample producers iscope graph is not an plotly object"

    assert iscope_producers.data != tuple(), "smpl_prod iscope graph empty."