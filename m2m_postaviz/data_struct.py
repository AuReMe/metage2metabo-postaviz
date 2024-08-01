import cProfile
import plotly.express as px
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
import os
import sys
import m2m_postaviz.data_utils as du
import plotly.graph_objects as go

class DataStorage:

    ID_VAR = "smplID"
    HAS_TAXONOMIC_DATA : bool = False
    HAS_ABUNDANCE_DATA : bool = False

    def __init__(self, save_path: str, data_container: dict, taxonomic_data_container: pd.DataFrame = None, abundance_data: pd.DataFrame = None, total_production_dataframe: pd.DataFrame = None, pcoa_dataframe: pd.DataFrame = None):

        self.main_data = data_container["main_dataframe"]
        self.sample_data = data_container["sample_data"]
        self.metadata = data_container["metadata"]
        
        self.producer_dataframe = data_container["producers_long_format"]
        self.total_production_dataframe = total_production_dataframe
        
        if du.is_valid_dir(save_path):
            self.output_path = save_path
        else:
            try:
                os.makedirs(save_path)
                self.output_path = save_path
            except Exception as e:
                print(e)
                sys.exit(1)

        if taxonomic_data_container is not None:
            self.long_taxonomic_data = taxonomic_data_container
            self.HAS_TAXONOMIC_DATA = True

        if abundance_data is not None:
            self.normalised_abundance_matrix = abundance_data
            self.HAS_ABUNDANCE_DATA = True
            self.melted_normalised_abundance_dataframe: pd.DataFrame = self.produce_long_abundance_dataframe()

        self.list_of_factor = list(self.metadata.columns)
        # self.factorize_metadata()

        self.current_pcoadf = pcoa_dataframe

    def performance_benchmark(self):
        cProfile.runctx("self.taxonomic_data_long_format()", globals(), locals())

    def get_total_production_by_sample(self, as_copy: bool = True):
        return self.total_production_dataframe.copy() if as_copy else self.total_production_dataframe

    def get_producer_dataframe(self, as_copy: bool = True):
        return self.producer_dataframe.copy() if as_copy else self.producer_dataframe

    def get_cpd_list(self):
        query = self.get_main_dataframe().columns.tolist()
        if "smplID" in query:
            query.remove("smplID")
        return query

    def get_all_sample_data(self, mode: str = "cscope", as_copy: bool = True):
        return self.sample_data.copy() if as_copy else self.sample_data

    def get_sample_data(self, smpl_id, mode: str = "cscope", as_copy: bool = True):
        return self.sample_data[smpl_id][mode].copy()

    def get_main_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.main_data.copy() if as_copy else self.main_data

    def get_main_metadata(self, as_copy: bool = True) -> pd.DataFrame:
        return self.metadata.copy() if as_copy else self.metadata

    def set_main_metadata(self, new_metadata):
        self.metadata = new_metadata

    def get_long_taxonomic_data(self, as_copy: bool = True) -> pd.DataFrame:
        return self.long_taxonomic_data.copy() if as_copy else self.long_taxonomic_data

    def get_norm_ab_matrix(self, as_copy: bool = True) -> pd.DataFrame:
        return self.normalised_abundance_matrix.copy() if as_copy else self.normalised_abundance_matrix

    def get_norm_ab_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.normalised_abundance_dataframe.copy() if as_copy else self.normalised_abundance_dataframe

    ### LATER STILL USEFULL ?
    def get_melted_norm_ab_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.melted_normalised_abundance_dataframe.copy() if as_copy else self.melted_normalised_abundance_dataframe

    def is_indexed(self, df: pd.DataFrame) -> bool:
        return True if df.index.name == "smplID" else False

    def get_bin_list_by_sample(self, sample_id: str, mode: str = "cscope"):
        return self.sample_data[sample_id][mode]["Name"].to_list()

    def get_cpd_label(self, sample_id: str, mode: str = "cscope"):
        return self.sample_data[sample_id][mode].columns

    def get_bin_list(self, mode: str = "cscope"):
        bin_list = {}
        for sample in self.sample_data.keys():
            if self.is_indexed(self.sample_data[sample][mode]):
                bin_list[sample] = self.sample_data[sample][mode].index.to_list()
            else:
                bin_list[sample] = self.sample_data[sample][mode][self.ID_VAR].to_list()
        return bin_list

    def get_factor_len(self):
        return len(self.metadata.columns)

    def get_metadata_label(self):
        print(self.metadata.columns)
        if self.is_indexed(self.metadata):
            return list(self.metadata.columns)
        print(self.metadata.columns[1:])
        return list(self.metadata.columns[1:])

    def get_factor_by_id(self, sample_id, factor):
        query = self.metadata.loc[self.metadata["smplID"] == sample_id][factor]
        print("get_factor : ", sample_id, factor, "found : ", query)
        return query

    def produce_long_abundance_dataframe(self):
        """Transform the wide format abundance dataframe into a long format.
        Usefull for anyplot.

        Returns:
            pd.dataframe: Long format abundance dataframe.
        """

        current_df = self.get_norm_ab_matrix()

        if self.is_indexed(current_df):
            current_df = current_df.reset_index()
        current_df = current_df.melt("smplID", var_name="Compound", value_name="Quantity")
        all_new_columns = {}

        metadata_label = self.get_metadata_label()
        for factor in metadata_label:
            all_new_columns[factor] = du.add_factor_column(self.get_main_metadata(), current_df["smplID"], factor)

        current_df = current_df.assign(**all_new_columns)
        current_df = current_df.loc[current_df["Quantity"] != 0]
        return current_df

    def factorize_metadata(self):
        for factor in self.get_metadata_label():
            self.metadata[factor] = self.metadata[factor].astype("category")

    def weird_way_to_do_it(self, id_value: str, metadata_col: str, metadata: pd.DataFrame):
        if self.is_indexed(metadata):
            metadata.reset_index(inplace=True)
        result = metadata.loc[metadata["smplID"] == id_value][metadata_col].values[0]
        return result

    def save_dataframe(self, df_to_save:pd.DataFrame, file_name: str, extension: str = "tsv"):
        path_to_save = self.output_path
        final_file_path = path_to_save + "/" + file_name + "." + extension
        if os.path.isfile(final_file_path):
            final_file_path = self.check_and_rename(final_file_path)
        try:
            df_to_save.to_csv(final_file_path)
            logs = f"Saved in :\n{final_file_path}"
        except Exception as e:
            logs = e
        return logs
    
    def check_and_rename(self, file_path: str, add: int = 0) -> str:
        original_file_path = file_path
        print(original_file_path)
        if add != 0:
            file_name, extension = os.path.splitext(file_path)
            file_name = file_name + "_" + str(add)
            file_path = file_name + extension
        if not os.path.isfile(file_path):
            return file_path
        else:
            return self.check_and_rename(original_file_path, add + 1)
        
    def make_figure_from_df(df: pd.DataFrame):
        fig = go.Figure()
        fig.add_trace(go.scatter(
                df,
                x='PC1',
                y='PC2',
                color=df.index
            )
        )
        return fig
    
    def plot(df: pd.DataFrame, compounds_input: str, first_input:str = "None", second_input = "None"):
        """Return a plotly express figure from the dataframe and the input(s).
        Checks for y single value case for barplot or boxplot.

        Args:
            df (pd.DataFrame): dataframe to plot
            compounds_input (str): list of compounds to show
            first_input (str, optional): first axis of the plot, must be a label of a column of the dataframe. Defaults to "None".
            second_input (str, optional): second axis of the plot, must be a column's label of the datafrale. Defaults to "None".

        Returns:
            px.figure: plotly express figure.
        """
        if first_input == "None":
            return px.box(df,y=compounds_input, notched=True)

        if second_input == "None":
            has_unique_value = True
            for value_input1 in df[first_input].unique():   
                print(len(df.loc[df[first_input] == value_input1][compounds_input]))
                if len(df.loc[df[first_input] == value_input1][compounds_input]) > 1:
                    has_unique_value = False
                    break

            return px.bar(df,x=first_input,y=compounds_input,color=first_input) if has_unique_value else px.box(df,x=first_input,y=compounds_input,color=first_input, notched=True)
            
        has_unique_value = True
        for value_input1 in df[first_input].unique():
                for value_input2 in df[second_input].unique():   
                    if len(df.loc[(df[second_input] == value_input2) & (df[first_input] == value_input1)][compounds_input]) > 1:
                        has_unique_value = False
                        break

        return px.bar(df, x=first_input, y=compounds_input,color=second_input) if has_unique_value else px.box(df, x=first_input, y=compounds_input,color=second_input, boxmode="group")

        