import cProfile

import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa

import m2m_postaviz.data_utils as du


class DataStorage:
    column_identifier = "smplID"
    HAS_TAXONOMIC_DATA : bool = False
    HAS_ABUNDANCE_DATA : bool = False
    def __init__(self, data_container: dict, taxonomic_data_container: pd.DataFrame = None, abundance_data: pd.DataFrame = None,):

        self.main_data = data_container["main_dataframe"]
        self.sample_data = data_container["sample_data"]
        self.metadata = data_container["metadata"]
        self.cpd_producers_long = data_container["producers_long_format"]

        if taxonomic_data_container is not None:
            self.long_taxonomic_data = taxonomic_data_container
            self.HAS_TAXONOMIC_DATA = True

        if abundance_data is not None:
            self.normalised_abundance_matrix = abundance_data
            self.HAS_ABUNDANCE_DATA = True

        self.melted_normalised_abundance_dataframe: pd.DataFrame = self.produce_long_abundance_dataframe()

        self.list_of_factor = list(self.metadata.columns)
        self.factorize_metadata()

        
        self.compound_production_by_sample = du.total_production_by_sample(self.get_main_dataframe(), self.get_all_sample_data())

    def performance_benchmark(self):
        cProfile.runctx("self.taxonomic_data_long_format()", globals(), locals())

    def get_total_cpd_production_by_sample(self, as_copy: bool = True):
        return self.compound_production_by_sample.copy() if as_copy else self.compound_production_by_sample

    def get_producer_long_dataframe(self, as_copy: bool = True):
        return self.cpd_producers_long.copy() if as_copy else self.cpd_producers_long

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

    def get_long_taxonomic_data(self, as_copy: bool = True) -> pd.DataFrame:
        return self.long_taxonomic_data.copy() if as_copy else self.long_taxonomic_data

    # def get_taxonomic_data(self, as_copy: bool = True) -> pd.DataFrame:
    #     return self.taxonomic_data.copy() if as_copy else self.taxonomic_data

    def get_norm_ab_matrix(self, as_copy: bool = True) -> pd.DataFrame:
        return self.normalised_abundance_matrix.copy() if as_copy else self.normalised_abundance_matrix

    def get_norm_ab_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.normalised_abundance_dataframe.copy() if as_copy else self.normalised_abundance_dataframe

    def get_melted_norm_ab_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.melted_normalised_abundance_dataframe.copy() if as_copy else self.melted_normalised_abundance_dataframe

    def get_melted_ab_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.melted_abundance_dataframe.copy() if as_copy else self.melted_abundance_dataframe

    def is_indexed(self, df: pd.DataFrame) -> bool:
        return True if df.index.name == self.column_identifier else False

    def get_bin_list_by_sample(self, sample_id: str, mode: str = "cscope"):
        return self.sample_data[sample_id][mode]["Name"].to_list()

    def get_cpd_label(self, sample_id: str, mode: str = "cscope"):
        return self.sample_data[sample_id][mode].columns

    def get_taxonomic_data_by_sample(self, sample_id: str) -> pd.DataFrame:
        return self.taxonomic_data.loc[self.taxonomic_data["mgs"].isin(self.get_bin_list_by_sample(sample_id))]

    def get_bin_list(self, mode: str = "cscope"):
        bin_list = {}
        for sample in self.sample_data.keys():
            if self.is_indexed(self.sample_data[sample][mode]):
                bin_list[sample] = self.sample_data[sample][mode].index.to_list()
            else:
                bin_list[sample] = self.sample_data[sample][mode][self.column_identifier].to_list()
        return bin_list

    def get_factor_len(self):
        return len(self.metadata.columns)

    def get_metadata_label(self):
        if self.is_indexed(self.metadata):
            return list(self.metadata.columns)
        return list(self.metadata.columns[1:])

    def get_factor_by_id(self, sample_id, factor):
        query = self.metadata.loc[self.metadata["smplID"] == sample_id][factor]
        print("get_factor : ", sample_id, factor, "found : ", query)
        return query

    def produce_long_abundance_dataframe(self):
        """Transform the wide format abundance dataframe into a long format.
        Usefull for anyplot.
        Args:
            with_normalisation (bool, optional): Which dataframe to transform. Defaults to True.

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

    def run_pcoa(self, main_df: pd.DataFrame, metadata: pd.DataFrame, distance_method: str = "jaccard"):
        """Calculate Principal Coordinate Analysis with the dataframe given in args.
        Use metadata's drataframe as second argument to return the full ordination result plus
        all metadata column inserted along Ordination.samples dataframe.
        Ready to be plotted.

        Args:
            main_df (pd.DataFrame): Main dataframe of compound production
            metadata (pd.DataFrame): Metadata's dataframe

        Returns:
            _type_: Ordination results object from skbio's package.
        """
        df = self.get_main_dataframe()
        # Need the matrix version for distance calculation.
        if not self.is_indexed(main_df):
            main_df.set_index("smplID", inplace=True)

        # Add metadata columns with the accurate value in temporary dataframe.
        for col in metadata.columns:
            if col == "smplID":
                continue
            df[col] = df["smplID"].apply(lambda row: self.weird_way_to_do_it(row, col, metadata))

        # Normalisation
        main_df = main_df.apply(lambda x: x / x.sum(), axis=0)
        # Calculate distance matrix with Bray-Curtis method.
        dist_m = pdist(main_df, distance_method)
        # Transform distance matrix into squareform.
        squaref_m = squareform(dist_m)
        # Run the PCOA with the newly generated distance matrix.
        pcoa_results = pcoa(squaref_m, number_of_dimensions=main_df.shape[0] - 1)

        # Verify if PCOA_results's samples dataframe is aligned with df dataframe.
        df = df.set_index(pcoa_results.samples.index)
        # Put each metadata column in the pcoa results's dataframe for PLOT.
        for col in metadata.columns:
            if col == "smplID":
                continue
            pcoa_results.samples[col] = df[col].astype("category")

        return pcoa_results
