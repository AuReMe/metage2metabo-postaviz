import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
import timeit
from m2m_postaviz.data_utils import benchmark_decorator
import cProfile


class DataStorage:
    column_identifier = "smplID"

    def __init__(
        self,
        main_data_container: dict,
        sample_data_container: dict,
        taxonomic_data_container: pd.DataFrame,
        norm_abundance_data: pd.DataFrame,
        abundance_data: pd.DataFrame,
    ):
        self.main_data = main_data_container
        self.sample_data = sample_data_container

        self.taxonomic_data = taxonomic_data_container
        self.long_taxo_data = self.taxonomic_data_long_format()

        self.normalised_abundance_matrix = norm_abundance_data
        self.abundance_matrix = abundance_data

        self.abundance_dataframe: pd.DataFrame = self.main_data["metadata"].merge(self.abundance_matrix.reset_index(), how="outer")
        self.normalised_abundance_dataframe: pd.DataFrame = self.main_data["metadata"].merge(
            self.normalised_abundance_matrix.reset_index(), how="outer"
        )

        self.melted_abundance_dataframe: pd.DataFrame = self.produce_long_abundance_dataframe(with_normalisation=False)
        self.melted_normalised_abundance_dataframe: pd.DataFrame = self.produce_long_abundance_dataframe()

        self.list_of_factor = list(self.main_data["metadata"].columns)
        self.factorize_metadata()

    def performance_benchmark(self):
        cProfile.runctx("self.taxonomic_data_long_format()", globals(), locals())

    def get_sample_data(self, smpl_id, mode: str = "cscope", as_copy: bool = True):
        return self.sample_data[smpl_id][mode].copy()

    def get_main_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.main_data["main_dataframe"].copy() if as_copy else self.main_data["main_dataframe"]

    def get_main_metadata(self, as_copy: bool = True) -> pd.DataFrame:
        return self.main_data["metadata"].copy() if as_copy else self.main_data["metadata"]

    def get_long_taxonomic_data(self, as_copy: bool = True) -> pd.DataFrame:
        return self.long_taxo_data.copy() if as_copy else self.long_taxo_data

    def get_taxonomic_data(self, as_copy: bool = True) -> pd.DataFrame:
        return self.taxonomic_data.copy() if as_copy else self.taxonomic_data

    def get_norm_ab_matrix(self, as_copy: bool = True) -> pd.DataFrame:
        return self.normalised_abundance_matrix.copy() if as_copy else self.normalised_abundance_matrix

    def get_norm_ab_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.normalised_abundance_dataframe.copy() if as_copy else self.normalised_abundance_dataframe

    def get_ab_matrix(self, as_copy: bool = True) -> pd.DataFrame:
        return self.abundance_matrix.copy() if as_copy else self.abundance_matrix

    def get_ab_dataframe(self, as_copy: bool = True) -> pd.DataFrame:
        return self.abundance_dataframe.copy() if as_copy else self.abundance_dataframe

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
        return len(self.main_data["metadata"].columns)

    def get_metadata_label(self):
        if self.is_indexed(self.main_data["metadata"]):
            return list(self.main_data["metadata"].columns)
        return list(self.main_data["metadata"].columns)[1:]

    def get_factor_by_id(self, sample_id, factor):
        query = self.main_data["metadata"].loc[self.main_data["metadata"]["smplID"] == sample_id][factor]
        print("get_factor : ", sample_id, factor, "found : ", query)
        return query

    def list_to_boolean_serie(self, model_list: list, with_quantity: bool = True):
        results = {}
        value = 1
        for model in model_list:
            if not model in results.keys():
                results[model] = value
            else:
                if with_quantity:
                    results[model] += value
        return pd.Series(results)

    def taxonomy_matrix_build(self):
        rank = "s"
        id_col = "mgs"
        all_series = {}
        taxonomic_df = self.get_taxonomic_data()
        binlist_by_id = self.get_bin_list()

        for sample in binlist_by_id.keys():
            res = taxonomic_df.loc[taxonomic_df[id_col].isin(binlist_by_id[sample])][rank]
            all_series[sample] = self.list_to_boolean_serie(res)

        matrix = pd.DataFrame(all_series)
        matrix.fillna(0, inplace=True)
        matrix = matrix.T
        matrix.index.name = "smplID"
        matrix.reset_index(inplace=True)
        # !!! NAME OF ID isnt SMPLID !!! CAN LEAD TO DRAMA
        return matrix

    def taxonomic_data_long_format(self):
        """
        Produce long format taxonomic dataframe for plot purpose.
        Returns:
            Datarame: Dataframe in long format
        """

        df = self.taxonomy_matrix_build()
        # df_matrix = df.set_index('smplID')
        # print(df_matrix.astype(int).agg('sum',axis=1))
        # quit()
        if self.is_indexed(df):
            df.reset_index()

        df = df.melt("smplID", var_name="Taxa", value_name="Quantity")
        brand_new_df = {}

        for factor in self.get_metadata_label():
            brand_new_df[factor] = self.add_factor_column(df["smplID"], factor)

        df = df.assign(**brand_new_df)
        df = df.astype(str)
        df["smplID"] = df["smplID"].astype('category')
        df['Quantity'] = df['Quantity'].astype(float)
        df['Nb_taxon'] = df["smplID"].apply(lambda row: self.search_long_format(row,df))
        return df
    
    def search_long_format(self, id_value, df):
        value = df.loc[(df['smplID'] == id_value) & (df['Quantity'] != 0)]["Taxa"].unique()
        return len(value)

    def produce_long_abundance_dataframe(self, with_normalisation: bool = True):
        """Transform the wide format abundance dataframe into a long format.
        Usefull for anyplot.
        Args:
            with_normalisation (bool, optional): Which dataframe to transform. Defaults to True.

        Returns:
            pd.dataframe: Long format abundance dataframe.
        """
        if with_normalisation:
            current_df = self.get_norm_ab_matrix()
        else:
            current_df = self.get_ab_matrix()
        if self.is_indexed(current_df):
            current_df = current_df.reset_index()
        current_df = current_df.melt("smplID", var_name="Compound", value_name="Quantity")
        all_new_columns = {}
        for factor in self.get_metadata_label():
            all_new_columns[factor] = self.add_factor_column(current_df["smplID"], factor)

        current_df = current_df.assign(**all_new_columns)
        return current_df

    def add_factor_column(self, serie_id, factor_id):
        metadata = self.get_main_metadata()
        metadata.set_index("smplID", inplace=True, drop=True)
        new_col = []
        for value in serie_id:
            new_col.append(str(metadata.at[value, factor_id]))
        return new_col

    def factorize_metadata(self):
        for factor in self.get_metadata_label():
            self.main_data["metadata"][factor] = self.main_data["metadata"][factor].astype("category")

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
