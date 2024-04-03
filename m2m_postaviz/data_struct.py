import pandas as pd


class DataStorage:
    def __init__(self, main_data_container: dict, sample_data_container: dict, taxonomic_data_container: pd.DataFrame, abundance_data: pd.DataFrame):
        self.main_data = main_data_container
        self.sample_data = sample_data_container
        self.taxonomic_data = taxonomic_data_container
        self.abundance_matrix = abundance_data
        self.list_of_factor = list(self.main_data["metadata"].columns)
        self.abundance_dataframe : pd.DataFrame = self.main_data["metadata"].merge(self.abundance_matrix.reset_index(), how="outer")
        ### Created on the roll
        self.melted_abundance_dataframe: pd.DataFrame = self.produce_long_abundance_dataframe()
        self.factorize_metadata()

    def _is_indexed(self, df: pd.DataFrame):
        if df.index.name == "smplID":
            return True
        else:
            return False

    def get_bin_list_by_sample(self, sample_id: str, mode: str = "cscope"):
        return self.sample_data[sample_id][mode]["Name"].to_list()

    def get_columns_label(self, sample_id: str, mode: str = "cscope"):
        return self.sample_data[sample_id][mode].columns

    def get_taxonomic_data_by_sample(self, sample_id: str) -> pd.DataFrame:
        return self.taxonomic_data.loc[self.taxonomic_data["mgs"].isin(self.get_bin_list_by_sample(sample_id))]

    def get_bin_list(self, mode: str = "cscope"):
        bin_list = {}
        for sample in self.sample_data.keys():
            bin_list[sample] = self.sample_data[sample][mode].index.to_list()
        return bin_list
        
    def get_factor_len(self):
        return len(self.main_data["metadata"].columns)
    
    def get_metadata_label(self):
        if self._is_indexed(self.main_data["metadata"]):
            return list(self.main_data["metadata"].columns)
        return list(self.main_data["metadata"].columns)[1:]
    
    def get_factor_by_id(self, sample_id, factor):
        query = self.main_data["metadata"].loc[self.main_data['metadata']["smplID"] == sample_id][factor]
        print("get_factor : ", sample_id, factor, "found : ",query)
        return query
    
    def produce_long_abundance_dataframe(self):
        ### Input
        current_df = self.abundance_matrix
        current_df = current_df.reset_index()
        current_df = current_df.melt('smplID',var_name='Compound',value_name='Quantity')
        all_new_columns = {}
        for factor in self.get_metadata_label():
            all_new_columns[factor] = self.truc(current_df['smplID'],factor)
        
        current_df = current_df.assign(**all_new_columns)
        self.main_data["metadata"].reset_index()
        return current_df
    
    def truc(self,serie_id,factor_id):
        metadata = self.main_data["metadata"]
        if not self._is_indexed(metadata):
            metadata.set_index('smplID',inplace=True,drop=True)
        new_col = []
        for value in serie_id:
            new_col.append(str(metadata.at[value,factor_id]))
        return new_col
    
    def factorize_metadata(self):
        for factor in self.get_metadata_label():
            self.main_data["metadata"][factor] = self.main_data["metadata"][factor].astype("category")