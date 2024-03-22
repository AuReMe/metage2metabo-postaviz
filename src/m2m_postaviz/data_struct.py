import pandas as pd


class DataStorage:
    def __init__(self, main_data_container: dict, sample_data_container: dict, taxonomic_data_container: pd.DataFrame, abundance_data: pd.DataFrame):
        self.main_data = main_data_container
        self.sample_data = sample_data_container
        self.taxonomic_data = taxonomic_data_container
        self.abundance_data = abundance_data

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

    def get_factor_list(self):
        query = list(self.main_data["metadata"].columns)
        # query.insert(0, "None")
        return query

    def get_factor_len(self):
        return len(self.main_data["metadata"].columns)
    
    def get_relevant_metadata(self):
        return list(self.main_data["metadata"].columns)[1:]