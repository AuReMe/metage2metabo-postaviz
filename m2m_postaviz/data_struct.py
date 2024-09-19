import cProfile
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

    def __init__(self, save_path: str, hdf5_file_path: str, taxonomy_provided: bool, abundance_provided: bool):

        self.hdf5_file = hdf5_file_path

        self.HAS_TAXONOMIC_DATA = taxonomy_provided

        self.HAS_ABUNDANCE_DATA = abundance_provided

        self.output_path = save_path

    def open_hdf5(self, key:str):
        """Read HDF5 file containing all relevants dataframes and return the dataframe corresponding to the key given as input.

        Args:
            key (str): Key of specific dataframe.

        Returns:
            Dataframe: Pandas dataframe or None if wrong key.
        """
        dataframe = None
        print(f"Looking for key {key}...")
        with pd.HDFStore(path=self.hdf5_file,mode='r') as storage:
            for k in storage.keys():
                print(k)
                if k == key:
                    print(f'{key} key called.')
                    try:
                        return storage[k]

                    except Exception as e:
                        print(e)
                        return None

        if dataframe is None:
            print(f'No key {key} found in hdf5 file, None value returned.')

        return dataframe

    def get_global_production_dataframe(self) -> pd.DataFrame:
        return self.open_hdf5(key='/global_production_dataframe')

    def get_metabolite_production_dataframe(self) -> pd.DataFrame:
        return self.open_hdf5(key='/metabolite_production_dataframe')

    def get_main_dataframe(self) -> pd.DataFrame:
        return self.open_hdf5(key='/main_dataframe')

    def get_metadata(self) -> pd.DataFrame:
        return self.open_hdf5(key='/metadata')

    def get_pcoa_dataframe(self) -> pd.DataFrame:
        return self.open_hdf5(key='/pcoa_dataframe')
    # def set_main_metadata(self, new_metadata):
    #     self.metadata = new_metadata

    def get_taxonomic_dataframe(self) -> pd.DataFrame:
        return self.open_hdf5(key='/taxonomic_dataframe')

    def get_normalised_abundance_dataframe(self) -> pd.DataFrame:
        return self.open_hdf5(key='/normalised_abundance_dataframe')

    def is_indexed(self, df: pd.DataFrame) -> bool:
        return True if df.index.name == "smplID" else False

    def get_factors(self) -> list:
        return self.get_metadata().columns.tolist()
        
    def get_compound_list(self):
        query = self.get_main_dataframe().columns.tolist()
        if "smplID" in query:
            query.remove("smplID")

        return query

    def get_factor_by_id(self, sample_id, factor):
        query = self.metadata.loc[self.metadata["smplID"] == sample_id][factor]
        print("get_factor : ", sample_id, factor, "found : ", query)
        return query

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
        
        