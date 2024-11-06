import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
import os
import pickle
import numpy as np
from json import load


class DataStorage:

    ID_VAR = "smplID"
    HAS_TAXONOMIC_DATA : bool = False
    HAS_ABUNDANCE_DATA : bool = False
    SAMPLES_DIRNAME = "all_samples_dataframe_postaviz"
    JSON_FILENAME = 'sample_info.json'
    HDF5_FILENAME = "postaviz_dataframes.h5"
    ABUNDANCE_FILE = 'abundance_file.tsv'
    DF_KEYS = ["metadata_dataframe_postaviz.tsv", "main_dataframe_postaviz.tsv", "normalised_abundance_dataframe_postaviz.tsv",
               "taxonomic_dataframe_postaviz.tsv", "producers_dataframe_postaviz.tsv", "total_production_dataframe_postaviz.tsv",
                "pcoa_dataframe_postaviz.tsv", "abundance_file.tsv"]

    def __init__(self, save_path: str, file_format:str, taxonomy_provided: bool, abundance_provided: bool):

        self.file_format = file_format

        self.HAS_TAXONOMIC_DATA = taxonomy_provided

        self.HAS_ABUNDANCE_DATA = abundance_provided

        if save_path is not None:
            self.output_path = save_path
        else:
            save_path = None
            print("No save path provided, dataframe save option disabled.")


    def open_pickle_file(self, file_path):
        with open(file_path, "rb") as f:
                obj = pickle.load(f)
        return obj


    def open_tsv(self, key: str):
        """Return the dataframe corresponding to the key given as input.

        Args:
            key (str): name of dataframe's file

        Returns:
            pd.Dataframe: Pandas dataframe
        """
        if not key in self.DF_KEYS:
            print(key, "not in keys: ", self.DF_KEYS)
            return
        
        for root, dirname, filename in os.walk(self.output_path):
            if key in filename:
                return pd.read_csv(os.path.join(root,key),sep="\t")


    def open_hdf5(self, key: str):
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
                # print(k)
                if k == key:
                    # print(f'{key} key called.')
                    try:
                        return storage[k]

                    except Exception as e:
                        print(e)
                        return None

        if dataframe is None:
            print(f'No key {key} found in hdf5 file, None value returned.')

        return dataframe


    def from_bin_get_dataframe(self, bin_id, factor) -> pd.DataFrame:
        return self.get_bin_dataframe(self.get_sample_list(bin_id), bin_id, factor)


    def get_bin_dataframe(self, sample_list, bin_id, factor) -> pd.DataFrame:
        cscope_path = os.path.join(self.output_path, self.SAMPLES_DIRNAME)
        res = []
        for sample in sample_list:
            
            tmp_df = pd.read_pickle(os.path.join(cscope_path, sample+"_cscope.pkl"))
            bin_row = tmp_df.loc[tmp_df.index == bin_id]
            bin_row.insert(0, "smplID", sample)
            res.append(bin_row)

        res = pd.concat(res)
        res.fillna(0,inplace=True)
        res.set_index("smplID",inplace=True)

        s = res.apply(lambda row: self.get_production_list_from_bin_dataframe(row), axis=1)
        count = res.apply(np.sum,axis=1)


        res.reset_index(inplace=True)

        abundance_matrix = self.get_raw_abundance_file()
        # metadata = self.get_metadata()
        # metadata = metadata.loc[metadata["smplID"].isin(res["smplID"])]
        # res = res[["smplID","Count","Production"]].merge(metadata, 'inner', "smplID")

        iscope_production = self.get_iscope_production(bin_id)
        iscope = [iscope_production for _ in range(len(res))]

        if abundance_matrix is not None:

            abundance_matrix.columns.values[0] = "binID"
            abundance_matrix.set_index("binID",inplace=True)
            abundance_matrix_normalised = abundance_matrix.apply(lambda x: x / x.sum(), axis=0)
            
            abundance = res.apply(lambda row: abundance_matrix_normalised.at[bin_id,row["smplID"]],axis=1)
            res = pd.concat([s,count,iscope,abundance], axis=1)

        else:

            res = pd.concat([s,count,iscope],axis=1)
            
        print(res)
        return res


    def get_production_list_from_bin_dataframe(self, serie: pd.Series) -> list:

        list_of_cpd_produced = []

        for label, value in serie.items():
            if value > 0:
                list_of_cpd_produced.append(label)

        return list_of_cpd_produced


    def get_iscope_production(self, bin_id) -> list:
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["iscope"][bin_id]


    def get_sample_list(self, bin_id) -> list:
        """Open JSON file in sve path to retrieve dictionnary containing all bins as keys and the list of sample as value. 

        Args:
            bin_id (str): bin_name (shiny call)

        Returns:
            list: List of sample where the bin is present.
        """
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["bins_sample_list"][bin_id]


    def get_bins_list(self) -> list:
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["bins_list"]


    def get_bins_count(self) -> int:
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["bins_count"]


    def get_raw_abundance_file(self):
        return self.open_tsv(key='abundance_file.tsv') if self.HAS_ABUNDANCE_DATA else None


    def get_global_production_dataframe(self) -> pd.DataFrame:
        if self.file_format == "hdf":
            return self.open_hdf5(key='/global_production_dataframe')
        else:
            return self.open_tsv(key="total_production_dataframe_postaviz.tsv")


    def get_metabolite_production_dataframe(self) -> pd.DataFrame:
        if self.file_format == "hdf":
            return self.open_hdf5(key='/metabolite_production_dataframe')
        else:
            return self.open_tsv(key="producers_dataframe_postaviz.tsv")
        

    def get_main_dataframe(self) -> pd.DataFrame:
        if self.file_format == "hdf":
            return self.open_hdf5(key='/main_dataframe')
        else:
            return self.open_tsv(key="main_dataframe_postaviz.tsv")
        

    def get_metadata(self) -> pd.DataFrame:
        if self.file_format == "hdf":
            return self.open_hdf5(key='/metadata')
        else:
            return self.open_tsv(key="metadata_dataframe_postaviz.tsv")


    def get_pcoa_dataframe(self) -> pd.DataFrame:
        if self.file_format == "hdf":
            return self.open_hdf5(key='/pcoa_dataframe')
        else:
            return self.open_tsv(key="pcoa_dataframe_postaviz.tsv")


    # def set_main_metadata(self, new_metadata):
    #     self.metadata = new_metadata

    def get_taxonomic_dataframe(self) -> pd.DataFrame:
        if not self.HAS_TAXONOMIC_DATA:
            return None
        if self.file_format == "hdf":
            return self.open_hdf5(key='/taxonomic_dataframe')
        else:
            return self.open_tsv(key="taxonomic_dataframe_postaviz.tsv")

    def get_normalised_abundance_dataframe(self) -> pd.DataFrame:
        if self.file_format == "hdf":
            return self.open_hdf5(key='/normalised_abundance_dataframe')
        else:
            return self.open_tsv(key="normalised_abundance_dataframe_postaviz.tsv")

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
        

    def get_taxonomy_rank(self) -> list:

        taxonomy = self.get_taxonomic_dataframe()
        if taxonomy is None:
            return ["Taxonomy not provided"]
        
        taxonomy.drop(taxonomy.columns.values[0], axis=1, inplace=True)

        return taxonomy.columns.tolist()
        

    def get_bin_list_from_taxonomic_rank(self, rank, choice):

        taxonomy = self.get_taxonomic_dataframe()
        taxonomy.drop(taxonomy.columns.values[0], axis=1, inplace=True)

        mgs_col_label = taxonomy.columns.values[0]

        return taxonomy.loc[taxonomy[rank] == choice][mgs_col_label].tolist()