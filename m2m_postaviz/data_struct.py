import os
from json import load
from typing import Optional

import pandas as pd


class DataStorage:

    ID_VAR = "smplID"
    HAS_TAXONOMIC_DATA : bool = False
    HAS_ABUNDANCE_DATA : bool = False
    SAMPLES_DIRNAME = "all_samples_dataframe_postaviz"
    JSON_FILENAME = "sample_info.json"
    ABUNDANCE_FILE = "abundance_file.tsv"
    DF_KEYS = ("metadata_dataframe_postaviz.tsv", "main_dataframe_postaviz.tsv", "normalised_abundance_dataframe_postaviz.tsv",
               "taxonomic_dataframe_postaviz.tsv", "producers_dataframe_postaviz.tsv", "total_production_dataframe_postaviz.tsv",
                "pcoa_dataframe_postaviz.tsv", "abundance_file.tsv", "sample_info.json")

    BIN_DATAFRAME_PARQUET_FILE = "bin_dataframe.parquet.gzip"

    def __init__(self, save_path: str):

        loaded_files = self.load_files(save_path)

        self.HAS_TAXONOMIC_DATA = loaded_files["taxonomic_dataframe_postaviz.tsv"]

        self.HAS_ABUNDANCE_DATA = loaded_files["abundance_file.tsv"]

        if save_path is not None:
            self.output_path = save_path

        print(f"Taxonomy provided : {self.HAS_TAXONOMIC_DATA}\nAbundance provided: {self.HAS_ABUNDANCE_DATA}")


    def open_tsv(self, key: str):
        """Return the dataframe corresponding to the key given as input.

        Args:
            key (str): name of dataframe's file

        Returns:
            pd.Dataframe: Pandas dataframe
        """
        if key not in self.DF_KEYS:
            print(key, "not in keys: \n", self.DF_KEYS)
            return

        for root, _dirname, filename in os.walk(self.output_path):
            if key in filename:
                return pd.read_csv(os.path.join(root,key),sep="\t")


    def read_parquet_with_pandas(self, path, col: Optional[list] = None, condition: Optional[list] = None) -> pd.DataFrame:

        kargs = {"path": path}

        if col is not None:

            kargs["columns"] = col

        if condition is not None:

            kargs["filters"] = condition

        df = pd.read_parquet(**kargs)

        return df


    def get_bin_dataframe(self, columns = None, condition = None) -> pd.DataFrame:

        files = []
        for i in os.listdir(self.output_path):
            if os.path.isfile(os.path.join(self.output_path,i)) and "bin_dataframe_chunk" in i:
                files.append(i)

        if len(files) == 0:
            print("No chunk of bin_dataframe has been found in directory.")
            return None

        all_df = []

        for file in files:

            df = self.read_parquet_with_pandas(os.path.join(self.output_path, file), col=columns, condition=condition)

            if len(df) == 0:
                continue

            all_df.append(df)

        if len(all_df) == 0:
            return None

        return pd.concat(all_df)


    def get_iscope_production(self, bin_id) -> list:
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["iscope"][bin_id]


    def get_bins_list(self) -> list:
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["bins_list"]


    def get_bins_count(self) -> int:
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["bins_count"]


    def get_raw_abundance_file(self):
        return self.open_tsv(key="abundance_file.tsv") if self.HAS_ABUNDANCE_DATA else None


    def get_global_production_dataframe(self) -> pd.DataFrame:
        return self.open_tsv(key="total_production_dataframe_postaviz.tsv")


    def get_metabolite_production_dataframe(self) -> pd.DataFrame:
        return self.open_tsv(key="producers_dataframe_postaviz.tsv")


    def get_main_dataframe(self) -> pd.DataFrame:
        return self.open_tsv(key="main_dataframe_postaviz.tsv")


    def get_metadata(self) -> pd.DataFrame:
        return self.open_tsv(key="metadata_dataframe_postaviz.tsv")


    def get_pcoa_dataframe(self) -> pd.DataFrame:
        return self.open_tsv(key="pcoa_dataframe_postaviz.tsv")


    # def set_main_metadata(self, new_metadata):
    #     self.metadata = new_metadata

    def get_taxonomic_dataframe(self) -> pd.DataFrame:
        if not self.HAS_TAXONOMIC_DATA:
            return None
        else:
            return self.open_tsv(key="taxonomic_dataframe_postaviz.tsv")

    def get_normalised_abundance_dataframe(self) -> pd.DataFrame:
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

        taxonomy_col = pd.read_csv(os.path.join(self.output_path, "taxonomic_dataframe_postaviz.tsv"),sep="\t").columns.tolist()

        if taxonomy_col is None:
            return ["Taxonomy not provided"]

        return taxonomy_col


    def associate_bin_taxonomy(self, bin_list:list) -> list:

        taxonomic_df = self.get_taxonomic_dataframe()

        first_col_value = taxonomic_df.columns.values[0]

        taxo_df_indexed = taxonomic_df.set_index(first_col_value)

        res = []

        for bin in bin_list:

            taxonomy = taxo_df_indexed.loc[taxo_df_indexed.index == bin].values[0].tolist()

            for i, value in enumerate(taxonomy):

                if type(value) is not str:
                    taxonomy[i] = ""

            new_bin_name = bin + " "

            res.append(new_bin_name + ";".join(taxonomy)) # .values return double list (in case of several lines selected which is not the case here)

        return res


    def get_bin_list_from_taxonomic_rank(self, rank, choice):

        taxonomy = self.get_taxonomic_dataframe()

        mgs_col_label = taxonomy.columns.values[0]

        return taxonomy.loc[taxonomy[rank] == choice][mgs_col_label].tolist()


    def load_files(self, load_path):

        all_files = {}

        for _root, _dir ,filenames in os.walk(load_path):

            for df_files in self.DF_KEYS:

                if df_files in all_files:

                    continue

                if df_files in filenames:

                    all_files[df_files] = True

                else:

                    all_files[df_files] = False

                print(df_files, "IS \t", all_files[df_files])

        required_files = ["metadata_dataframe_postaviz.tsv", "main_dataframe_postaviz.tsv", "producers_dataframe_postaviz.tsv", "total_production_dataframe_postaviz.tsv", "pcoa_dataframe_postaviz.tsv", "sample_info.json"]

        # Check if necessary files arent' True
        for file in required_files:
            if file in all_files and all_files[file] is True:
                continue
            else:
                print(file)
                raise RuntimeError("Required files are missing when directly loading from directory.")

        return all_files
