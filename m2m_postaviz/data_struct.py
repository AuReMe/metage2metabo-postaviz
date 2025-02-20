import os
from json import load
from typing import Optional
from m2m_postaviz.lineage import Lineage
import time
import pandas as pd


class DataStorage:

    ID_VAR = "smplID"
    HAS_TAXONOMIC_DATA : bool = False
    HAS_ABUNDANCE_DATA : bool = False
    USE_METACYC_PADMET : bool = False
    # SAMPLES_DIRNAME = "all_samples_dataframe_postaviz"
    JSON_FILENAME = "sample_info.json"
    ABUNDANCE_FILE = "abundance_file.tsv"
    ALL_FILE_NAMES = ("metadata_dataframe_postaviz.tsv", "main_dataframe_postaviz.tsv", "normalised_abundance_dataframe_postaviz.tsv",
               "taxonomic_dataframe_postaviz.tsv", "producers_dataframe_postaviz.tsv", "producers_iscope_dataframe_postaviz.tsv", "total_production_dataframe_postaviz.tsv",
                "pcoa_dataframe_postaviz.tsv", "abundance_file.tsv", "sample_info.json", "padmet_compounds_category_tree.json")

    BIN_DATAFRAME_PARQUET_FILE = "bin_dataframe.parquet.gzip"
    SRC_DIR = os.path.dirname(os.path.abspath(__file__))
    PADMET_DIR = os.path.join(os.path.dirname(SRC_DIR), "padmet_data")

    def __init__(self, save_path: str):

        if os.path.isdir(save_path) is not True:
            raise FileNotFoundError(f"{save_path} is not a directory.")

        loaded_files = self.load_files(save_path)

        self.HAS_TAXONOMIC_DATA = loaded_files["taxonomic_dataframe_postaviz.tsv"]

        self.HAS_ABUNDANCE_DATA = loaded_files["abundance_file.tsv"]

        self.USE_METACYC_PADMET = loaded_files["padmet_compounds_category_tree.json"]

        if save_path is not None:

            self.output_path = save_path

        print(f"Taxonomy provided : {self.HAS_TAXONOMIC_DATA}\nAbundance provided: {self.HAS_ABUNDANCE_DATA}\nMetacyc database in use: {self.USE_METACYC_PADMET}")


    def open_tsv(self, key: str):
        """Return the dataframe corresponding to the key given as input.

        Args:
            key (str): name of dataframe's file

        Returns:
            pd.Dataframe: Pandas dataframe
        """

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


    # def get_minimal_cpd_dataframe(self, compound) -> pd.DataFrame:

    #     cpd_conditon = [(compound, "=", 1.0)]
    #     col = ["binID", "smplID", compound]

    #     files = []
    #     for i in os.listdir(self.output_path):
    #         if os.path.isfile(os.path.join(self.output_path,i)) and "cpd_" in i:
    #             files.append(i)

    #     if len(files) == 0:
    #         print("No chunk of cpd_dataframe has been found in directory.")
    #         return None

    #     cscope_dataframe = []
    #     iscope_dataframe = []

    #     # Get separate dataframe for each cpd.

    #     for file in files:

    #         df = self.read_parquet_with_pandas(os.path.join(self.output_path, file), col=col, condition=cpd_conditon)

    #         if len(df) == 0:
    #             continue

    #         if "cpd_cscope" in file:

    #             cscope_dataframe.append(df)

    #         if "cpd_iscope" in file:

    #             iscope_dataframe.append(df)

    #     if len(cscope_dataframe) == 0:
    #         return None

    #     dfc = pd.concat(cscope_dataframe).reset_index()

    #     dfi = pd.concat(iscope_dataframe).reset_index()

    #     return dfc, dfi


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


    def get_total_unique_bins_count(self) -> str:
        with open(os.path.join(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return str(len(sample_info["bins_list"]))


    def get_raw_abundance_file(self):
        return self.open_tsv(key="abundance_file.tsv") if self.HAS_ABUNDANCE_DATA else None


    def get_global_production_dataframe(self) -> pd.DataFrame:
        return self.open_tsv(key="total_production_dataframe_postaviz.tsv")


    def get_metabolite_production_dataframe(self, with_metadata = True) -> pd.DataFrame:

        df = self.open_tsv(key="producers_dataframe_postaviz.tsv")

        if with_metadata:

            metadata = self.get_metadata()
            df = df.merge(metadata,"inner","smplID")

        return df


    def get_iscope_metabolite_production_dataframe(self, with_metadata = True) -> pd.DataFrame:

        df = self.open_tsv(key="producers_iscope_dataframe_postaviz.tsv")

        if with_metadata:

            metadata = self.get_metadata()
            df = df.merge(metadata,"inner","smplID")

        return df


    def get_main_dataframe(self) -> pd.DataFrame:
        return self.open_tsv(key="main_dataframe_postaviz.tsv")


    def get_metadata(self) -> pd.DataFrame:
        return self.open_tsv(key="metadata_dataframe_postaviz.tsv")


    def get_pcoa_dataframe(self) -> pd.DataFrame:
        return self.open_tsv(key="pcoa_dataframe_postaviz.tsv")


    def get_list_of_tests(self):
        return ["bonferroni","sidak","holm-sidak","holm","simes-hochberg","hommel","fdr_bh","fdr_by","fdr_tsbh","fdr_tsbky"]

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


    def get_factors(self, remove_smpl_col = False, insert_none = False) -> list:

        result = self.get_metadata().columns.tolist()

        if remove_smpl_col:

             result.remove("smplID")

        if insert_none:

            result.insert(0, "None")

        return result


    def get_sample_list(self) -> list:
        return self.get_main_dataframe()["smplID"].tolist()


    def get_compound_list(self, without_compartment: Optional[bool] = False):
        query = self.get_main_dataframe().columns.tolist()
        if "smplID" in query:
            query.remove("smplID")
        if without_compartment:
            new_query = []
            for cpd in query:
                new_query.append(cpd[:-3])
            return new_query
        return query


    def save_dataframe(self, df_to_save:pd.DataFrame, file_name: str, extension: str = "tsv"):
        """Save the dataframe in input. Check for already saved file and change the name accordingly.

        Args:
            df_to_save (pd.DataFrame): _description_
            file_name (str): _description_
            extension (str, optional): _description_. Defaults to "tsv".

        Returns:
            _type_: _description_
        """
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
        """Associate for each bins in the list a taxonomic rank separated by <;>.

        Args:
            bin_list (list): _description_

        Returns:
            list: _description_
        """
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

            for df_files in self.ALL_FILE_NAMES:

                if df_files in all_files:

                    continue

                if df_files in filenames:

                    all_files[df_files] = True

                else:

                    all_files[df_files] = False

                # print(df_files, "IS \t", all_files[df_files])

        required_files = ["metadata_dataframe_postaviz.tsv", "main_dataframe_postaviz.tsv",
                        "producers_dataframe_postaviz.tsv", "total_production_dataframe_postaviz.tsv",
                        "pcoa_dataframe_postaviz.tsv", "sample_info.json"]

        # Check if necessary files are not True
        for file in required_files:
            if file in all_files and all_files[file] is True:
                continue
            else:
                print(file)
                raise RuntimeError(f"Required {file} is missing when directly loading from directory.")

        return all_files


    def get_added_value_dataframe(self, cpd_input = None, sample_filter_enabled = False, sample_filter_mode = "", sample_filter_value = []):

        cscope_df = self.get_metabolite_production_dataframe(False).sort_values("smplID")

        iscope_df = self.get_iscope_metabolite_production_dataframe(False).sort_values("smplID")

        col_diff = cscope_df.columns.difference(iscope_df.columns)

        col_diff_dict = dict.fromkeys(col_diff, 0.0)

        temp_df = pd.DataFrame(col_diff_dict, index=iscope_df.index)

        iscope_df = pd.concat([iscope_df, temp_df], axis=1)

        if sample_filter_enabled is not False:

            if sample_filter_mode == "Include":

                cscope_df = cscope_df.loc[cscope_df.smplID.isin(sample_filter_value)]
                iscope_df = iscope_df.loc[iscope_df.smplID.isin(sample_filter_value)]

            if sample_filter_mode == "Exclude":

                cscope_df = cscope_df.loc[~cscope_df.smplID.isin(sample_filter_value)]
                iscope_df = iscope_df.loc[~iscope_df.smplID.isin(sample_filter_value)]

        cscope_df.set_index("smplID",inplace=True)
        iscope_df.set_index("smplID",inplace=True)

        cscope_df.sort_index(axis=1,inplace=True)
        iscope_df.sort_index(axis=1,inplace=True)

        if cpd_input is not None:

            cscope_df = cscope_df[[*cpd_input]]
            iscope_df = iscope_df[[*cpd_input]]

        return cscope_df, iscope_df, cscope_df - iscope_df


    def get_cpd_category_tree(self) -> dict:
        with open(os.path.join(self.output_path, "padmet_compounds_category_tree.json")) as fp:
            tree = load(fp)

        return tree


    def get_all_tree_keys(self, tree = None):

        if tree is None:

            tree = self.get_cpd_category_tree()

        lin=Lineage()
        lin.construct_dict(tree,0)
        list_final=list()
        for k,v in lin.level_dict.items():
            list_final.extend(v)

        list_final.insert(0,list(tree.keys())[0])

        return list_final


    def get_sub_tree_recursive(self, data, id, results):
        """Search throught the tree for a match between key and id.
        Return only the part of the tree with the node id as the root.

        Args:
            data (dict): original Tree.
            id (str): ID of the node.
            results (list, optional): List used as transport of results between recursive. Ignore and let it to default. Defaults to [].

        Returns:
            list: list containing the dictionary of the node.
        """
        if len(results) > 0:

            return results

        for key, child in data.items():

            if id == key:

                results.append(data[id])

                return

            else:

                self.get_sub_tree_recursive(child, id, results)

                if len(results) > 0:

                    return


    def find_compounds_from_category(self, data, results):
        """Find and return in a list all the leaf of the tree. each leaf is a compounds
        A compounds has not children, but work need te bo done to be sure that category node
        that do not have any children (not supposed to) will be in the result list.

        Args:
            data (dict): Tree
            results (list, optional): List used as transport of results between recursive. Ignore and let it to default. Defaults to [].

        Returns:
            list: List of childless node found in tree (compounds).
        """

        for key, child in data.items():

            if not bool(child):

                results.append(key)

            else:

                self.find_compounds_from_category(child, results)

        return


    def get_metacyc_category_list(self, tree = None):
        """Return the category list of the metacyc database. By default it return thel list of the category
        of the whole tree. If any sub tree is given it return only the sub category of that tree.

        Args:
            tree (Dict, optional): Sub tree from to get the keys from if None takes the whole tree. Defaults to None.

        Returns:
            List: _description_
        """
        if tree is None:

            tree = self.get_cpd_category_tree()

        res = self.get_all_tree_keys(tree)

        final_res = []

        data_cpd_list = self.get_compound_list(without_compartment=True)

        for key in res:

            cpd_in_category = []

            sub_tree = []

            self.get_sub_tree_recursive(tree, key,sub_tree)

            sub_tree = sub_tree[0]

            self.find_compounds_from_category(sub_tree, cpd_in_category)

            final_cpd_list = [cpd for cpd in data_cpd_list if cpd in cpd_in_category]

            new_key = key+" "+"("+f"{len(final_cpd_list)}"+"/"+f"{len(cpd_in_category)}"+")"

            final_res.append(new_key)
        # start = time.time()

        # shiny_dict_level = {}
        # for key, value in level_dict.items():
        
        #     tmp_dict = {}

        #     for val in value:

        #         tmp_dict[val] = val

        #     key_integer = int(key)

        #     new_key = " "

        #     for i in range(key_integer):

        #         new_key += " "

        #     shiny_dict_level[new_key] = tmp_dict
        # print(f"Took {time.time() - start} sec. --Metacyc_category_list Getter.")
        return final_res
        

