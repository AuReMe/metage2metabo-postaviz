from json import load
from pathlib import Path
from typing import Optional

import pandas as pd
import polars as pl

from m2m_postaviz.lineage import Lineage


class DataStorage:

    ID_VAR = "smplID"
    HAS_TAXONOMIC_DATA : bool = False
    HAS_ABUNDANCE_DATA : bool = False
    USE_METACYC_PADMET : bool = False
    JSON_FILENAME = "sample_info.json"
    ABUNDANCE_FILE = "abundance_file.tsv"
    ALL_FILE_NAMES = ("metadata_dataframe_postaviz.parquet.gzip", "main_dataframe_postaviz.tsv", "normalised_abundance_dataframe_postaviz.tsv",
               "taxonomic_dataframe_postaviz.tsv", "producers_cscope_dataframe.parquet.gzip", "producers_iscope_dataframe.parquet.gzip", "total_production_dataframe_postaviz.tsv",
                "pcoa_dataframe_postaviz.tsv", "abundance_file.tsv", "sample_info.json", "padmet_compounds_category_tree.json")

    def __init__(self, save_path: Path):

        if save_path.is_dir() is not True:
            raise FileNotFoundError(f"{save_path} is not a directory.")

        loaded_files = self.load_files(save_path)

        self.HAS_TAXONOMIC_DATA = loaded_files["taxonomic_dataframe_postaviz.tsv"]

        self.HAS_ABUNDANCE_DATA = loaded_files["abundance_file.tsv"]

        self.USE_METACYC_PADMET = loaded_files["padmet_compounds_category_tree.json"]

        self.current_working_dataframe = {}

        if save_path is not None:

            self.output_path = save_path

        self.raw_data_path = Path(save_path, "raw_dataframes")

        if not self.raw_data_path.is_dir():
            self.raw_data_path.mkdir(parents=True, exist_ok=True)

        print(f"Taxonomy provided : {self.HAS_TAXONOMIC_DATA}\nAbundance provided: {self.HAS_ABUNDANCE_DATA}\nMetacyc database in use: {self.USE_METACYC_PADMET}")

    def keep_working_dataframe(self, tab_id, dataframe, on_RAM_only = False):

        if on_RAM_only:
            self.current_working_dataframe[tab_id] = dataframe
            return

        full_path = Path(self.raw_data_path, tab_id+".tsv")
        
        if isinstance(dataframe, pd.DataFrame):
            dataframe.to_csv(full_path, sep="\t")
        if isinstance(dataframe, pl.DataFrame):
            dataframe.write_csv(full_path, separator="\t")

        print(f"saved {tab_id} at: {full_path}")

    def get_working_dataframe(self, tab_id):

        try:
            current_dataframe = self.current_working_dataframe[tab_id]
        except KeyError:
            print("No current dataframe stored, It may happen when no plot have been made during the session.")
            return

        return current_dataframe

    def open_tsv(self, key: str):
        """Return the dataframe corresponding to the key given as input.

        Args:
            key (str): name of dataframe's file

        Returns:
            pd.Dataframe: Pandas dataframe
        """

        for root, _dirname, filename in self.output_path.walk():
            if key in filename:
                return pl.read_csv(Path(root,key),separator="\t")


    def read_parquet_with_pandas(self, path: Path, col: Optional[list] = None, condition: Optional[list] = None) -> pd.DataFrame:
        """Transfer the column choice and condition as keyword-arguments to the pandas read parquet function.

        Args:
            path (str): path of the parquet file.
            col (Optional[list], optional): Label of the column to open. Defaults to None.
            condition (Optional[list], optional): Nested tuple used to select only rows who matches the conditions. Defaults to None.

        Returns:
            pd.DataFrame: _description_
        """
        kargs = {"path": path}

        if col is not None:

            kargs["columns"] = col

        # Don't use filters parameter due to string_view incompatibility with is_in kernel
        # Read the full dataframe and apply filters in pandas instead
        df = pd.read_parquet(**kargs)

        if condition is not None:
            # Apply filters in pandas to avoid PyArrow filter pushdown issues
            for col_name, op, val in condition:
                if op == "in":
                    df = df[df[col_name].isin(val)]
                elif op == "==":
                    df = df[df[col_name] == val]
                elif op == "!=":
                    df = df[df[col_name] != val]
                elif op == ">":
                    df = df[df[col_name] > val]
                elif op == "<":
                    df = df[df[col_name] < val]
                elif op == ">=":
                    df = df[df[col_name] >= val]
                elif op == "<=":
                    df = df[df[col_name] <= val]

        return df


    def get_bin_dataframe(self, columns = None, condition = None, scope_mode = "cscope") -> pd.DataFrame:
        """Find the bin_dataframe file in the save_path of DataStorage object and read it with the condition given in args.

        Args:
            columns (str, optional): Columns label. Defaults to None.
            condition (Tuple, optional): Tuple of conditions. Defaults to None.

        Returns:
            pd.DataFrame: Resulting bin_dataframe.
        """

        for file in self.output_path.iterdir():
            if file.is_file() and file.name == f"bin_dataframe_{scope_mode}.parquet.gzip":
                return self.read_parquet_with_pandas(file, col=columns, condition=condition)


    def get_iscope_production(self, bin_id) -> list:
        with open(Path(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["iscope"][bin_id]


    def get_bins_list(self) -> list:
        with open(Path(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["bins_list"]


    def get_bins_count(self) -> int:
        with open(Path(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return sample_info["bins_count"]


    def get_total_unique_bins_count(self) -> str:
        with open(Path(self.output_path, self.JSON_FILENAME)) as f:
            sample_info = load(f)

        return str(len(sample_info["bins_list"]))


    def get_raw_abundance_file(self):
        return self.open_tsv(key="abundance_file.tsv") if self.HAS_ABUNDANCE_DATA else None


    def get_global_production_dataframe(self) -> pl.DataFrame:
        return self.open_tsv(key="total_production_dataframe_postaviz.tsv")


    def get_cscope_producers_dataframe(self, col = None, with_metadata = True) -> pl.DataFrame:

        df = pl.read_parquet(Path(self.output_path, "producers_cscope_dataframe.parquet.gzip"), columns=col)

        if with_metadata:

            metadata = self.get_metadata()
            df = df.join(metadata,"smplID",how="left")

        return df


    def get_iscope_producers_dataframe(self, col = None, with_metadata = True) -> pl.DataFrame:

        df = pl.read_parquet(Path(self.output_path, "producers_iscope_dataframe.parquet.gzip"), columns=col)

        if with_metadata:

            metadata = self.get_metadata()
            df = df.join(metadata,"smplID",how="left")

        return df


    def get_main_dataframe(self) -> pl.DataFrame:
        return self.open_tsv(key="main_dataframe_postaviz.tsv")


    def get_metadata(self) -> pl.DataFrame:

        metadata_path = Path(self.output_path, "metadata_dataframe_postaviz.parquet.gzip")

        return pl.read_parquet(metadata_path)


    def get_pcoa_dataframe(self) -> pl.DataFrame:
        return self.open_tsv(key="pcoa_dataframe_postaviz.tsv")


    def get_list_of_tests(self):
        return ["bonferroni","sidak","holm-sidak","holm","simes-hochberg","hommel","fdr_bh","fdr_by","fdr_tsbh","fdr_tsbky"]


    def set_metadata(self, new_metadata):

        metadata_path = Path(self.output_path, "metadata_dataframe_postaviz.parquet.gzip")

        if isinstance(new_metadata, pl.DataFrame):

            new_metadata.write_parquet(metadata_path, compression="gzip")

        if isinstance(new_metadata, pd.DataFrame):

            new_metadata.to_parquet(metadata_path, compression="gzip")


    def get_taxonomic_dataframe(self) -> pl.DataFrame:
        if not self.HAS_TAXONOMIC_DATA:
            return None
        else:
            return self.open_tsv(key="taxonomic_dataframe_postaviz.tsv")


    def get_normalised_abundance_dataframe(self, with_metadata = False) -> pl.DataFrame:
        df = self.open_tsv(key="normalised_abundance_dataframe_postaviz.tsv")

        if with_metadata:

            metadata = self.get_metadata()
            df = df.join(metadata,"smplID",how="left")

        return df


    def get_normalised_iscope_abundance_dataframe(self, with_metadata = False) -> pl.DataFrame:
        df = self.open_tsv(key="normalised_iscope_abundance_dataframe_postaviz.tsv")

        if with_metadata:

            metadata = self.get_metadata()
            df = df.join(metadata,"smplID",how="left")

        return df


    def get_factors(self, remove_smpl_col = False, insert_none = False, with_dtype = False) -> list:

        result = self.get_metadata().columns

        if remove_smpl_col:

             result.remove("smplID")

        if insert_none:

            result.insert(0, "None")

        if with_dtype:

            metadata = self.get_metadata()
            new_name = []
            for col in metadata.columns:

                new_name.append(str(col) + "/" + str(metadata[col].dtype))

            return new_name

        return result


    def get_sample_list(self) -> list:
        return self.get_main_dataframe()["smplID"].to_list()


    def get_compound_list(self, without_compartment: Optional[bool] = False):
        query = self.get_main_dataframe().columns
        if "smplID" in query:
            query.remove("smplID")
        if without_compartment:
            new_query = []
            for cpd in query:
                new_query.append(cpd[:-3])
            return new_query
        return query


    def save_dataframe(self, df_to_save , file_name: str, extension: str = ".tsv"):
        """Save the dataframe in input. Check for already saved file and change the name accordingly.

        Args:
            df_to_save (pd.DataFrame): _description_
            file_name (str): _description_
            extension (str, optional): _description_. Defaults to "tsv".

        Returns:
            _type_: _description_
        """

        final_file_path = Path(self.output_path, file_name).with_suffix(extension)

        if final_file_path.is_file():
            final_file_path = self.check_and_rename(final_file_path)

        if isinstance(df_to_save, pl.DataFrame):
            try:
                df_to_save.write_csv(final_file_path, separator="\t")
                logs = f"Saved in :\n{final_file_path}"
            except Exception as e:
                logs = e
            return logs

        if isinstance(df_to_save, pd.DataFrame):
            try:
                df_to_save.to_csv(final_file_path, sep="\t")
                logs = f"Saved in :\n{final_file_path}"
            except Exception as e:
                logs = e
            return logs


    def save_seaborn_plot(self, sns_obj, file_name):

        if Path(self.output_path, file_name).is_file():
            new_file_name = self.check_and_rename(Path(self.output_path, file_name))
            sns_obj.savefig(new_file_name)
            return f"Filed saved in: {Path(self.output_path, new_file_name)}"

        else:
            sns_obj.savefig(Path(self.output_path, file_name))
            return f"Filed saved in: {Path(self.output_path, file_name)}"


    def check_and_rename(self, file_path: Path, add: int = 0) -> Path:
        original_file_path = file_path
        if add != 0:
            file_name = file_path.stem
            extension = file_path.suffix
            file_name = file_name + "_" + str(add)
            file_path = file_path.with_name(file_name + extension)
        if not file_path.is_file():
            return file_path
        else:
            return self.check_and_rename(original_file_path, add + 1)


    def get_taxonomy_rank(self) -> list:

        taxonomy_col = pl.read_csv(Path(self.output_path, "taxonomic_dataframe_postaviz.tsv"),separator="\t").columns

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

        mgs_column = taxonomic_df.columns[0]

        # taxo_df_indexed = taxonomic_df.set_index(first_col_value)

        res = []

        for bin in bin_list:

            taxonomy = list(taxonomic_df.row(by_predicate=(pl.col(mgs_column) == bin)))

            for i, value in enumerate(taxonomy):

                if type(value) is not str:
                    taxonomy[i] = ""

            new_bin_name = bin + " "

            res.append(new_bin_name + ";".join(taxonomy)) # .values return double list (in case of several lines selected which is not the case here)

        return res


    def get_bin_list_from_taxonomic_rank(self, rank, choice):
        """Return a list of bins corresponding to the taxonomic rank given in input.

        EXAMPLE : taxonomic rank = order, choice = Clostideria.

        Args:
            rank (str): Taxonomic rank
            choice (str): one of the unique choice in taxonomic rank

        Returns:
            list: list of bins in the taxonomic scope
        """
        taxonomy = self.get_taxonomic_dataframe()

        mgs_col_label = taxonomy.columns[0]

        taxonomy = taxonomy.filter(pl.col(rank) == choice)

        return taxonomy.get_column(mgs_col_label).to_list()


    def load_files(self, load_path: Path):
        """Loop through files in save directory and return a dictionnary of True/false for each files.

        If necessary files are not present RaiseRuntimeError

        Args:
            load_path (_type_): _description_

        Raises:
            RuntimeError: If required files are absent.

        Returns:
            dict: _description_
        """
        all_files = {}

        for _root, _dir ,filenames in load_path.walk():

            for df_files in self.ALL_FILE_NAMES:

                if df_files in all_files:

                    continue

                if df_files in filenames:

                    all_files[df_files] = True

                else:

                    all_files[df_files] = False

        required_files = ["metadata_dataframe_postaviz.parquet.gzip", "main_dataframe_postaviz.tsv",
                        "producers_cscope_dataframe.parquet.gzip", "producers_iscope_dataframe.parquet.gzip", "total_production_dataframe_postaviz.tsv",
                        "pcoa_dataframe_postaviz.tsv", "sample_info.json"]

        # Check if necessary files are not True
        for file in required_files:
            if file in all_files and all_files[file] is True:
                continue
            else:
                print(file)
                raise RuntimeError(f"Required file {file} is missing when directly loading from directory.")

        return all_files


    def get_added_value_dataframe(self, cpd_input = None, sample_filter_mode = "", sample_filter_value = None):
        """Return Cscope producers dataframe, Iscope producers dataframe and the difference of these two dataframes.

        Args:
            cpd_input (list, optional): list of compounds of intereset. Defaults to None.
            sample_filter_mode (str, optional): Filter by sample mode. Defaults to "".
            sample_filter_value (list, optional): Filter by sample list of value. Defaults to [].

        Returns:
            pd.DataFrame: Tuple of three (producers) dataframes
        """
        cscope_df = self.get_cscope_producers_dataframe(with_metadata=False).sort("smplID")

        iscope_df = self.get_iscope_producers_dataframe(with_metadata=False).sort("smplID")

        # ccol = cscope_df.columns
        # icol = iscope_df.columns

        # columns_difference = list(set(ccol) - set(icol))

        # for col in columns_difference:

        #     iscope_df = iscope_df.with_columns(pl.lit(0).alias(col)) NOT NEEDED ANYMORE

        if sample_filter_mode != "All" and sample_filter_value is not None:

            if sample_filter_mode == "Include":

                cscope_df = cscope_df.filter(pl.col("smplID").is_in(sample_filter_value))
                iscope_df = iscope_df.filter(pl.col("smplID").is_in(sample_filter_value))

            if sample_filter_mode == "Exclude":

                cscope_df = cscope_df.filter(~pl.col("smplID").is_in(sample_filter_value))
                iscope_df = iscope_df.filter(~pl.col("smplID").is_in(sample_filter_value))

        cscope_df = cscope_df.sort("smplID")
        iscope_df = iscope_df.sort("smplID")

        if cpd_input is not None:

            cscope_df = cscope_df.select(["smplID", *cpd_input])
            iscope_df = iscope_df.select(["smplID", *cpd_input])

        smplid_column = cscope_df["smplID"]

        added_value_dataframe = cscope_df.select(pl.exclude("smplID")) - iscope_df.select(pl.exclude("smplID"))

        added_value_dataframe = added_value_dataframe.with_columns(smplid_column)

        return cscope_df, iscope_df, added_value_dataframe


    def get_cpd_category_tree(self) -> dict:
        with open(Path(self.output_path, "padmet_compounds_category_tree.json")) as fp:
            tree = load(fp)

        return tree


    def get_all_tree_keys(self, tree = None):

        if tree is None:

            tree = self.get_cpd_category_tree()

        # Handle empty tree case
        if not tree or len(tree) == 0:
            return []

        lin=Lineage()
        lin.construct_dict(tree,0)
        list_final= []
        for _k,v in lin.level_dict.items():
            list_final.extend(v)

        list_final.insert(0,list(tree.keys())[0])  # noqa: RUF015

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


    def get_compounds_from_category(self, data, results):
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

                self.get_compounds_from_category(child, results)

        return


    def get_metacyc_category_list(self, tree = None):
        """Return the category list of the metacyc database. By default it return the list of the category
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

            self.get_compounds_from_category(sub_tree, cpd_in_category)

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


    def get_outsider_cpd(self):
        """Return the compounds found in data but doesnt fit in OTHERS category

        Returns:
            Tuple: cpd list / category names
        """
        cpd_in_db = []
        self.get_compounds_from_category(self.get_cpd_category_tree(), cpd_in_db)
        cpd_in_data = self.get_compound_list(True)

        diff = [cpd for cpd in cpd_in_data if cpd not in cpd_in_db]

        return diff, str("Others " + str(len(diff)) + "/" + str(len(cpd_in_db)))


    def get_cpd_label(cpd_index: pd.Series, cpd_list_index):

        return cpd_index[cpd_list_index]
