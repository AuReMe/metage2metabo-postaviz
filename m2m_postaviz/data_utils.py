import json
import os
import os.path
import tarfile
import time
import sys

import numpy as np
import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci
from scipy import stats
from multiprocessing import Pool
from multiprocessing import cpu_count
from functools import partial
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import squareform,pdist
from statsmodels.stats.multitest import multipletests


def is_valid_dir(dirpath):
    """Return True if directory exists or not
    
    Args:
        dirpath (str): path of directory
    Returns:
        bool: True if dir exists, False otherwise
    """
    if os.path.isdir(dirpath) == True:
        return True
    else:
        return False


def extract_tarfile(tar_file, outdir):
    file = tarfile.open(tar_file, "r:gz")

    file.extractall(outdir, filter="data")
    # if sys.version_info >= (3, 12):
    # else:
    #     tar.extractall(outdir)


def benchmark_decorator(func):
    def wrapper(*args, **kwargs):
        results = list()
        n_repeats = 3
        for i in range(n_repeats):
            time_start = time.perf_counter()
            result = func(*args, **kwargs)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            results.append(time_duration)
            print(f">run {i+1} took {time_duration} seconds")
        avg_duration = sum(results) / n_repeats
        print(f"Took {avg_duration} seconds on average for {func.__name__} function.")
        return result

    return wrapper


def has_only_unique_value(dataframe: pd.DataFrame, input1, input2: str = "None"):
    """
    Return True if df column value contain only one value per input, False otherwise.

    Args:
        dataframe (pd.DataFrame): _description_
        column_value (_type_): _description_
        input1 (_type_): _description_
        input2 (str, optional): _description_. Defaults to "None".
    """
    nb_row = len(dataframe)

    if input2 == "None":
        return True if nb_row == len(dataframe[input1].unique()) else False
        
    else:
        return True if nb_row == len(dataframe[input1].unique()) and nb_row == len(dataframe[input2].unique()) else False


def relative_abundance_calc(sample_data: dict, abundance_path: str, save_path) -> pd.DataFrame:
    """Generate a second main_dataframe with the production based on weight from the abundance matrix.

    Args:
        abundance_matrix (pd.DataFrame): abundance matrix given in input.
        sample_data (dict): Dictionnary of sample's cscopes.

    Raises:
        RuntimeError: If more than one column of type other than INT.

    Returns:
        Dataframe: production dataframe with sample in rows and compounds in column. Weighted by abundance.
    """
    
    if "normalised_abundance_dataframe_postaviz.tsv" in os.listdir(save_path) and "abundance_file.tsv" in os.listdir(save_path):
        print("normalised_abundance_dataframe already exist in save directory.")
        return

    print("Building main dataframe weighted with abundance...")

    abundance_matrix = open_tsv(abundance_path)

    smpl_norm_abundance = []
    smpl_norm_index = []

    # Checking if all column are INT type, if one is not its used as index, if more than 1 raise RuntimeERROR.
    str_filter = abundance_matrix.select_dtypes(include=["string","object","category"])

    if len(str_filter.columns) == 1:
        index_column = str_filter.columns.values[0]
        print(index_column, "column used as index")
        abundance_matrix.set_index(index_column,drop=True,inplace=True)
        abundance_matrix.index.name = "binID"
        
    elif len(str_filter.columns) > 1:
        raise RuntimeError("More than one non-numeric columns in abundance dataframe.")
    
    # Normalisation
    abundance_matrix_normalised = abundance_matrix.apply(lambda x: x / x.sum(), axis=0)

    # For all sample's cscopes, multiply each row (bin's production) by the normalised abundance matrix. 
    for sample in sample_data.keys():

        sample_matrix = sample_data[sample]["cscope"].copy()

        if not is_indexed_by_id(sample_matrix):
            sample_matrix.set_index("smplID", inplace=True)
        
        sample_matrix = sample_matrix.apply(lambda row: row * abundance_matrix_normalised.at[row.name, sample], axis=1)
        sample_matrix = sample_matrix.apply(lambda col: col.to_numpy().sum(), axis=0)
        sample_matrix.name = sample

        smpl_norm_abundance.append(sample_matrix)
        smpl_norm_index.append(str(sample))

    main_dataframe_abundance_weigthed = pd.concat(smpl_norm_abundance,axis=1)
    main_dataframe_abundance_weigthed.fillna(0, inplace=True)
    main_dataframe_abundance_weigthed = main_dataframe_abundance_weigthed.T
    main_dataframe_abundance_weigthed.index.name = "smplID"

    main_dataframe_abundance_weigthed.to_csv(os.path.join(save_path,"normalised_abundance_dataframe_postaviz.tsv"),sep="\t")
    abundance_matrix.to_csv(os.path.join(save_path,"abundance_file.tsv"),sep="\t")
    
    abundance_matrix_normalised.to_csv(os.path.join(save_path,"abundance_file_normalised.tsv"),sep="\t")

    print("Abundance dataframes done and saved.")

    return


def open_added_value(file_name, path):
    for root, dirs, files in os.walk(path):
        if file_name in files:
            added_file = sbml_to_classic(open_json(os.path.join(root, file_name))["addedvalue"])
            return added_file


def open_json(file_path):
    with open(file_path) as file:
        file_data = json.load(file)
    return file_data


def open_tsv(file_name: str, convert_cpd_id: bool = False, rename_columns: bool = False, first_col: str = "smplID"):
    """Open tsv file as a pandas dataframe.

    Args:
        file_name (str): Path of the file
        rename_columns (bool, optional): Rename the first column and decode the metabolites names in sbml format into readable format. Defaults to False.
        first_col (str, optional): Label of the first col if rename_columns is True. Defaults to "smplID".

    Returns:
        Dataframe: Pandas dataframe
    """
    data = pd.read_csv(file_name, sep="\t")
    if rename_columns:
        data.columns.values[0] = first_col
    if convert_cpd_id:
        data.set_index(first_col,inplace=True,drop=True)
        data.columns = sbml_to_classic(data.columns.values)
    return data


def get_scopes(file_name, path) -> pd.DataFrame:
    for root, dirs, files in os.walk(path):
        if file_name in files:
            scope_matrix = open_tsv(os.path.join(root, file_name), convert_cpd_id=True, rename_columns=True)
            return scope_matrix


def is_indexed_by_id(df: pd.DataFrame):
    if df.index.name == "smplID":
        return True
    else:
        return False


def sbml_to_classic(compounds_list):
    uncoded = []
    for coded in compounds_list:
        id, id_type, compart = cfci(coded)
        new_value = str(id)+"["+str(compart)+"]"
        uncoded.append(new_value)
    return uncoded


def contribution_processing(file_opened: dict):
    for key in file_opened.keys():
        for second_key in file_opened[key].keys():
            file_opened[key][second_key] = sbml_to_classic(file_opened[key][second_key])
    return file_opened


def get_contributions(file_name, path):
    for root, dirs, files in os.walk(path):
        if file_name in files:
            contributions_file = open_json(os.path.join(root, file_name))
            contributions_file = contribution_processing(contributions_file)
            return contributions_file


def retrieve_all_sample_data(sample, path):
    """Retrieve iscope, cscope, added_value and contribution_of_microbes files in the path given using os.listdir().

    Args:
        path (str): Directory path

    Returns:
        dict: Return a nested dict object where each key is a dictionnary of a sample. The key of those second layer dict [iscope, cscope, advalue, contribution] give acces to these files.
    """
    sample_directory_path = os.path.join(path, sample)
    if os.path.isdir(sample_directory_path):

        cscope_dataframe = get_scopes("rev_cscope.tsv", sample_directory_path)
        if cscope_dataframe is None:
            return None, sample

    else:
        return None, sample

    return cscope_dataframe, sample


def producers_by_compounds_and_samples_multi(sample_data: dict, save_path):
    """Create and save a dataframe which sum all the compounds produced by each bins in sample cscope for each sample. 

    Args:
        sample_data (dict): Sample's cscope.
        save_path (_type_): Save path given in CLI.

    Raises:
        Exception: Empty sample directory (-d given in CLI).
    """

    if "producers_dataframe_postaviz.tsv" in os.listdir(save_path):
        print("producers dataframe already in save directory.")
        return

    if not bool(sample_data):
        raise Exception("Sample data empty.")

    metadata = open_tsv(os.path.join(save_path,"metadata_dataframe_postaviz.tsv"))

    cpu_available = cpu_count() - 1

    if not type(cpu_available) == int or cpu_available < 1:
        cpu_available = 1

    pool = Pool(cpu_available)
    all_producers = pool.starmap(individual_producers_processing,[(sample_data[sample]["cscope"], sample) for sample in sample_data.keys()])
    pool.close()
    pool.join()

    res = pd.concat(all_producers,axis=1).T
    res.fillna(0,inplace=True)
    res.index.name = "smplID"
    res.reset_index(inplace=True)
    res = res.merge(metadata,'inner',"smplID")

    res.to_csv(os.path.join(save_path,"producers_dataframe_postaviz.tsv"),sep="\t",index=False)

    return


def individual_producers_processing(sample_cscope: pd.DataFrame , sample: str):
    """sums the count of metabolites produced in the sample's cscope.

    Args:
        sample_cscope (pd.DataFrame): Sample cscope dataframe
        sample (str): Sample ID

    Returns:
        pd.Series: Pandas serie with all metabolites columns sum
    """
    serie_value = []
    serie_index = []

    for i in range(len(sample_cscope.columns)):
        serie_index.append(sample_cscope.columns[i])
        serie_value.append(sample_cscope[sample_cscope.columns[i]].to_numpy().sum())

    return pd.Series(serie_value,index=serie_index,name=sample)


def multiprocess_retrieve_data(path):
    """Open all directories given in -d path input. Get all cscopes tsv and load them in memory as pandas
    dataframe. 

    Args:
        path (str): Path of directory

    Returns:
        dict: sample_data dictionnary
    """
    retrieve_data = partial(retrieve_all_sample_data, path=path)

    nb_cpu = cpu_count() - 1
    if not type(nb_cpu) == int or nb_cpu < 1:
        nb_cpu = 1
    pool = Pool(nb_cpu)
    results_list = pool.map(retrieve_data,[sample for sample in os.listdir(path)])
    
    pool.close()
    pool.join()

    all_data = {}
    for df, smpl in results_list:
        if not df is None: 
            all_data[smpl] = {}
            all_data[smpl]["cscope"] = df

    sample_info = {}
    sample_info["bins_list"] = []
    sample_info["bins_count"] = {}
    sample_info["bins_sample_list"] = {}

    for sample in all_data.keys():

        dataframe = all_data[sample]["cscope"]
        all_bins_in_sample = dataframe.index.tolist()

        sample_info["bins_list"] = sample_info["bins_list"] + all_bins_in_sample

        for bin in all_bins_in_sample:
            
            if not bin in sample_info["bins_count"]:
                sample_info["bins_count"][bin] = 0
            
            sample_info["bins_count"][bin] += 1

            if not bin in sample_info["bins_sample_list"]:
                sample_info["bins_sample_list"][bin] = []
            
            sample_info["bins_sample_list"][bin].append(str(sample))

    # Remove duplicate from list
    sample_info["bins_list"] = list(dict.fromkeys(sample_info["bins_list"]))

    return sample_info, all_data


def build_main_dataframe(sample_data: dict, save_path):
    """Create and save the main dataframe. Samples in rows and compounds in columns.
    It takes the compounds production in each samples cscope and return a pandas Series with 1 produced or 0 absent for each compounds.
    Merge all the series returned into a dataframe. 

    Args:
        sample_data (dict): Samples's cscope.
        save_path (_type_): Save path given in CLI.
    """

    if "main_dataframe_postaviz.tsv" in os.listdir(save_path):
        print("main dataframe already in save directory.")
        return

    print("Building main dataframe...")

    all_series = []
    for sample in sample_data.keys():
        current_sample_df = sample_data[sample]["cscope"]
        serie_index = current_sample_df.columns.values
        serie_data = []
        for i in range(len(serie_index)):
            serie_data.append(1)
        all_series.append(pd.Series(data=serie_data,index=serie_index,name=sample))

    results = pd.concat(all_series, axis=1).T
    results.fillna(0,inplace=True)
    results = results.astype(int)
    results.index.name = "smplID"

    results.to_csv(os.path.join(save_path, "main_dataframe_postaviz.tsv"),sep="\t")
    
    print("Main dataframe done and saved.")

    return


def build_dataframes(dir_path, metadata_path: str, abundance_path: str = None, taxonomic_path: str = None, save_path: str = None):
    """
    Main function that build all major dataframes used in shiny. Execpt for the sample's data, all dataframes are saved in parquet format
    in save_path given in CLI.

    Args:
        dir_path (str): Directory path with all samples directory produced by Metage2metabo
        metadata (tsv file): A TSV metadata file.

    Returns:
        None
    """
    if not is_valid_dir(dir_path):
        print(dir_path, "Sample directory path is not a valid directory")
        sys.exit(1)

    if save_path is None:
        sys.exit(1)

    if not is_valid_dir(save_path):
        os.makedirs(save_path)

    metadata_processing(metadata_path, save_path)

    sample_info, sample_data = multiprocess_retrieve_data(dir_path)

    # sample_info = iscope_production(dir_path=dir_path, sample_info_dict=sample_info)

    producers_by_compounds_and_samples_multi(sample_data, save_path) 

    build_main_dataframe(sample_data, save_path)

    relative_abundance_calc(sample_data, abundance_path, save_path)

    taxonomy_processing(taxonomic_path, save_path)

    total_production_by_sample(save_path, abundance_path)

    build_pcoa_dataframe(save_path)

    bin_dataframe_build(sample_info, sample_data, abundance_path, taxonomic_path, save_path)

    # Sample_info JSON

    if os.path.isfile(os.path.join(save_path,"sample_info.json")):
        print(os.path.join(save_path,"sample_info.json"), "directory already exist in save_path.")

    else:
        with open(os.path.join(save_path,"sample_info.json"), 'w') as f:
            json.dump(sample_info, f)

    return 


def metadata_processing(metadata_path, save_path) -> pd.DataFrame:

    if "metadata_dataframe_postaviz.tsv" in os.listdir(save_path):
        return
    else:
        metadata = open_tsv(metadata_path)
        metadata.to_csv(os.path.join(save_path,"metadata_dataframe_postaviz.tsv"),sep="\t", index= True if is_indexed_by_id(metadata) else False)

    return


def list_to_boolean_serie(model_list: list, with_quantity: bool = True):
    results = {}
    value = 1
    for model in model_list:
        if model not in results.keys():
            results[model] = value
        else:
            if with_quantity:
                results[model] += value
    return pd.Series(results)


def total_production_by_sample(save_path, abundance_path: str = None):
    """Create and save the total production dataframe. This dataframe contain all samples in row and all compounds in columns.
    For each samples the compounds produced by each bins is add up to get the estimated total production of compound by samples
    and the number of bins who produced those compounds.

    If the abundance is provided, each production (1) of bins is multiplied by their abundance in their sample which gives an
    estimated production of compounds weighted by the abundance of the bin producer.

    Args:
        save_path (_type_): Save path given in CLI
        abundance_path (str, optional): Abundance file path fiven in CLI. Defaults to None.
    """
    if "total_production_dataframe_postaviz.tsv" in os.listdir(save_path):
        print("total_production_dataframe_postaviz already in save directory")
        return
    
    print("Building total production dataframe...")

    main_dataframe = open_tsv(os.path.join(save_path, "main_dataframe_postaviz.tsv"))
    metadata_dataframe = open_tsv(os.path.join(save_path, "metadata_dataframe_postaviz.tsv"))

    if not is_indexed_by_id(main_dataframe):
        main_dataframe.set_index("smplID",inplace=True,drop=True)

    main_dataframe["Total_production"] = main_dataframe.apply(lambda row: row.to_numpy().sum(), axis=1)
    results = pd.DataFrame(main_dataframe["Total_production"])

    if abundance_path is not None:
        abundance_dataframe = open_tsv(os.path.join(save_path, "normalised_abundance_dataframe_postaviz.tsv"))

        if not is_indexed_by_id(abundance_dataframe):
            abundance_dataframe.set_index("smplID",inplace=True,drop=True)

        abundance_dataframe["Total_abundance_weighted"] = abundance_dataframe.apply(lambda row: row.to_numpy().sum(), axis=1)
        abundance_dataframe = abundance_dataframe["Total_abundance_weighted"]

        results = pd.concat([results,abundance_dataframe], axis=1)

    results.reset_index(inplace=True)
    results = results.merge(metadata_dataframe,'inner','smplID')

    results.to_csv(os.path.join(save_path, "total_production_dataframe_postaviz.tsv"),sep="\t", index=False)

    print("Total production dataframe done and saved.")

    return


def preprocessing_for_statistical_tests(dataframe: pd.DataFrame, y_value, input1, input2 = None, multipletests: bool = False, multipletests_method: str = "bonferroni"):
    """Create dataframe for each y_value in the list, to separate them and use wilcoxon_man_whitney function.
    Concat all results into one dataframe.

    Args:
        dataframe (pd.DataFrame): Dataframe to test.
        y_value (_type_): list of columns labels to separata into several dataframe. Must be of lenght 1 at least. 
        input1 (_type_): First user's input.
        input2 (_type_, optional): Second user's input. Defaults to None.

    Returns:
        Dataframe: Dataframe of statistical test.
    """
    all_results = []

    for y in y_value:

        if input2 is None:
            all_results.append(wilcoxon_man_whitney(dataframe[[y, input1]], y, input1, None, multipletests, multipletests_method))
        else:
            all_results.append(wilcoxon_man_whitney(dataframe[[y, input1, input2]], y , input1, input2, multipletests, multipletests_method))
            
    return pd.concat(all_results)


def wilcoxon_man_whitney(dataframe: pd.DataFrame, y, first_factor: str, second_factor: str = None, multiple_correction: bool = False, correction_method: str = "hs"):
    """ 
    Takes one dataframe with only one value column y and return a dataframe of statistical tests.
    First all sub arrays by the first input then the second input are made and convert to numpy array.
    Then Wilcoxon or Mann Whitney test are run on each pair without doublon.
    If pairs array have the same lenght -> Wilcoxon, if not -> Mann Whitney
    
    Args:
    dataframe (pd.Dataframe): Pandas dataframe
    y (str): Column label containing the values to test.
    first_factor (str): Column label of the first user's input.
    second_factor (str): Column label of the second user's input. Default to None
    

    Returns:
        Dataframe: Dataframe of test's results.
    """
    # Array sorting by the first input and the second input if NOT None.
    sub_dataframes = {}

    for first_factor_array in dataframe[first_factor].unique():

        if second_factor is None:

            sub_dataframes[first_factor_array] = dataframe.loc[dataframe[first_factor] == first_factor_array][y].to_numpy()
            continue

        sub_dataframes[first_factor_array] = {}

        for second_factor_array in dataframe[second_factor].unique():

            sub_dataframes[first_factor_array][second_factor_array] = dataframe.loc[(dataframe[first_factor] == first_factor_array) & (dataframe[second_factor] == second_factor_array)][y].to_numpy()

    # Dataframe's structure declaration, Axis column added if second input is NOT None.
    if second_factor is None:

        results = pd.DataFrame(columns=["Compound", "Factor1", "Sample size1", "Factor2", "Sample size2", "Method", "Statistic", "Pvalue", "Significance"])

    else:

        results = pd.DataFrame(columns=["Compound", "Axis", "Factor1", "Sample size1", "Factor2", "Sample size2", "Method", "Statistic", "Pvalue", "Significance"])

    # Test each pairs avoiding useless duplicates and Array of lenght <= 1.
    for name in sub_dataframes.keys():

        if second_factor is None: # One input selected

            for name2 in sub_dataframes.keys():
                
                if name == name2:
                    continue

                if len(sub_dataframes[name]) < 1 or len(sub_dataframes[name2]) < 1:
                    continue

                if name2 in results["Factor1"].tolist():
                    continue

                if len(sub_dataframes[name]) == len(sub_dataframes[name2]):

                    test_value, pvalue = stats.wilcoxon(sub_dataframes[name], sub_dataframes[name2])
                    test_method = "Wilcoxon"

                else:

                    test_value, pvalue = stats.mannwhitneyu(sub_dataframes[name], sub_dataframes[name2])
                    test_method = "Mann-Whitney"

                if pvalue >= 0.05:
                    symbol = "ns"
                elif pvalue >= 0.01:
                    symbol = "*"
                elif pvalue >= 0.001:
                    symbol = "**"
                else:
                    symbol = "***"

                results.loc[len(results)] = {"Compound": y,
                                            "Factor1": name, "Sample size1": len(sub_dataframes[name]),
                                              "Factor2": name2, "Sample size2": len(sub_dataframes[name2]),
                                                "Statistic": test_value, "Pvalue": pvalue, "Significance": symbol, "Method": test_method}

        else: # Two inputs selected 

            for name2 in sub_dataframes[name].keys():
                
                for name3 in sub_dataframes[name].keys():
                    
                    if name2 == name3:
                        continue

                    if len(sub_dataframes[name][name2]) < 1 or len(sub_dataframes[name][name3]) < 1:
                        continue

                    if len(results.loc[(results["Factor1"] == str(second_factor+": "+str(name3))) & (results["Axis"] == name)]) > 0: # Avoid duplicate
                        continue

                    if len(sub_dataframes[name][name2]) == len(sub_dataframes[name][name3]):

                        test_value, pvalue = stats.wilcoxon(sub_dataframes[name][name2], sub_dataframes[name][name3])
                        test_method = "Wilcoxon"

                    else:

                        test_value, pvalue = stats.mannwhitneyu(sub_dataframes[name][name2], sub_dataframes[name][name3])
                        test_method = "Mann-Whitney"

                    if pvalue >= 0.05:
                        symbol = "ns"
                    elif pvalue >= 0.01:
                        symbol = "*"
                    elif pvalue >= 0.001:
                        symbol = "**"
                    else:
                        symbol = "***"
                    
                    results.loc[len(results)] = {"Compound": y, "Axis": name,
                                                "Factor1": str(second_factor+": "+str(name2)), "Sample size1": len(sub_dataframes[name][name2]),
                                                "Factor2": str(second_factor+": "+str(name3)), "Sample size2": len(sub_dataframes[name][name3]),
                                                    "Statistic": test_value, "Pvalue": pvalue, "Significance": symbol, "Method": test_method}

    if multiple_correction:

        pvals_before_correction = results["Pvalue"].to_numpy()
        reject, pvals_after_correction, _, __ = multipletests(pvals_before_correction, method = correction_method)
        
        results["Pvalue corrected"] = pvals_after_correction
        results["Significance corrected"] = results["Pvalue corrected"].apply(lambda x:get_significance_symbol(x))
        results["Correction method"] = correction_method

    return results


def get_significance_symbol(pval: float) -> str:
    """Return Significance symbol depending on pvalue given.

    Args:
        pval (float): Pvalue of the test

    Returns:
        str: Significance's symbol
    """
    if pval >= 0.05:
        return "ns"
    elif pval >= 0.01:
        return "*"
    elif pval >= 0.001:
        return "**"
    else:
        return "***"


def build_pcoa_dataframe(save_path) -> pd.DataFrame:
    """Comptute Principal Coordinate Analysis from the main_dataframe given in input. Merge with metadata from the smplID column or index.

    Args:
        main_dataframe (pd.DataFrame): Dataframe from which the pcoa will be made.
        metadata (pd.DataFrame): Metadata dataframe (Must have smplID identifer column for the merge to work)

    Returns:
        pd.DataFrame: Pcoa dataframe with sample ID as index, PC1 and PC2 results and all metadata.
    """

    if "pcoa_dataframe_postaviz.tsv" in os.listdir(save_path):
        print("Pcoa dataframe already in save directory.")
        return

    main_dataframe = open_tsv(os.path.join(save_path, "main_dataframe_postaviz.tsv"))

    metadata = open_tsv(os.path.join(save_path, "metadata_dataframe_postaviz.tsv"))

    if not is_indexed_by_id(main_dataframe):
        main_dataframe = main_dataframe.set_index("smplID")

    if is_indexed_by_id(metadata):
        metadata = metadata.reset_index("smplID")

    dmatrix = main_dataframe.to_numpy()
    dist_m = pdist(dmatrix, "jaccard")
    square_df = squareform(dist_m)
    pcoa_result = pcoa(square_df,number_of_dimensions=2)
    coordinate = pcoa_result.samples

    df_pcoa = coordinate[['PC1','PC2']]
    df_pcoa['smplID'] = main_dataframe.index.to_numpy()

    df_pcoa = df_pcoa.merge(metadata, "inner", "smplID")
    df_pcoa.set_index("smplID",inplace=True)

    df_pcoa.to_csv(os.path.join(save_path,"pcoa_dataframe_postaviz.tsv"),sep="\t")

    return


def correlation_test(value_array, factor_array, factor_name, method:str = "pearson"):
    
    if method == "pearson":
        res = stats.pearsonr(value_array, factor_array)
    else:
        res = stats.spearmanr(value_array, factor_array)

    if res.pvalue >= 0.05:
        symbol = "ns"
    elif res.pvalue >= 0.01:
        symbol = "*"
    elif res.pvalue >= 0.001:
        symbol = "**"
    else:
        symbol = "***"

    return pd.DataFrame([[factor_name, len(value_array), res.statistic, res.pvalue, symbol, method]],columns=["Factor", "Sample size", "Statistic", "Pvalue", "Significance", "Method"])


def serie_is_float(ser: pd.Series):

    if np.issubdtype(ser.dtype, np.integer) or np.issubdtype(ser.dtype, np.floating):
        return True
    
    return False


def taxonomy_processing(taxonomy_filepath, save_path):
    """Open taxonomy file and process it if in txt format.

    Args:
        taxonomy_filepath (str): TSV or TXT format

    Raises:
        RuntimeError: Wrong file's format

    Returns:
        pd.DataFrame: Pandas dataframe
    """

    if "taxonomic_dataframe_postaviz.tsv" in os.listdir(save_path):
        print("Taxonomic dataframe already exist in save directory.")
        return

    print("Building taxonomic dataframe...")

    if taxonomy_filepath.endswith(".tsv"):

        df = open_tsv(taxonomy_filepath)
        df = df.rename(columns={df.columns[0]: 'mgs'})
        df.to_csv(os.path.join(save_path,"taxonomic_dataframe_postaviz.tsv"), sep="\t", index=False)

        return

    if not taxonomy_filepath.endswith(".txt"):
        raise RuntimeError("Taxonomy file must be either a txt file or tsv file.")
    
    with open(taxonomy_filepath, 'r') as f:
        lines = f.readlines()

    df = pd.DataFrame(columns=["mgs","kingdom","phylum","class","order","family","genus"])

    del lines[0] # Delete header line

    for row in lines:

        mgs = row.split("\t")[0]
        genus = row.split("\t")[1:]

        k, p, c, o, f, g = genus[0].strip("\n").split(";")

        df.loc[len(df)] = [mgs,k,p,c,o,f,g]


    df.to_csv(os.path.join(save_path,"taxonomic_dataframe_postaviz.tsv"), sep="\t", index=False)

    print("Taxonomic dataframe done and saved.")

    return


def iscope_production(dir_path: str, sample_info_dict: dict):

    indiv_scope_path = "indiv_scopes/rev_iscope.tsv"
    sample_info_dict["iscope"] = {}
    # Takes first of like of sample where bin is present then get iscope production via file rev_iscope.tsv 
    for bin in sample_info_dict["bins_sample_list"].keys():

        if bin in sample_info_dict["iscope"]:
            continue

        sample_used = sample_info_dict["bins_sample_list"][bin][0]
        file_path = os.path.join(os.path.join(dir_path, sample_used), indiv_scope_path)
        df = open_tsv(file_path,True,True)
        bin_row = df.loc[bin]

        sample_info_dict["iscope"][bin] = bin_row.index.tolist()

    return sample_info_dict


def bin_dataframe_build(sample_info: dict, sample_data: dict, abundance_path = None, taxonomy_path = None, savepath = None):
    """Build a large dataframe with all the bins of the different samples as index, the dataframe contain the list of production, abundance, count,
    the metadata and the taxonomic rank associated.

    Args:
        sample_info (dict): _description_
        sample_data (dict): _description_
        metadata (Dataframe): _description_
        abundance_file (Dataframe, optional): _description_. Defaults to None.
        taxonomy_path (Dataframe, optional): _description_. Defaults to None.

    Returns:
        pd.DataFrame: Pandas dataframe
    """
    
    if "bin_dataframe_chunk_1.parquet.gzip" in os.listdir(savepath):
        print("Chunk of bin_dataframe already in save directory")
        return

    print("Building bin dataframe...")

    start = time.time()

    metadata = open_tsv(os.path.join(savepath,"metadata_dataframe_postaviz.tsv"))

    ##### Abundance normalisation, give percentage of abundance of bins in samples.
    if abundance_path is not None:

        abundance_file_normalised = pd.read_csv(os.path.join(savepath,"abundance_file_normalised.tsv"),sep="\t",index_col=0)

    if taxonomy_path is not None:

        taxonomy_file = open_tsv(os.path.join(savepath,"taxonomic_dataframe_postaviz.tsv"))

        ##### Checks if taxonomy file has default index which mean it is not indexed. TEMPORARY until found better way to deal with open/save from -t option OR load taxonomic_df option which return non indexed / indexed df

        if not pd.Index(np.arange(0, len(taxonomy_file))).equals(taxonomy_file.index):
                taxonomy_file = taxonomy_file.reset_index()

        mgs_col_taxonomy = taxonomy_file.columns[0]

    ##### Iterate thought the bins key in sample_info_dict to get list of sample where they are present. 

    sample_list = []
    bin_list = []
    for bin in sample_info["bins_sample_list"].keys():
        sample_list += sample_info["bins_sample_list"][bin]
        bin_list.append(bin)

    #   Delete replicate in list
    sample_unique_list = list(dict.fromkeys(sample_list))

    ##### Create Generator to process data by chunks.
    print("Making chunk of sample list...")

    chunk_generator = chunks(sample_unique_list, 250)

    ##### Loop throught generator

    chunk_index = 0
    for current_chunk in chunk_generator:
        start_chunk = time.time()
        chunk_index += 1
        list_of_dataframe = []
        logs = []

    ##### Loop in sample unique list, get their index (where bins are listed) then select row with isin(bins)

        for sample in current_chunk:

            try:
                df = sample_data[sample]["cscope"]
                rows = df.loc[df.index.isin(bin_list)]
                rows.insert(0 , "smplID", sample)
                rows.index.name = "binID"
                list_of_dataframe.append(rows)

            except:
                logs.append(f'No dataframe named {sample} in sample_data dictionnary')

        results = pd.concat(list_of_dataframe)
        results.fillna(0,inplace=True)
        print(f'Chunk {chunk_index} first concat with {sys.getsizeof(results)/1000000000} Gb memory size')

        s = results.apply(lambda row: get_production_list_from_bin_dataframe(row), axis=1)
        s.name = "Production"
        print(f'Chunk {chunk_index} production serie produced with {sys.getsizeof(s)/1000000000} Gb memory size of production serie')

        count = results.drop("smplID", axis=1).apply(np.sum,axis=1,raw=True)
        count.name = "Count"
        print(f'Chunk {chunk_index} count serie produced with {sys.getsizeof(count)/1000000000} Gb memory size')

        results = pd.concat([results["smplID"],count,s],axis=1)
        del count
        del s

        if abundance_path is not None: # If abundance is provided, multiply each Count column with the relative abundance of the bins in their samples.

            abundance = results.apply(lambda row: abundance_file_normalised.at[row.name,row["smplID"]],axis=1)
            abundance.name = "Abundance"

            final_result = pd.concat([results,abundance], axis=1)
            del results

            final_result["Count_with_abundance"] = final_result["Count"] * final_result["Abundance"]

        else: # Seems quite useless but felt cute might deleted later.
            
            final_result = results
            del results
            
        final_result = final_result.reset_index().merge(metadata, "inner", "smplID")

        if taxonomy_path is not None: # If taxonomy is provided, merge the dataframe with the taxonomic_dataframe.

            final_result = final_result.merge(taxonomy_file, "inner", left_on="binID", right_on=mgs_col_taxonomy)

        ##### Save current chunks into parquet file

        filename = 'bin_dataframe_chunk_'+str(chunk_index)+'.parquet.gzip'
        filepath = os.path.join(savepath,filename)

        if len(final_result) == 0:
            print(f'Chunks {chunk_index} is empty !')

        final_result.to_parquet(filepath, compression='gzip')

        print(f'Chunk {chunk_index} done in {time.time() - start_chunk} with {sys.getsizeof(final_result) / 1000000000} Gb memory size.')
        del final_result


    print("Took: ", time.time() - start, "Before saving")

    if len(logs) != 0:
        with open(os.path.join(savepath,'bin_dataframe_logs.txt'),'w') as log_file:
            for line in logs:
                log_file.write("%s\n" % line)

    return


def get_production_list_from_bin_dataframe(serie: pd.Series) -> list:

    list_of_cpd_produced = []

    for label, value in serie.items():
        if label == "smplID":
            continue
        if value > 0:
            list_of_cpd_produced.append(label)

    return list_of_cpd_produced


def chunks(lst, n):
    """Yield successive n-sized chunks from list."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]