import json
import os
import os.path
import time
import tarfile

import pandas as pd
from padmet.utils.sbmlPlugin import convert_from_coded_id as cfci
from scipy import stats
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
# from m2m_postaviz.time_decorator import timeit
# import cProfile
# import threading



def weird_way_to_do_it(id_value: str, metadata_col: str, metadata: pd.DataFrame):
    if is_indexed_by_id(metadata):
        metadata.reset_index(inplace=True)
    result = metadata.loc[metadata['smplID'] == id_value][metadata_col].values[0]
    return result


def run_pcoa(main_df: pd.DataFrame, metadata: pd.DataFrame, distance_method: str = "jaccard"):
    """Calculate Principal Coordinate Analysis with dataframe given in first arg.
    Use metadata's drataframe as second argument to return the full ordination result plus
    all metadata column inserted along Ordination.samples dataframe.
    Ready to be plotted.

    Args:
        main_df (pd.DataFrame): Main dataframe of compound production
        metadata (pd.DataFrame): Metadata's dataframe

    Returns:
        _type_: Ordination results object from skbio's package.
    """
    df = main_df.copy()
    # Need the matrix version for distance calculation.
    if not is_indexed_by_id(main_df):
        main_df.set_index("smplID", inplace=True)

    # Add metadata columns with the accurate value in temporary dataframe.
    for col in metadata.columns:
        if col == "smplID":
            continue
        df[col] = df["smplID"].apply(lambda row: weird_way_to_do_it(row,col,metadata))

    # Normalisation
    main_df = main_df.apply(lambda x: x / x.sum(), axis=0)
    # Calculate distance matrix with Bray-Curtis method.
    dist_m = pdist(main_df, distance_method)
    # Transform distance matrix into squareform.
    squaref_m = squareform(dist_m)
    # Run the PCOA with the newly generated distance matrix.
    pcoa_results = pcoa(squaref_m, number_of_dimensions=main_df.shape[0]-1)

    # Verify if PCOA_results's samples dataframe is aligned with df dataframe.
    df = df.set_index(pcoa_results.samples.index)
    # Put each metadata column in the pcoa results's dataframe for PLOT.
    for col in metadata.columns:
        if col == "smplID":
            continue
        pcoa_results.samples[col] = df[col].astype('category')

    if is_indexed_by_id(main_df):
        main_df.reset_index(inplace=True)

    return pcoa_results

def list_to_boolean_serie(model_list: list, with_quantity: bool = True):
    results = {}
    for model in model_list:
        if not model in results.keys():
            value = 1
            results[model] = value
        elif with_quantity:
            results[model] += value
    return pd.Series(results)


def intest_taxonomy_matrix_build(binlist_by_id: dict, taxonomic_df: pd.DataFrame):
    # temporary solution, should be species if available.
    rank = 'Genus'
    id_col = 'mgs'
    all_series = {}
    # for each sample get the row of the taxo_df.
    for sample in binlist_by_id.keys():
        res = taxonomic_df.loc[taxonomic_df[id_col].isin(binlist_by_id[sample])][rank]
        all_series[sample] = list_to_boolean_serie(res)
    # Concatenation of the series into a dataframe.
    matrix = pd.DataFrame(all_series)
    matrix.fillna(0,inplace=True)
    # matrix = matrix.T
    matrix.reset_index(inplace=True)

    # print(matrix)
    return matrix


def wid_to_long_format(df: pd.DataFrame):
    return df.melt('smplID',var_name='Compound',value_name='Quantity')


def extract_tarfile(tar_file, outdir):
    file = tarfile.open(tar_file, 'r:gz')

    file.extractall(outdir, filter='data')
    # if sys.version_info >= (3, 12):
    # else:
    #     tar.extractall(outdir)


def multiply_production_abundance(df_row: pd.Series, abundance_matrix: pd.DataFrame,sample_id):
    df_row = df_row.astype(float)
    abundance_value = abundance_matrix.at[df_row.name,sample_id]
    df_row *= abundance_value
    return df_row


def generate_normalized_stoichiometric_matrix(binary_matrix: pd.DataFrame, abundance_matrix: pd.DataFrame, sample_id: str):
    """
    Produce a stoichiometric matrix from a binary matrix (presence or absence) and an abundance matrix.

    Args:
        binary_matrix (pandas.DataFrame): Binary dataframe containing the presence or absence of compounds.
        abundance_matrix (pandas.DataFrame): Dataframe containing the abundance of each bin in sample.
        sample_id (str): Name of the sample.

    Returns:
        pandas Dataframe: Dataframe with the theorical quantity of compounds produced by each bin. 
    """
    if not is_indexed_by_id(binary_matrix):
        binary_matrix.set_index("smplID",inplace=True)
    # Normalisation
    normalized_abundance = abundance_matrix.apply(lambda x: x / x.sum(), axis=0)

    # Multiplication par abondance
    binary_matrix = binary_matrix.apply(lambda row: multiply_production_abundance(row, normalized_abundance,sample_id),axis=1) 
    # Retourne une matrice type CSCOPE avec un abondance relative a la quantité de bin.
    return binary_matrix


def relative_abundance_calc(abundance_file_path: str, sample_data: dict):
    """Use the abundance matrix given in input with the dataframe of the sample's production in community to create an 2 abundance dataframe.
    1 with normalised bin quantity with x / x.sum(), the other without any transformation.

    Args:
        abundance_file_path (str): path to the abundance file.
        sample_data (dict): dictionnary of the cscope, iscope of each sample.

    Returns:
        Tuple: (Normalised abundance dataframe, Abundance dataframe)
    """
    abundance_matrix = open_tsv(abundance_file_path)
    smpl_norm_abundance = []
    smpl_norm_index = []
    smpl_abundance = []
    smpl_index = []

    for sample in sample_data.keys():
        smpl_abundance.append(sum_abundance_table(generate_stoichiometric_matrix(sample_data[sample]["cscope"], abundance_matrix, sample),sample))
        smpl_index.append(str(sample))

    for sample in sample_data.keys():
        smpl_norm_abundance.append(sum_abundance_table(generate_normalized_stoichiometric_matrix(sample_data[sample]["cscope"], abundance_matrix, sample),sample))
        smpl_norm_index.append(str(sample))

    normalized_abundance = pd.concat(smpl_norm_abundance, join="outer", ignore_index=True)
    normalized_abundance.fillna(0,inplace=True)
    normalized_abundance.insert(0,"smplID",smpl_norm_index)
    normalized_abundance.set_index("smplID",inplace=True,drop=True)

    vanilla_abundance = pd.concat(smpl_abundance, join="outer", ignore_index=True)
    vanilla_abundance.fillna(0,inplace=True)
    vanilla_abundance.insert(0,"smplID",smpl_index)
    vanilla_abundance.set_index("smplID",inplace=True,drop=True)

    return normalized_abundance, vanilla_abundance


def sum_abundance_table(abundance_table: pd.DataFrame, sample_id: str):
    # Prend la nouvelle matrice d'abondance du sample crée
    new_dataframe =  {}
    # Flip avec métabolites en index et bin en column
    abundance_table = abundance_table.T
    for index,row in abundance_table.iterrows():
        # Pour chaque métabolites, calcul le total crée par l'ensemble des bin de l'échantillon
        new_dataframe[index] = row.values.sum()
    return pd.DataFrame(new_dataframe, index=[sample_id])


def taxonomic_overview(list_bin_id, taxonomic_dataframe, metadata, mode: str = "cscope",with_count: bool = True):
    current_selection = []
    current_rank = taxonomic_dataframe.columns[-1]

    for sample in list_bin_id.keys():
        tmp_df = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(list_bin_id[sample])]
        tmp_df = tmp_df[["mgs", current_rank]]
        if with_count:
            tmp_df = tmp_df.groupby([current_rank]).count()
            tmp_df = tmp_df.reset_index()
            tmp_df.columns = [current_rank, "Count"]
        value, label = get_metadata(sample, metadata)
        for i in range(len(label)):
            tmp_df.insert(0, label[i], value[i])
        tmp_df.insert(0, 'nb_taxon', tmp_df[current_rank].nunique())

        current_selection.append(tmp_df)

    return pd.concat(current_selection, join="outer", ignore_index=True)


def get_metadata(sample_id: str, metadata: pd.DataFrame):
    return tuple(
        [metadata.loc[metadata["smplID"] == sample_id].values.tolist()[0], metadata.loc[metadata["smplID"] == sample_id].columns.to_list()]
    )


def taxonomy_groupby(
    metadata: pd.DataFrame,
    current_sample: str,
    bin_id_by_sample: dict,
    taxonomic_dataframe: pd.DataFrame,
    target_rank: str = "s",
    taxonomic_choice: list = [],
):
    """Generate a taxonomic count dataframe from a sample id, his dataframe and the rank choosen. Return only the selected taxonomic choice.

    Args:
        current_sample (str): Sample's id
        taxonomic_dataframe (pd.DataFrame): The taxonomic dataframe
        target_rank (str, optional): Selected rank in shiny's input. Defaults to "Genus".
        taxonomic_choice (list, optional): The list of taxonomic selection from shiny's input to keep in the returned Dataframe. Defaults to [].

    Returns:
        Dataframe: Pandas dataframe containing the selected taxonomic choice and rank.
    """
    df = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(bin_id_by_sample[current_sample])]
    for choice in taxonomic_choice:
        if not choice in df[target_rank].unique():
            taxonomic_choice.remove(choice)
            print(choice, " Removed.")

    df = df[["mgs", target_rank]]
    df = df.groupby([target_rank]).count()
    df = df.reset_index()
    df.columns = [target_rank, "Count"]
    df = df.loc[df[target_rank].isin(taxonomic_choice)]
    value, label = get_metadata(current_sample, metadata)
    for i in range(len(label)):
        df.insert(0, label[i], value[i])
    return df


def taxonomic_dataframe_from_input(
    taxonomic_rank: str, bin_id_by_sample: dict, taxonomic_choice: list, taxonomic_dataframe: pd.DataFrame, metadata: pd.DataFrame
):
    results = []
    taxonomic_choice = list(taxonomic_choice)
    if len(taxonomic_choice) == 0:
        print("The taxonomic choice list is empty")
        return
    for sample in bin_id_by_sample.keys():
        results.append(taxonomy_groupby(metadata, sample, bin_id_by_sample, taxonomic_dataframe, taxonomic_rank, taxonomic_choice))
    final_results = pd.concat(results, join="outer", ignore_index=True)
    return final_results


def wilcoxon_test(group1, group2):
    try:
        res = stats.wilcoxon(group1, group2)
    except Exception as e:
        res = e
    return res


def student_test(group1, group2):
    try:
        res = stats.ttest_ind(group1, group2)
    except Exception as e:
        res = e
    return res


def get_added_value_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["advalue"])


def get_individual_production_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["iscope"].columns)


def get_total_production_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["cscope"].columns)


def get_taxonomy_size(sample_data: pd.DataFrame, taxonomic_dataframe: pd.DataFrame, only_metabolic_model_size: bool = False):
    """Return the numbers of different species in one sample.

    Args:
        sample_data (pd.DataFrame): dataframe of one sample.
        taxonomic_dataframe (pd.DataFrame): taxonomic dataframe of all samples.
        only_metabolic_model_size (bool, optional): Return only the number of individual in the sample. Defaults to False.

    Returns:
        int: Number of individual or number of different species.
    """
    if only_metabolic_model_size:
        taxonomy_size = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(sample_data["smplID"])]
        return len(taxonomy_size)
    taxonomy_size = taxonomic_dataframe.loc[taxonomic_dataframe["mgs"].isin(sample_data["smplID"])]
    taxonomy_size = taxonomy_size[["mgs", "Genus"]]
    taxonomy_size = taxonomy_size.groupby(["Genus"]).count()
    taxonomy_size = taxonomy_size.reset_index()
    return len(taxonomy_size)


def get_metabolic_info(sample_data: dict, metadataframe: pd.DataFrame, taxonomic_dataframe: pd.DataFrame):
    """Really ugly and probably useless need rework.
    Return a dataframe with global information of all the sample in input.

    Args:
        sample_data (dict): _description_
        metadataframe (pd.DataFrame): _description_
        taxonomic_dataframe (pd.DataFrame): _description_

    Returns:
        _type_: _description_
    """
    tot_size = []
    ind_size = []
    ad_size = []
    taxo_size = []
    model_size = []
    for sample in metadataframe["smplID"]:
        tot_size.append(get_total_production_size(sample, sample_data))
        ind_size.append(get_individual_production_size(sample, sample_data))
        ad_size.append(get_added_value_size(sample, sample_data))
        taxo_size.append(get_taxonomy_size(sample_data[sample]["iscope"], taxonomic_dataframe))
        model_size.append(get_taxonomy_size(sample_data[sample]["iscope"], taxonomic_dataframe, only_metabolic_model_size=True))
    new_df = metadataframe
    new_df["prod_community"] = tot_size
    new_df["prod_individual"] = ind_size
    new_df["added_value_total"] = ad_size
    new_df["Numbers of models"] = model_size
    new_df["Numbers of species"] = taxo_size
    return new_df


def get_metabolic_model_prod_size(sample: str, metadataframe: dict):
    return len(metadataframe[sample]["cscope"].columns)


def remove_metadata(df: pd.DataFrame) -> pd.DataFrame:
    new_df = new_df.drop("Time_rel", axis=1)
    # df = df.drop("Name", axis=1)
    new_df = new_df.drop("Antibiotics", axis=1)
    new_df = new_df.set_index("smplID")
    return new_df


def add_row(df: pd.DataFrame, row: list):
    df.loc[len(df)] = row


def deal_with_duplicated_row(df: pd.DataFrame, dropcol=None):
    if not dropcol == None:
        df = df.drop(dropcol, axis=1)
    df = pd.get_dummies(df, prefix="", prefix_sep="", columns=["category"])
    df = df.groupby(["compound_id"]).sum()
    return df


def get_threshold_value(df: pd.DataFrame, threshold: int = 0, transpose: bool = False) -> dict:
    """Take a matrix with percentage as value instead of a count and return only value above threshold.
    value must be row indexed and sample in columns.

    Args:
        df (pd.DataFrame): Dataframe to extract values from.
        threshold (int, optional): The threshold used to extract value (by default all percentage above but not equal 0 are exctrated). Defaults to 0.
        transpose (bool, optional): True if the dataframe needs to be transpose (when value in columns and sample in rows). Defaults to False.

    Returns:
        dict: A dictionary, for which key correspond to a sample and value in a nested list of value above threshold and their index. return{sample1 : [[value1,value2],[index_value1,index_value2]]}
    """
    if transpose:
        df = df.T
    results = {}
    for sample in df:
        res = df.loc[df[sample] > 1]
        res = res[sample]
        results[sample] = res
    return results


def get_percent_value(df: pd.DataFrame, transpose: bool = False) -> pd.DataFrame:
    """Calculate the total of each columns and replace value by their percentage.

    Args:
        df (pd.DataFrame): Count matrix as dataframe value must be row indexed and sample in columns.
        transpose (bool, optional): Transpose the matrix. Defaults to False.

    Returns:
        pd.DataFrame: ERROR 404
    """
    if transpose:
        df = df.T
    total = df.sum()
    percentage = (df / total) * 100
    return percentage


def get_files(file_name: str, path: str, with_directory_name: bool = True):
    """Retrieve files from the name within the directory path given in CLI.

    Args:
        file_name (str): Name of the file
        path (str): directory path (-d from cli) from which os.walk begin the search.
        with_directory_name (bool, optional): If True, also return the directory name containing the file. Defaults to True.

    Returns:
        _type_: _description_
    """
    result = []
    for root, dirs, files in os.walk(path):
        if file_name in files:
            if with_directory_name:
                result.append(
                    [os.path.join(root, file_name), os.path.basename(os.path.dirname(os.path.dirname(os.path.join(root, file_name))))]
                )
            else:
                result.append(os.path.join(root, file_name))
    return result


def open_added_value(file_name, path):
    for root, dirs, files in os.walk(path):
        if file_name in files:
            added_file = sbml_to_classic(open_json(os.path.join(root, file_name))["addedvalue"])
            return added_file


def open_json(file_path):
    with open(file_path) as file:
        file_data = json.load(file)
    return file_data


def open_tsv(file_name: str, rename_columns: bool = False, first_col: str = "smplID"):
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
        data.columns.values[1:] = sbml_to_classic(data.columns.values[1:])
    return data


def get_scopes(file_name, path):
    for root, dirs, files in os.walk(path):
        if file_name in files:
            iscope_dataframe = open_tsv(os.path.join(root, file_name), rename_columns=True)
            return iscope_dataframe


def convert_to_dict(list_of_produced_cpd):
    decoded_cpd_list = sbml_to_classic(list_of_produced_cpd)
    data = {}
    for i in decoded_cpd_list:
        value = [1]
        data[i] = value
    return data


def is_indexed_by_id(df: pd.DataFrame):
    if df.index.name == "smplID":
        return True
    else:
        return False


def get_column_size(df: pd.DataFrame):
    size = len(df.columns)
    return size


def get_row_size(df: pd.DataFrame):
    size = len(df)
    return size


def get_columns_index(df: pd.DataFrame, key_list):
    index_list = [0]
    for k in key_list:
        index_list.append(df.columns.get_loc(k))
    return index_list


def convert_to_dataframe(all_file):
    all_dataframe = []
    for file in all_file:
        current_file = open_json(file)["com_scope"]
        all_dataframe.append(pd.DataFrame.from_dict(convert_to_dict(current_file)))
    return all_dataframe


def sbml_to_classic(list_of_metabolites):
    uncoded = []
    for coded in list_of_metabolites:
        id, id_type, compart = cfci(coded)
        uncoded.append(id)
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


def retrieve_all_sample_data(path):
    """Retrieve iscope, cscope, added_value and contribution_of_microbes files in the path given using os.listdir().

    Args:
        path (str): Directory path given in CLI.

    Returns:
        dict: Return a nested dict object where each key is a dictionnary of a sample. The key of those second layer dict [iscope, cscope, advalue, contribution] give acces to these files. 
    """
    data_by_sample = {}
    for sample in os.listdir(path):
        if os.path.isdir(os.path.join(path, sample)):
            data_by_sample[sample] = {}
            # data_by_sample[sample]["iscope"] = get_scopes("rev_iscope.tsv", os.path.join(path, sample))
            data_by_sample[sample]["cscope"] = get_scopes("rev_cscope.tsv", os.path.join(path, sample))
            # data_by_sample[sample]["advalue"] = open_added_value("addedvalue.json", os.path.join(path, sample))
            # data_by_sample[sample]["contribution"] = get_contributions("contributions_of_microbes.json", os.path.join(path, sample))
    return data_by_sample


# def retrieve_sample_data(sample_directory_path, sample_name, sample_dictionnary: dict):
#     sample_dictionnary[sample_name] = {}
#     sample_dictionnary[sample_name]["iscope"] = get_scopes("rev_iscope.tsv", sample_directory_path)
#     sample_dictionnary[sample_name]["cscope"] = get_scopes("rev_cscope.tsv", sample_directory_path)
#     sample_dictionnary[sample_name]["advalue"] = open_added_value("addedvalue.json", sample_directory_path)
#     sample_dictionnary[sample_name]["contribution"] = get_contributions("contributions_of_microbes.json", sample_directory_path)
#     return


def merge_metadata_with_df(main_dataframe, metadata):
    return pd.merge(metadata, main_dataframe, how="left")


def merge_df(left_df, right_df, how: str = "left"):
    data = left_df.iloc[:, 0]
    filter = right_df.loc[right_df["mgs"].isin(data)]
    return filter


def multiply_production_abundance(row: pd.Series, abundance_matrix: pd.DataFrame, sample_id):
    row = row.astype(float)
    abundance_value = abundance_matrix.at[row.name, sample_id]
    row *= abundance_value
    return row


def generate_stoichiometric_matrix(binary_matrix: pd.DataFrame, abundance_matrix: pd.DataFrame, sample_id: str):
    """Produce a stoichiometric matrix from a binary matrix (presence or absence) and an abundance matrix.

    Args:
        binary_matrix (pandas.DataFrame): Binary dataframe containing the presence or absence of compounds.
        abundance_matrix (pandas.DataFrame): Dataframe containing the abundance of each bin in sample.
        sample_id (str): Name of the sample.

    Returns:
        pandas Dataframe: Dataframe with the theorical quantity of compounds produced by each bin.
    """
    if not is_indexed_by_id(binary_matrix):
        binary_matrix.set_index("smplID", inplace=True)

    binary_matrix = binary_matrix.apply(lambda row: multiply_production_abundance(row, abundance_matrix, sample_id), axis=1)
    return binary_matrix


# @timeit(repeat=3,number=10)
def build_df(dir_path, metadata, abundance_path):
    """
    Extract community scopes present in directory from CLI then build a single dataframe from the metabolites produced by each comm_scopes.

    Args:
        dir_path (str): Directory path containing comm scopes
        metadata (tsv file): tsv file containing the metadata of the scopes. The number of row must be equal to the number of comm_scopes given in dir_path.

    Returns:
        global_data: dict
        sample_data: dict
        abundance_data: pandas dataframe
    """
    data_by_sample = {}
    sample_data = retrieve_all_sample_data(dir_path)

    ### Multi-Threading for building data, slower than vanilla (for NOW).

    # all_task = []
    # data_by_sample = {}
    # for sample_directory in os.listdir(dir_path):
    #     if os.path.isdir(os.path.join(dir_path, sample_directory)):
    #         task = threading.Thread(target=retrieve_sample_data,args=(os.path.join(dir_path, sample_directory), sample_directory, data_by_sample))
    #         task.start()
    #         all_task.append(task)

    # for task in all_task:
    #     task.join()

    ### Main dataframe and metadata dataframe. ###

    global_data = {}
    dir_files = get_files("comm_scopes.json", dir_path)
    file_list = []
    dir_list = []
    for file, dir in dir_files:
        file_list.append(file)
        dir_list.append(dir)

    all_df = convert_to_dataframe(file_list)

    main_df = pd.concat(all_df, join="outer", ignore_index=True)
    main_df = main_df.fillna(0)
    main_df = main_df.astype(int)
    main_df.insert(0, "smplID", dir_list)

    global_data["metadata"] = open_tsv(metadata)

    global_data["main_dataframe"] = main_df

    norm_abundance_data, abundance_data = relative_abundance_calc(abundance_path, sample_data)

    return global_data, sample_data, norm_abundance_data, abundance_data


def build_test_data(test_dir_path):
    """Load, open and return the test data to be used by the DataStorage object.

    Args:
        test_dir_path (_type_): Path to the archive

    Returns:
        dict, dict, pd.dataframe: 2 dict, 1 pandas dataframe
    """
    global_data = {}
    sample_data = {}
    for file in os.listdir(test_dir_path):
        filename = os.fsdecode(file)
        if not filename.endswith('.tsv'):
            continue
        if not filename.startswith("sample"):
            if filename == 'main_table.tsv':
                global_data["main_dataframe"] = open_tsv(os.path.join(test_dir_path,filename))
            if filename == 'metadata_table.tsv':
                global_data["metadata"] = open_tsv(os.path.join(test_dir_path,filename))
            if filename == 'abundance_table.tsv':
                abundance_table = open_tsv(os.path.join(test_dir_path,filename))
        else:
            current_file = os.path.splitext(filename)[0]
            current_file = current_file.split("_")
            sample_data[current_file[2]] = {}

    for sample_file in os.listdir(test_dir_path):
        filename = os.fsdecode(sample_file)
        if filename.startswith('sample'):
            current_file = os.path.splitext(filename)[0]
            current_file = current_file.split("_")
            sample_data[current_file[2]][current_file[1]] = open_tsv(os.path.join(test_dir_path,filename))

    return global_data, sample_data, abundance_table

def performance_test(dir, meta):
    start_time = time.perf_counter()

    # pr = cProfile.Profile()
    # pr.enable()

    build_df(dir, meta)

    end_time = time.perf_counter()

    Total_time = end_time - start_time
    print("TEST TIME: ", Total_time)
    # pr.disable()
    # pr.print_stats(sort='time')

    # cProfile.runctx('build_df(dir,meta)', globals(), locals(), sort=1)

def produce_test_data(global_data, sample_data, abundance_data):
    global_data["main_dataframe"].to_csv("/home/lbrindel/output/postaviz/data_test/main_table.tsv",sep='\t', index=False)
    global_data["metadata"].to_csv("/home/lbrindel/output/postaviz/data_test/metadata_table.tsv",sep='\t', index=False)
    abundance_data.to_csv("/home/lbrindel/output/postaviz/data_test/abundance_table.tsv",sep='\t', index=False)
    for sample in sample_data.keys():
        sample_data[sample]["cscope"].to_csv("/home/lbrindel/output/postaviz/data_test/sample_cscope_"+sample+".tsv", sep='\t', index=False)
        sample_data[sample]["iscope"].to_csv("/home/lbrindel/output/postaviz/data_test/sample_iscope_"+sample+".tsv", sep='\t', index=False)
    quit()

def unit_test_1():
    mock_dataframe = pd.DataFrame(
        {
            "SmplID": ["CPD1", "CPD2", "CPD3", "CPD4", "CPD5", "CPD6", "CPD7"],
            "A": [1,0,0,0,1,1,0],
            "B": [0,0,1,0,1,1,0],
            "C": [1,0,0,1,1,1,1],
            "D": [0,1,0,0,0,0,0]
        }
    )
    print(mock_dataframe)
    return