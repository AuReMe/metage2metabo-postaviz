import time

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

# import plotly.graph_objects as gogit stat
import polars as pl
import seaborn as sns
from matplotlib.patches import Patch
from pandas.api.types import is_numeric_dtype

# from plotly.subplots import make_subplots
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
from scipy import stats
from statsmodels.stats.multitest import multipletests


from m2m_postaviz import data_utils as du
from m2m_postaviz.data_struct import DataStorage


def bin_exploration_processing(data: DataStorage, factor, factor_choice, rank, rank_choice, with_abundance, color, group_by_metadata = False, save_raw_data = False):
    """Takes inputs from shiny application to return 3 ploty objects:
    - hist plot of the unique production of metabolites by selected bins, weighted by abundance or not.
    - box plot of production of metabolites by bin selected.
    - bar plot of the abundance of each bin by samples.

    Each plot can be customised by the metadata from the input selected by user.

    A pre-processing is needed first to get only the bins of interest from the chunks of bins_dataframe from hard drive.

    Args:
        data (DataStorage): Data object giving access to the dataframe in disk.
        factor (str): Column of the metadata selected for filtering.
        factor_choice (str): One or several unique value from the column factor selected.
        rank (str): The taxonomic rank selected.
        rank_choice (str): The unique value of the taxonomic rank selected.
        with_abundance (bool): If the production value of the bins should be weighted by their abundance in their sample.
        color (str): Column of the metadata selected to group result by color.

    Returns:
        tuple: (Tuple(bin_production_plot_cscope, bin_production_plot_iscope), Abundance_plot, time)
    """
    start_timer = time.time()

    if not data.HAS_ABUNDANCE_DATA:
        with_abundance = False

    list_of_bins = data.get_bins_list()

    if rank == "mgs":
        rank_choice = rank_choice.split(" ")[0]

    if rank == "all":
        list_of_bin_in_rank = list_of_bins

    else:
        list_of_bin_in_rank = data.get_bin_list_from_taxonomic_rank(rank, rank_choice)

    # Taxonomic dataframe can contain MORE information and MORE bin than the data who can be a subset of the whole data. Filtering is needed.

    if rank == "all":

        filtered_list_of_bin = list_of_bin_in_rank

    else:

        filtered_list_of_bin = []

        for x in list_of_bin_in_rank:
            if x in list_of_bins:
                filtered_list_of_bin.append(x)

    if len(filtered_list_of_bin) == 0:
        print("The lenght of the list of bin in selected input is zero. Possibly because the select input list come from the taxonomic dataframe while the sample in bin_dataframe does not contain those bins.")
        return

    filter_condition=[("binID", "in", filtered_list_of_bin)]
    # if factor != "None" and len(factor_choice) > 0:
    #     filter_condition.append((factor, "in", factor_choice)) ### Metadata filter applied directly on parquet dataframe / DISABLED because the bin_dataframe no longer hold metadata. MAY CHANGE.

    figures1 = []
    raw_save = []
    for mode in ["cscope", "iscope"]:

        df = data.get_bin_dataframe(condition=filter_condition, scope_mode=mode)
        # METADATA filter applied here instead.

        metadata = data.get_metadata().to_pandas()

        df = df.merge(metadata, "inner", "smplID")

        if factor_choice != "None" and len(factor_choice) > 0:

            if is_numeric_dtype(df[factor]):

                factor_choice = list(map(float, factor_choice))

            df = df.loc[df[factor].isin(factor_choice)]

        unique_sample_in_df = df["smplID"].unique()

        new_serie_production = pd.DataFrame(columns=["smplID", "unique_production_count"])

        for sample in unique_sample_in_df:

            tmp_df = df.loc[df["smplID"] == sample][["binID","smplID","Production"]]
            all_production = tmp_df["Production"].values

            tmp_production = []

            for prod_list in all_production:

                tmp_production += list(prod_list)

            unique_production_count = len(set(tmp_production))
            new_serie_production.loc[len(new_serie_production)] = {"smplID": sample, "unique_production_count": unique_production_count}

        df = df.merge(new_serie_production, how="inner", on="smplID")

        df.sort_index(inplace=True)
        raw_save.append(df)

        if group_by_metadata and factor != "None":

            figures1.append(px.box(df, x=factor, y="Count_with_abundance" if with_abundance else "unique_production_count",
                                        color=None if color == "None" else color,
                                        hover_data="binID",
                                        labels={
                                                "unique_production_count": "unique metabolites produced",
                                                "Count_with_abundance": "abundance-weighted production of metabolites"
                                                }))
        
        else:

            figures1.append(px.histogram(df, x="smplID", y="Count_with_abundance" if with_abundance else "unique_production_count",
                                        color="smplID" if color =="None" else color,
                                        hover_data="binID",
                                        text_auto="Count_with_abundance" if with_abundance else "unique_production_count",
                                        labels={
                                                "unique_production_count": "unique metabolites produced",
                                                "Count_with_abundance": "abundance-weighted production of metabolites"
                                                }))

    # If only one bin selected do not make boxplot.

    # if len(filtered_list_of_bin) > 1:
    #     fig3 = px.box(df, x="smplID", y="Count_with_abundance" if with_abundance else "Count", color="smplID" if color =="None" else color, hover_data="binID")

    # else:
    #     fig3 = None
    if with_abundance:
        fig2 = px.bar(df, x="smplID", y="Abundance", color="Abundance", hover_data="binID")
    else:
        fig2 = None


    if save_raw_data:
        data.keep_working_dataframe("bin_production_cscope", raw_save[0])
        data.keep_working_dataframe("bin_production_iscope", raw_save[1])
        data.keep_working_dataframe("bin_abundance", df)

    return figures1, fig2, time.time() - start_timer


def global_production_statistical_dataframe(data: DataStorage, user_input1, user_input2, multiple_test_correction, correction_method, with_abundance):

    x1, x2 = user_input1, user_input2

    # No input selected
    if x1 == "None":
        return

    if multiple_test_correction:
        multipletests_method = correction_method
    else:
        multipletests_method = "hs"

    if with_abundance and data.HAS_ABUNDANCE_DATA:
        column_value = "Total_abundance_weighted"
    else:
        column_value = "Total_production"

    df = data.get_global_production_dataframe().to_pandas()

    # At least first axis selected
    if x2 == "None":
        df = df[[column_value,x1]]
        df = df.dropna()

        if is_numeric_dtype(df[x1]):

            res = du.correlation_test(df[column_value].to_numpy(), df[x1].to_numpy(), x1)

            return res

        res = du.preprocessing_for_statistical_tests(df, [column_value], x1, multipletests=multiple_test_correction, multipletests_method=multipletests_method)
        # all_dataframe["global_production_test_dataframe"] = res

        return res

    # Both axis have been selected and are not equal.
    if x1 != x2:

        df = df[[column_value,x1,x2]]
        df = df.dropna()

        if is_numeric_dtype(df[x1]):

            if is_numeric_dtype(df[x2]):

                # Double cor
                res1 = du.correlation_test(df[column_value].to_numpy(), df[x1].to_numpy(), x1)
                res2 = du.correlation_test(df[column_value].to_numpy(), df[x2].to_numpy(), x2)
                return pd.concat([res1,res2])

            else:

                # cor filtered by second categorical factor .loc
                all_results = []
                for unique_x2_value in df[x2].unique():

                    value_array = df.loc[df[x2] == unique_x2_value][column_value]
                    factor_array = df.loc[df[x2] == unique_x2_value][x1]

                    all_results.append(du.correlation_test(value_array, factor_array, unique_x2_value))

                return pd.concat(all_results)

        res = du.preprocessing_for_statistical_tests(df, [column_value], x1, x2, multipletests=multiple_test_correction, multipletests_method=multipletests_method)
        # all_dataframe["global_production_test_dataframe"] = res

        return res

    return


def metabolites_production_statistical_dataframe(data: DataStorage, metabolites_choices, user_input1, user_input2, multiple_test_correction, correction_method, save_raw_data):
    y1, x1, x2 = metabolites_choices, user_input1, user_input2

    if len(y1) == 0:
        return

    if x1 == "None":
        return

    if multiple_test_correction:
        correction_method = correction_method
    else:
        correction_method = "hs"

    if x2 == "None":
        df = data.get_cscope_producers_dataframe().to_pandas()
        df = df[[*y1,x1]]
        df = df.dropna()

        if is_numeric_dtype(df[x1]):

            correlation_results = []

            for y_value in y1:
                correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), x1))

            if save_raw_data:
                data.keep_working_dataframe("cpdtab_stats_raw_dataframe", df)

            return pd.concat(correlation_results)

        res = du.preprocessing_for_statistical_tests(df, y1, x1, multipletests = multiple_test_correction, multipletests_method= correction_method)

        if save_raw_data:
            data.keep_working_dataframe("cpdtab_stats_raw_dataframe", df)

        return res

    if x1 != x2:

        df = data.get_cscope_producers_dataframe().to_pandas()
        df = df[[*y1,x1,x2]]
        df = df.dropna()

        if is_numeric_dtype(df[x1]): # First input is Float type

            if is_numeric_dtype(df[x2]): # Second input is Float type

                correlation_results = []

                for y_value in y1:

                    correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), str(x1+" "+y_value)))
                    correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), str(x2+" "+y_value)))

                if save_raw_data:
                    data.keep_working_dataframe("cpdtab_stats_raw_dataframe", df)

                return pd.concat(correlation_results)

            else: # Second input is not Float type

                correlation_results = []

                for y_value in y1:

                    for x2_unique_value in df[x2].unique():

                        factor_array = df.loc[df[x2] == x2_unique_value][x1]
                        value_array = df.loc[df[x2] == x2_unique_value][y_value]

                        correlation_results.append(du.correlation_test(value_array.to_numpy(), factor_array.to_numpy(), str(x2_unique_value)+" "+y_value))

                if save_raw_data:
                    data.keep_working_dataframe("cpdtab_stats_raw_dataframe", df)

                return pd.concat(correlation_results)

        else:

            res = du.preprocessing_for_statistical_tests(df, y1, x1, x2, multipletests = multiple_test_correction, multipletests_method= correction_method)
            # all_dataframe["metabolites_production_test_dataframe"] = res

        if save_raw_data:
            data.keep_working_dataframe("cpdtab_stats_raw_dataframe", df)

        return res


def make_pcoa(data: DataStorage, column, choices, abundance, color):
    """Produce a Principal Coordinate Analysis with data. The Pcoa can be customized
    by filtering on specific column, using the abundance data and color the resulting plot.

    Args:
        data (DataStorage): DataStorage Object.
        column (str): Column label used for filtering.
        choices (list): Choice of the unique of the column input to use.
        abundance (bool): Option to use the column with abundance values instead of the {0 not produced ,1 produced} values.
        color (str): Column label used for the color option of the plot.

    Returns:
        px.scatter: Plotly scatter figure.
    """

    print(f"Abundance: {abundance}")

    if abundance:
        df = data.get_normalised_abundance_dataframe().to_pandas()
    else:
        df = data.get_main_dataframe().to_pandas()

    metadata = data.get_metadata().to_pandas()

    if du.is_indexed_by_id(df):
        df.reset_index(inplace=True)

    if du.is_indexed_by_id(metadata):
        metadata.reset_index(inplace=True)

    if is_numeric_dtype(metadata[column]):

        selected_sample = metadata.loc[(metadata[column] >= choices[0]) & (metadata[column] <= choices[1])]["smplID"]
        df = df.loc[df["smplID"].isin(selected_sample)]
        metadata = metadata.loc[metadata["smplID"].isin(selected_sample)]

    else:

        selected_sample = metadata.loc[metadata[column].isin(choices)]["smplID"]
        df = df.loc[df["smplID"].isin(selected_sample)]
        metadata = metadata.loc[metadata["smplID"].isin(selected_sample)]

    plot_df = run_pcoa(df, metadata)

    fig = px.scatter(plot_df, x="PC1", y="PC2", title="PCoA of reached metabolites in samples (Jaccard distance)",
                        color= color
                        )

    return fig


def run_pcoa(main_dataframe: pd.DataFrame, metadata: pd.DataFrame, distance_method: str = "jaccard"):
    """Calculate Principal Coordinate Analysis with the dataframe given in args.
    Use metadata's drataframe as second argument to return the full ordination result plus
    all metadata column inserted along Ordination.samples dataframe.
    Ready to be plotted.

    Args:
        main_df (pd.DataFrame): Main dataframe of compound production
        metadata (pd.DataFrame): Metadata's dataframe

    Returns:
        pd.DataFrame: Ordination results object from skbio's package.
    """

    if not du.is_indexed_by_id(main_dataframe):
        main_dataframe = main_dataframe.set_index("smplID")

    if du.is_indexed_by_id(metadata):
        metadata = metadata.reset_index("smplID")

    dmatrix = main_dataframe.to_numpy()
    dist_m = pdist(dmatrix, "jaccard")
    square_df = squareform(dist_m)
    pcoa_result = pcoa(square_df,number_of_dimensions=2)
    coordinate = pcoa_result.samples

    df_pcoa = coordinate[["PC1","PC2"]]
    df_pcoa["smplID"] = main_dataframe.index.to_numpy()

    df_pcoa = df_pcoa.merge(metadata, "inner", "smplID")
    df_pcoa.set_index("smplID",inplace=True)

    return df_pcoa


def render_reactive_total_production_plot(data: DataStorage, user_input1, user_input2, with_abundance):
    """Produce and return a plotly figure object. Barplot or Boxplot if there is only unique value in columns.

    Args:
        data (DataStorage): DataStorage object.
        user_input1 (_type_): Column label for metadata filtering.
        user_input2 (_type_): Column label for metadata filtering.
        with_abundance (bool): Option to use the column with abundance values instead of the {0 not produced ,1 produced} values.

    Returns:
        px.box: Plotly express object.
        pd.Dataframe: dataframe used for the plot.
    """
    if with_abundance and data.HAS_ABUNDANCE_DATA:
        column_value = "Total_abundance_weighted"
    else:
        column_value = "Total_production"

    df = data.get_global_production_dataframe().to_pandas()

    if user_input1 == "None":
        fig = px.box(df, y=column_value, title="Numbers of unique metabolites reached by sample.")
        return fig, df

    elif user_input2 == "None" or user_input1 == user_input2:

        if is_numeric_dtype(df[user_input1]):
            df.sort_values(user_input1, inplace=True)

        df = df[[column_value,user_input1]]
        df = df.dropna()

        if du.has_only_unique_value(df, user_input1):

            return px.bar(df, x=user_input1 , y=column_value, color=user_input1, title=f"Numbers of unique compounds produced by samples filtered by {user_input1}")

        else:

            fig = px.box(df, x=user_input1 , y=column_value, color=user_input1, title=f"Numbers of unique compounds produced by samples filtered by {user_input1}")
            fig.update_xaxes(type="category")
            return fig, df

    else:

        if is_numeric_dtype(df[user_input1]):
            df.sort_values(user_input1, inplace=True)

        df = df[[column_value,user_input1,user_input2]]
        df = df.dropna()
        has_unique_value = du.has_only_unique_value(df, user_input1, user_input2)
        if has_unique_value:
            fig = px.bar(df,x=user_input1,y=column_value,color=user_input2)
        else:
            fig = px.box(df,x=user_input1,y=column_value,color=user_input2)
            fig.update_xaxes(type="category")
        return fig, df


def render_reactive_metabolites_production_plot(data: DataStorage, compounds_input, user_input1, color_input = "None", sample_filter_button = "All", sample_filter_value = [], with_abundance = None, save_raw_data = False):

    if len(compounds_input) == 0:
        return

    if not data.HAS_ABUNDANCE_DATA:
        with_abundance = False

    if with_abundance:
        cscope_df = data.get_normalised_abundance_dataframe(with_metadata=True)
        iscope_df = data.get_normalised_iscope_abundance_dataframe(with_metadata=True)
    else:
        cscope_df = data.get_cscope_producers_dataframe()
        iscope_df = data.get_iscope_producers_dataframe()
    # cscope_df_iscope = data.get_iscope_metabolite_production_dataframe()

    if sample_filter_button != "All" and len(sample_filter_value) != 0:

        if sample_filter_button == "Include":

            cscope_df = cscope_df.filter(pl.col("smplID").is_in(sample_filter_value))
            iscope_df = cscope_df.filter(pl.col("smplID").is_in(sample_filter_value))

        elif sample_filter_button == "Exclude":

            cscope_df = cscope_df.filter(~pl.col("smplID").is_in(sample_filter_value))
            iscope_df = cscope_df.filter(~pl.col("smplID").is_in(sample_filter_value))

    if user_input1 == "None":

        if color_input == "None":

            df = cscope_df.select([*compounds_input])
            dfi = iscope_df.select([*compounds_input])
            fig1 = px.box(df, y=compounds_input).update_layout(yaxis_title="Numbers of metabolic network producers for each sample")
            fig2 = px.box(dfi, y=compounds_input).update_layout(yaxis_title="Numbers of metabolic network producers for each sample")

            if save_raw_data:
                data.keep_working_dataframe("cpd_producers_cscope", df)
                data.keep_working_dataframe("cpd_producers_iscope", dfi)
            return fig1, fig2

        else:

            df = cscope_df.select([*compounds_input, color_input])
            dfi = iscope_df.select([*compounds_input, color_input])
            has_unique_value = du.has_only_unique_value(df, color_input)

            if has_unique_value:
                fig = px.bar(df, y=compounds_input, color=color_input).update_layout(yaxis_title="Numbers of metabolic network producers")
                fig2 = px.bar(dfi, y=compounds_input, color=color_input).update_layout(yaxis_title="Numbers of metabolic network producers")
            else:
                fig = px.box(df, y=compounds_input, color=color_input).update_layout(yaxis_title="Numbers of metabolic network producers")
                fig2 = px.box(dfi, y=compounds_input, color=color_input).update_layout(yaxis_title="Numbers of metabolic network producers")

            if save_raw_data:
                data.keep_working_dataframe("cpd_producers_cscope", df)
                data.keep_working_dataframe("cpd_producers_iscope", dfi)
            return fig, fig2

    if color_input == "None" or user_input1 == color_input: # If only the filtering by metadata has been selected (no color).

        df = cscope_df.select([*compounds_input,user_input1])
        df = df.drop_nulls()
        dfi = iscope_df.select([*compounds_input,user_input1])
        dfi = dfi.drop_nulls()

        if is_numeric_dtype(df[user_input1]):
            df.sort_values(user_input1, inplace=True)
            dfi.sort_values(user_input1, inplace=True)

        has_unique_value = du.has_only_unique_value(df, user_input1)

        if df.get_column(user_input1).dtype.is_numeric():

            df = df.sort(user_input1)
            df = df.with_columns(pl.col(user_input1).cast(pl.String))

            dfi = dfi.sort(user_input1)
            dfi = dfi.with_columns(pl.col(user_input1).cast(pl.String))

        if has_unique_value:
            fig = px.bar(df, x=user_input1, y=compounds_input, color=user_input1).update_layout(yaxis_title="Numbers of metabolic network producers")
            fig2 = px.bar(dfi, x=user_input1, y=compounds_input, color=user_input1).update_layout(yaxis_title="Numbers of metabolic network producers")
        else:
            fig = px.box(df, x=user_input1, y=compounds_input, color=user_input1).update_layout(yaxis_title="Numbers of metabolic network producers")
            fig.update_xaxes(type="category")
            fig2 = px.box(dfi, x=user_input1, y=compounds_input, color=user_input1).update_layout(yaxis_title="Numbers of metabolic network producers")
            fig2.update_xaxes(type="category")

        if save_raw_data:
            data.keep_working_dataframe("cpd_producers_cscope", df)
            data.keep_working_dataframe("cpd_producers_iscope", dfi)
        return fig, fig2

    df = cscope_df.select([*compounds_input,user_input1,color_input])
    df = df.drop_nulls()
    dfi = iscope_df.select([*compounds_input,user_input1,color_input])
    dfi = dfi.drop_nulls()

    if is_numeric_dtype(df[user_input1]):
        df.sort_values(user_input1, inplace=True)
        dfi.sort_values(user_input1, inplace=True)

    has_unique_value = du.has_only_unique_value(df, user_input1, color_input)

    if df.get_column(user_input1).dtype.is_numeric():
        df = df.sort(user_input1)
        df = df.with_columns(pl.col(user_input1).cast(pl.String))

        dfi = dfi.sort(user_input1)
        dfi = dfi.with_columns(pl.col(user_input1).cast(pl.String))

    if has_unique_value:
        fig = px.bar(df, x=user_input1, y=compounds_input,color=color_input).update_layout(yaxis_title="Numbers of metabolic network producers")
        fig.update_xaxes(type="category")
        fig2 = px.bar(dfi, x=user_input1, y=compounds_input,color=color_input).update_layout(yaxis_title="Numbers of metabolic network producers")
        fig2.update_xaxes(type="category")
    else:
        fig = px.box(df, x=user_input1, y=compounds_input, color=color_input, boxmode="group").update_layout(yaxis_title="Numbers of metabolic network producers")
        fig.update_xaxes(type="category")
        fig2 = px.box(dfi, x=user_input1, y=compounds_input, color=color_input, boxmode="group").update_layout(yaxis_title="Numbers of metabolic network producers")
        fig2.update_xaxes(type="category")
        if save_raw_data:
            data.keep_working_dataframe("cpd_producers_cscope", df)
            data.keep_working_dataframe("cpd_producers_iscope", dfi)
    return fig, fig2


def df_to_plotly(df):
    return {"z": df.values.tolist(),
            "x": df.columns.tolist(),
            "y": df.index.tolist()}


def percentage_smpl_producing_cpd(data: DataStorage, cpd_input: list, metadata_filter_input: str, sample_filter_button = "All", sample_filter_value = [], save_raw_data = False):
    """Produce two plotly figure barplot from the list of compounds and the column filter given in input.

    Args:
        data (DataStorage): DataStorage object
        cpd_input (list): List of compounds input
        metadata_filter_input (str): Column label of metadata filter

    Returns:
        Tuple: Tuple with cscope plot and iscope plot
    """
    cscope_df = data.get_cscope_producers_dataframe()
    iscope_df = data.get_iscope_producers_dataframe()

    cscope_df = cscope_df.to_pandas()
    iscope_df = iscope_df.to_pandas()

    # Check if the iscope dataframe contain all the cpd in cscope dataframe. IF not add them as column filled with 0 value.
    col_diff = cscope_df.columns.difference(iscope_df.columns)

    col_diff_dict = dict.fromkeys(col_diff, 0.0)

    temp_df = pd.DataFrame(col_diff_dict, index=iscope_df.index)

    iscope_df = pd.concat([iscope_df, temp_df], axis=1)

    # Select only the column of interest.
    cscope_df = cscope_df[["smplID", *cpd_input, metadata_filter_input]].dropna()
    iscope_df = iscope_df[["smplID", *cpd_input, metadata_filter_input]].dropna()

    # Samples filtering
    if sample_filter_button != "All":

        if sample_filter_button == "Include":

            cscope_df = cscope_df.loc[cscope_df["smplID"].isin(sample_filter_value)]
            iscope_df = iscope_df.loc[iscope_df["smplID"].isin(sample_filter_value)]

        if sample_filter_button == "Exclude":

            cscope_df = cscope_df.loc[~cscope_df["smplID"].isin(sample_filter_value)]
            iscope_df = iscope_df.loc[~iscope_df["smplID"].isin(sample_filter_value)]

    # Check for numeric dtype (boolean / int / unsigned / float / complex).
    if cscope_df[metadata_filter_input].dtype.kind in "biufc":
        cscope_df[metadata_filter_input] = cscope_df[metadata_filter_input].astype("str")

    if iscope_df[metadata_filter_input].dtype.kind in "biufc":
        iscope_df[metadata_filter_input] = iscope_df[metadata_filter_input].astype("str")

    # Set Id and metadata column as index to get the matrix.
    cscope_df.set_index(["smplID", metadata_filter_input], inplace=True)
    iscope_df.set_index(["smplID", metadata_filter_input], inplace=True)

    # Replace any value above 0 by 1.
    cscope_df.mask(cscope_df > 0.0, 1, inplace=True)
    iscope_df.mask(iscope_df > 0.0, 1, inplace=True)

    cscope_df.reset_index(inplace=True)
    iscope_df.reset_index(inplace=True)

    cscope_series = []
    # Loop throught sub dataframe of each unique value of metadata input.
    for metadata_value in cscope_df[metadata_filter_input].unique():

        current_rows = cscope_df.loc[cscope_df[metadata_filter_input] == metadata_value]

        current_rows.set_index(["smplID", metadata_filter_input], inplace=True)

        current_rows = current_rows.apply(col_value_to_percent, axis=0)

        current_rows.name = metadata_value

        cscope_series.append(current_rows)

    iscope_series = []

    for metadata_value in cscope_df[metadata_filter_input].unique():

        current_rows = iscope_df.loc[iscope_df[metadata_filter_input] == metadata_value]

        current_rows.set_index(["smplID", metadata_filter_input], inplace=True)

        current_rows = current_rows.apply(col_value_to_percent, axis=0)

        current_rows.name = metadata_value

        iscope_series.append(current_rows)

    cscope_df = pd.concat(cscope_series, axis=1)
    cscope_df = cscope_df.T

    iscope_df = pd.concat(iscope_series, axis=1)
    iscope_df = iscope_df.T

    cscope_df.reset_index(inplace=True)
    cscope_df.rename(columns={"index" : "metadata"}, inplace=True)
    cscope_df = cscope_df.melt("metadata")

    iscope_df.reset_index(inplace=True)
    iscope_df.rename(columns={"index" : "metadata"}, inplace=True)
    iscope_df = iscope_df.melt("metadata")

    #TODO legend title should be the metadata column label, not "matadata"

    fig1 = px.bar(cscope_df, x = "variable", y = "value", color = "metadata", barmode="group", title="Percentage of sample producing the selected metabolites (interacting community)", labels={"metadata": metadata_filter_input}).update_layout(bargap=0.2,xaxis_title="Compounds", yaxis_title="Percent of sample producing the compound")

    fig2 = px.bar(iscope_df, x = "variable", y = "value", color = "metadata", barmode="group", title="Percentage of sample producing the selected metabolites (non-interacting community)", labels={"metadata": metadata_filter_input}).update_layout(bargap=0.2,xaxis_title="Compounds", yaxis_title="Percent of sample producing the compound")

    if save_raw_data:
        data.keep_working_dataframe("percentage_producers_cscope", cscope_df)
        data.keep_working_dataframe("percentage_producers_iscope", iscope_df)

    return fig1, fig2


def col_value_to_percent(col: pd.Series):

    sum_val = col.sum()

    len_val = len(col.values)

    final_val = (sum_val / len_val) * 100

    return final_val


def sns_clustermap(data: DataStorage, cpd_input, metadata_input = None, row_cluster = False, col_cluster = False, filter_mode = None, filter_values = None, save_raw_data = False):
    """Produce a customizable Seaborn clustermap. Distance matrix use the jaccard method when clustering enabled.

    Args:
        data (DataStorage): DataStorage object.
        cpd_input (list): list of compounds input to filter.
        metadata_input (str, optional): Column label to filter sample by their metadata. Add a ROW color. Defaults to None.
        row_cluster (bool, optional): Dendogram for rows from distance matrix. Defaults to False.
        col_cluster (bool, optional): Dendogram for cols from distance matrix. Defaults to False.
        filter_mode (str, optional): Mode of sample's filter if enabled. Defaults to None.
        filter_values (list, optional): list of samples to filter. Defaults to None.

    Returns:
        list: List of three clustermap matrix object.
    """
    matplotlib.use("pdf")
    plots = []
    raw_dataframes = []

    for dataframe in data.get_added_value_dataframe(cpd_input, filter_mode, filter_values):

        dataframe = dataframe.to_pandas()
        dataframe.set_index("smplID", inplace = True)

        if metadata_input is not None and metadata_input != "None":

            metadata = data.get_metadata()
            metadata = metadata.to_pandas()
            metadata = metadata[["smplID",metadata_input]]
            metadata.set_index("smplID",inplace=True)

            dataframe = dataframe.merge(metadata, right_index=True, left_index=True)
            dataframe.dropna(inplace=True)
            dataframe.sort_values(metadata_input,inplace=True)

            lut_size = sns.husl_palette(len(dataframe[metadata_input].unique().tolist()), s=1, l=0.5)
            lut = dict(zip(dataframe[metadata_input].unique(), lut_size))
            row_colors = dataframe[metadata_input].map(lut)

            dataframe.drop(columns=metadata_input,inplace=True)

            g = sns.clustermap(dataframe, row_colors=row_colors, metric="jaccard", col_cluster=col_cluster, row_cluster=row_cluster, cbar_pos=(.9, .2, .03, .2), xticklabels=True)
            handles = [Patch(facecolor=lut[name]) for name in lut]
            plt.legend(handles, lut, title=metadata_input,
                       
            bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc="upper right")
            if save_raw_data:
                raw_dataframes.append(dataframe)
            plots.append(g)

        else :

            g = sns.clustermap(dataframe, metric="jaccard", col_cluster=col_cluster, row_cluster=row_cluster, cbar_pos=(.9, .2, .03, .2), xticklabels=True)

            if save_raw_data:
                raw_dataframes.append(dataframe)
            plots.append(g)

    if save_raw_data:
        for i in range(len(raw_dataframes)):
            name = "clustermap_seaborn_raw_dataframe_" + str(i)
            data.keep_working_dataframe(name, raw_dataframes[i])

    return plots


def cpd_reached_plot(data: DataStorage, metadata_input: str, multiple_correction, correction_method):
    """Produce and return a plotly.express boxplot of the compounds reached by the sample in community, individually or not reached.
    The plot can be grouped by the metadata.
    Args:
        data (DataStorage): DataStorage object.
        metadata_input (str): Metadata column label.

    Returns:
        plotly.express.boxplot: Plotly boxplot
    """
    # total_cpd = len(data.get_compound_list())

    # Comm_reach
    df = data.get_cscope_producers_dataframe(with_metadata=False)
    df = df.with_columns(pl.when(pl.exclude("smplID") > 0).then(1).otherwise(0).name.keep())
    df = df.unpivot(index=["smplID"],value_name="community-reached")
    df = df.with_columns(pl.exclude("smplID", "variable").sum().over("smplID").name.keep())
    df = df.select("smplID", "community-reached")
    df = df.group_by("smplID").agg(pl.mean("community-reached").cast(pl.Int32))

    # df = df.with_columns(pl.lit(total_cpd).alias("Unreached"))
    # df = df.with_columns(pl.col("Unreached").sub(pl.col("Comm_reached")))

    # Ind_reach
    dfi = data.get_iscope_producers_dataframe(with_metadata=False)
    dfi = dfi.with_columns(pl.when(pl.exclude("smplID") > 0).then(1).otherwise(0).name.keep())
    dfi = dfi.unpivot(index=["smplID"],value_name="individually-reached")
    dfi = dfi.with_columns(pl.exclude("smplID", "variable").sum().over("smplID").name.keep())
    dfi = dfi.select("smplID", "individually-reached")
    dfi = dfi.group_by("smplID").agg(pl.mean("individually-reached").cast(pl.Int32))

    df = df.join(dfi, on="smplID")

    if metadata_input == "None" or metadata_input == "smplID":

        df = df.unpivot(index=["smplID"])
        # return px.bar(df, x="smplID",y="value",color="variable",barmode="group")
        return

    else:

        metadata = data.get_metadata().select("smplID",metadata_input)
        df = df.join(metadata, on ="smplID")
        df = df.unpivot(index=["smplID",metadata_input])

        if df.get_column(metadata_input).dtype.is_numeric():
            df = df.sort(metadata_input)

        stat_df = reached_compounds_plot_stats_tests(df, metadata_input, multiple_correction, correction_method)
        # data.set_overview_reachplot_stats_dataframe(stat_df)

        fig = px.box(df, x="variable",y="value", title="Individually- and community-reached metabolites in samples", color=metadata_input)
        fig.update_xaxes(type="category")

        return fig, stat_df


def wilcoxon_mann_whitney(data, metadata_input, context, multiple_correction, multiple_correction_method):
    """Receive a dictionnary whose key/value be used to apply wilcoxon test if they value array are the same lenght. Mann-Whitney otherwise.

    Args:
        data (Dict): Dictionnary of the Metadata {unique metadata value : [ value array ]}
        metadata_input (_type_): _description_

    Returns:
        pd.Dataframe: Pandas dataframe of the resulting tests.
    """
    results = pd.DataFrame(columns=["Context", "Metadata_column", "Factor1", "Sample size1", "Factor2", "Sample size2", "Method", "Statistic", "Pvalue", "Significance"])

    for factor, value in data.items():

        for factor2, value2 in data.items():

            if factor == factor2:
                continue

            if len(value) < 1 or len(value2) < 1:
                continue

            if factor2 in results["Factor1"].tolist():
                continue

            if len(value) == len(value2):
                test_value, pvalue = stats.wilcoxon(value, value2)
                test_method = "Wilcoxon"

            else:
                test_value, pvalue = stats.mannwhitneyu(value, value2)
                test_method = "Mann-Whitney"

            if pvalue >= 0.05:
                symbol = "ns"
            elif pvalue >= 0.01:
                symbol = "*"
            elif pvalue >= 0.001:
                symbol = "**"
            else:
                symbol = "***"

            results.loc[len(results)] = {"Context": context,"Metadata_column": metadata_input,
                                        "Factor1": factor, "Sample size1": len(value),
                                        "Factor2": factor2, "Sample size2": len(value2),
                                        "Statistic": test_value, "Pvalue": pvalue, "Significance": symbol, "Method": test_method
                                        }

    if multiple_correction:

        pvals_before_correction = results["Pvalue"].to_numpy()
        reject, pvals_after_correction, _, __ = multipletests(pvals_before_correction, method = multiple_correction_method)

        results["Pvalue corrected"] = pvals_after_correction
        results["Significance corrected"] = results["Pvalue corrected"].apply(lambda x:get_significance_symbol(x))
        results["Correction method"] = multiple_correction_method

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


def split_value_column_with_metadata(df, metadata_column):
    """Split the metadata column of the dataframe into dictionnary whose keys are the unique value of the column
    and value as a list of value from "value" column corresponding to th metadata value.

    Example:

    ┌────────────┬──────┬───────┐
    │ smplID     ┆ Days ┆ value │   
    │ ---        ┆ ---  ┆ ---   │
    │ str        ┆ i64  ┆ i64   │
    ╞════════════╪══════╪═══════╡
    │ ERAS1d0    ┆ 0    ┆ 1027  │
    │ ERAS2d0    ┆ 0    ┆ 1021  │
    │ ERAS3d0    ┆ 0    ┆ 942   │
    │ ERAS4d0    ┆ 0    ┆ 1086  │
    │ ERAS5d0    ┆ 0    ┆ 1069  │
    │ ERAS8d180  ┆ 180  ┆ 1034  │
    │ ERAS9d180  ┆ 180  ┆ 1040  │
    │ ERAS10d180 ┆ 180  ┆ 1061  │
    │ ERAS11d180 ┆ 180  ┆ 1105  │
    │ ERAS12d180 ┆ 180  ┆ 1027  │
    └────────────┴──────┴───────┘

    Expected result: {0: [1027, 1021, 942, 1086, 1069, 1034], 180: [1034, 1034, 1040, 1061, 1105, 1027]}

    Args:
        df (pl.Dataframe): Polars dataframe.
    """

    data_to_test = {}

    for unique_value in df[metadata_column].unique().__iter__():

        data_to_test[unique_value] = df.filter(pl.col(metadata_column) == unique_value).select(pl.col("value")).to_series().to_list()

    return data_to_test


def reached_compounds_plot_stats_tests(df, metadata_input, multiple_correction, correction_method):
    """Takes the reached_cpd_plot dataframe to process the different combination for the Wilcoxon/Whitney test.

    Args:
        df (pl.Dataframe): Dataframe of the plot.
        metadata_input (_type_): Metadata column choosed in input.

    Raises:
        TypeError: The dataframe in input must be a Polars dataframe.

    Returns:
        Dataframe: Dataframe of the results
    """
    if isinstance(df, pl.DataFrame):

        df_comm = df.filter(pl.col("variable") == "community-reached").drop("variable")

        df_ind = df.filter(pl.col("variable") == "individually-reached").drop("variable")

    else:

        raise TypeError(f"Dataframe received for the reach_compounds_plot stat function should be a Polars dataframe. not {type(df)}")

    comm_data = split_value_column_with_metadata(df_comm, metadata_input)
    ind_data = split_value_column_with_metadata(df_ind, metadata_input)

    comm_res = wilcoxon_mann_whitney(comm_data, metadata_input, "community", multiple_correction, correction_method)
    ind_res = wilcoxon_mann_whitney(ind_data, metadata_input, "individually", multiple_correction, correction_method)

    res = pd.concat([comm_res, ind_res])

    return res