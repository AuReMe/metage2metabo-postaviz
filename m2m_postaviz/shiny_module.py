import time

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa

from m2m_postaviz import data_utils as du
from m2m_postaviz.data_struct import DataStorage


def bin_exploration_processing(data: DataStorage, factor, factor_choice, rank, rank_choice, with_abundance, color):
    """Takes inputs from shiny application to return 3 ploty objects:
    - hist plot of the unique production of metabolites by bins selected weighted by abundance or not.
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
        tuple: (Histplot, boxplot, the dataframe created from the filters, time, barplot)
    """
    start_timer = time.time()

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

    filter_condition=[("binID", "in", filtered_list_of_bin)]
    # if factor != "None" and len(factor_choice) > 0:
    #     filter_condition.append((factor, "in", factor_choice))

    df = data.get_bin_dataframe(condition=filter_condition)

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

    metadata = data.get_metadata()
    metadata = metadata.loc[metadata["smplID"].isin(unique_sample_in_df)]

    df = df.merge(metadata, "inner", "smplID")

    if factor_choice == "None":

        df = df.loc[df[factor].isin(factor_choice)]

    df.sort_index(inplace=True)

    fig1 = px.histogram(df, x="smplID", y="Count_with_abundance" if with_abundance else "unique_production_count", color="smplID" if color =="None" else color, hover_data="binID")

    # If only one bin selected do not make boxplot.

    if len(filtered_list_of_bin) > 1:
        fig3 = px.box(df, x="smplID", y="Count_with_abundance" if with_abundance else "Count", color="smplID" if color =="None" else color, hover_data="binID")

    else:
        fig3 = None

    fig2 = px.bar(df, x="smplID", y="Abundance", color="Abundance", hover_data="binID")

    return fig1, fig2, df, time.time() - start_timer, fig3


def global_production_statistical_dataframe(data: DataStorage, user_input1, user_input2, multiple_test_correction, correction_method, with_abundance):

            x1, x2 = user_input1, user_input2

            # No input selected
            if x1 == "None":
                return

            multipletests_correction = multiple_test_correction

            if multipletests_correction:
                multipletests_method = correction_method
            else:
                multipletests_method = "hs"

            with_normalised_data = with_abundance

            if with_normalised_data and data.HAS_ABUNDANCE_DATA:
                column_value = "Total_abundance_weighted"
            else:
                column_value = "Total_production"

            df = data.get_global_production_dataframe()

            # At least first axis selected
            if x2 == "None":
                df = df[[column_value,x1]]
                df = df.dropna()

                if du.serie_is_float(df[x1]):

                    res = du.correlation_test(df[column_value].to_numpy(), df[x1].to_numpy(), x1)

                    return res

                res = du.preprocessing_for_statistical_tests(df, [column_value], x1, multipletests=multipletests_correction, multipletests_method=multipletests_method)
                # all_dataframe["global_production_test_dataframe"] = res

                return res

            # Both axis have been selected and are not equal.
            if x1 != x2:

                df = df[[column_value,x1,x2]]
                df = df.dropna()

                if du.serie_is_float(df[x1]):

                    if du.serie_is_float(df[x2]):

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

                res = du.preprocessing_for_statistical_tests(df, [column_value], x1, x2, multipletests=multipletests_correction, multipletests_method=multipletests_method)
                # all_dataframe["global_production_test_dataframe"] = res

                return res

            return


def metabolites_production_statistical_dataframe(data: DataStorage, metabolites_choices, user_input1, user_input2, multiple_test_correction, correction_method, with_abundance = None):
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

        df = data.get_metabolite_production_dataframe()[[*y1,x1]]
        df = df.dropna()

        if du.serie_is_float(df[x1]):

            correlation_results = []

            for y_value in y1:

                correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), x1))

            return pd.concat(correlation_results)


        res = du.preprocessing_for_statistical_tests(df, y1, x1, multipletests = multiple_test_correction, multipletests_method= correction_method)
        # all_dataframe["metabolites_production_test_dataframe"] = res

        return res

    if x1 != x2:

        df = data.get_metabolite_production_dataframe()[[*y1,x1,x2]]
        df = df.dropna()

        if du.serie_is_float(df[x1]): # First input is Float type

            if du.serie_is_float(df[x2]): # Second input is Float type

                correlation_results = []

                for y_value in y1:

                    correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), str(x1+" "+y_value)))
                    correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), str(x2+" "+y_value)))

                return pd.concat(correlation_results)

            else: # Second input is not Float type

                correlation_results = []

                for y_value in y1:

                    for x2_unique_value in df[x2].unique():

                        factor_array = df.loc[df[x2] == x2_unique_value][x1]
                        value_array = df.loc[df[x2] == x2_unique_value][y_value]

                        correlation_results.append(du.correlation_test(value_array.to_numpy(), factor_array.to_numpy(), str(x2_unique_value)+" "+y_value))

                return pd.concat(correlation_results)

        else:

            res = du.preprocessing_for_statistical_tests(df, y1, x1, x2, multipletests = multiple_test_correction, multipletests_method= correction_method)
            # all_dataframe["metabolites_production_test_dataframe"] = res

        return res


def make_pcoa(data: DataStorage, column, choices, abundance, color):

    if abundance:
        df = data.get_normalised_abundance_dataframe()
    else:
        df = data.get_main_dataframe()

    metadata = data.get_metadata()

    if du.is_indexed_by_id(df):
        df.reset_index(inplace=True)

    if du.is_indexed_by_id(metadata):
        metadata.reset_index(inplace=True)

    if du.serie_is_float(metadata[column]):

        selected_sample = metadata.loc[(metadata[column] >= choices[0]) & (metadata[column] <= choices[1])]["smplID"]
        df = df.loc[df["smplID"].isin(selected_sample)]
        metadata = metadata.loc[metadata["smplID"].isin(selected_sample)]

    else:

        selected_sample = metadata.loc[metadata[column].isin(choices)]["smplID"]
        df = df.loc[df["smplID"].isin(selected_sample)]
        metadata = metadata.loc[metadata["smplID"].isin(selected_sample)]

    plot_df = run_pcoa(df, metadata)

    fig = px.scatter(plot_df, x="PC1", y="PC2",
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
        _type_: Ordination results object from skbio's package.
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

    if with_abundance and data.HAS_ABUNDANCE_DATA:
        column_value = "Total_abundance_weighted"
    else:
        column_value = "Total_production"

    df = data.get_global_production_dataframe()

    if user_input1 == "None":

        # all_dataframe["global_production_plot_dataframe"] = df
        return px.box(df, y=column_value, title="Numbers of unique compounds produced by sample.")

    elif user_input2 == "None" or user_input1 == user_input2:

        df = df[[column_value,user_input1]]
        df = df.dropna()
        # all_dataframe["global_production_plot_dataframe"] = df

        if du.has_only_unique_value(df, user_input1):

            return px.bar(df, x=user_input1 , y=column_value, color=user_input1, title=f"Numbers of unique compounds produced by samples filtered by {user_input1}")

        else:

            fig = px.box(df, x=user_input1 , y=column_value, color=user_input1, title=f"Numbers of unique compounds produced by samples filtered by {user_input1}")
            return fig

    else:

        df = df[[column_value,user_input1,user_input2]]
        df = df.dropna()
        # all_dataframe["global_production_plot_dataframe"] = df
        has_unique_value = du.has_only_unique_value(df, user_input1, user_input2)

        return px.bar(df,x=user_input1,y=column_value,color=user_input2) if has_unique_value else px.box(df,x=user_input1,y=column_value,color=user_input2)


def render_reactive_metabolites_production_plot(data: DataStorage, compounds_input, user_input1, color_input, sample_filtering_enabled = False, sample_filter_button = "None", sample_filter_value = [], with_abundance = None):

    if len(compounds_input) == 0:
        return

    producer_data = data.get_metabolite_production_dataframe()
    # producer_data_iscope = data.get_iscope_metabolite_production_dataframe()

    if sample_filtering_enabled and len(sample_filter_value) != 0:

        if sample_filter_button == "Include":

            producer_data = producer_data.loc[producer_data["smplID"].isin(sample_filter_value)]

        elif sample_filter_button == "Exclude":

            producer_data = producer_data.loc[~producer_data["smplID"].isin(sample_filter_value)]

    # producer_data = producer_data.set_index("smplID")

    if user_input1 == "None":

        if color_input == "None":

            df = producer_data[[*compounds_input]]
            return px.box(df, y=compounds_input).update_layout(yaxis_title="Numbers of mgs producers for each samples")

        else:

            df = producer_data[[*compounds_input, color_input]]
            has_unique_value = du.has_only_unique_value(df, color_input)

            if has_unique_value:

                fig = px.bar(df, y=compounds_input, color=color_input).update_layout(yaxis_title="Numbers of genomes producers")
            else:
                fig = px.box(df, y=compounds_input, color=color_input).update_layout(yaxis_title="Numbers of genomes producers")

            return fig

    if color_input == "None" or user_input1 == color_input: # If only the filtering by metadata has been selected (no color).

        df = producer_data[[*compounds_input,user_input1]]
        df = df.dropna()

        has_unique_value = du.has_only_unique_value(df, user_input1)

        if df[user_input1].dtypes == float or df[user_input1].dtypes == int:

            df[user_input1] = df[user_input1].astype(str)

        if has_unique_value:

            fig = px.bar(df, x=user_input1, y=compounds_input, color=user_input1).update_layout(yaxis_title="Numbers of genomes producers")
        else:
            fig = px.box(df, x=user_input1, y=compounds_input, color=user_input1).update_layout(yaxis_title="Numbers of genomes producers")

        return fig

    df = producer_data[[*compounds_input,user_input1,color_input]]
    df = df.dropna()

    has_unique_value = du.has_only_unique_value(df, user_input1, color_input)

    if df[user_input1].dtypes == float or df[user_input1].dtypes == int:

        df[user_input1] = df[user_input1].astype(str)

    if has_unique_value:

        fig = px.bar(df, x=user_input1, y=compounds_input,color=color_input).update_layout(yaxis_title="Numbers of genomes producers")
    else:
        fig = px.box(df, x=user_input1, y=compounds_input, color=color_input, boxmode="group").update_layout(yaxis_title="Numbers of genomes producers")

    return fig


def metabolites_heatmap(data: DataStorage, user_cpd_list, user_input1, sample_filtering_enabled, sample_filter_button, sample_filter_value, by = None):

    all_cscope_dataframe = []
    all_iscope_dataframe = []

    for cpd in user_cpd_list:

        try:

            tmp_df = data.get_minimal_cpd_dataframe(cpd)
            all_cscope_dataframe.append(tmp_df[0])
            all_iscope_dataframe.append(tmp_df[1])

        except Exception:

            print(f"{cpd} not in icscope dataframe")

    dfc = pd.concat(all_cscope_dataframe,axis=0)
    dfi = pd.concat(all_iscope_dataframe,axis=0)

    if sample_filtering_enabled and len(sample_filter_value) != 0:

        if sample_filter_button == "Include":

            dfc = dfc.loc[dfc["smplID"].isin(sample_filter_value)]
            dfi = dfi.loc[dfi["smplID"].isin(sample_filter_value)]

        elif sample_filter_button == "Exclude":

            dfc = dfc.loc[~dfc["smplID"].isin(sample_filter_value)]
            dfi = dfi.loc[~dfi["smplID"].isin(sample_filter_value)]

    if by == "taxonomy" and user_input1 != "None":

        taxonomy_file = data.get_taxonomic_dataframe()

        mgs_col_taxonomy = taxonomy_file.columns[0]

        dfc["binID"] = dfc["binID"].map(taxonomy_file.set_index(mgs_col_taxonomy)[user_input1])

        dfi["binID"] = dfi["binID"].map(taxonomy_file.set_index(mgs_col_taxonomy)[user_input1])

    if by == "metadata" and user_input1 != "None":

        metadata = data.get_metadata()

        dfc["smplID"] = dfc["smplID"].map(metadata.set_index("smplID")[user_input1])

        dfi["smplID"] = dfi["smplID"].map(metadata.set_index("smplID")[user_input1])

    dfc.fillna(0,inplace=True)
    # print(dfc)

    dfc.set_index(["binID","smplID"],inplace=True)
    dfc = dfc.groupby(level="binID").sum()

    dfi.fillna(0,inplace=True)
    # print(dfi)

    dfi.set_index(["binID","smplID"],inplace=True)
    dfi = dfi.groupby(level="binID").sum()



    fig = make_subplots(1,2, subplot_titles=["cscope", "iscope"])

    fig.add_trace(go.Heatmap(df_to_plotly(dfc), coloraxis="coloraxis", ),1,1)

    fig.add_trace(go.Heatmap(df_to_plotly(dfi), coloraxis="coloraxis", ),1,2)

    fig.update_layout(coloraxis=dict(colorscale = "ylgnbu"))

    return fig


def df_to_plotly(df):
    return {"z": df.values.tolist(),
            "x": df.columns.tolist(),
            "y": df.index.tolist()}


def added_value_heatmap(data: DataStorage, cpd_input : list, sample_filtering_enabled, sample_filter_mode, sample_filter_value):
    """Get three dataframe from the DataStorage object. cscope producers, icscope_producers and the difference between the two.
    The dataframe have been filtered by DataStorage previously by the 4 input given in args. 
    Make three plotly Heatmap from those dataframe and return them to be rendered. 
    Args:
        data (DataStorage): DataStorage object.
        cpd_input (list): List of compounds selected by user.
        sample_filtering_enabled (bool): Enable the filtering by sample
        sample_filter_mode (str): Either Inlcude or Exclude selected by user.
        sample_filter_value (list): List of sample selected in filter.

    Returns:
        _type_: _description_
    """
    cscope_df, iscope_df, added_value_df = data.get_added_value_dataframe(cpd_input, sample_filtering_enabled, sample_filter_mode, sample_filter_value)

    added_value_fig = go.Figure(data=go.Heatmap(df_to_plotly(added_value_df), coloraxis="coloraxis"))

    cscope_fig = go.Figure(data=go.Heatmap(df_to_plotly(cscope_df), coloraxis="coloraxis"))

    iscope_fig = go.Figure(data=go.Heatmap(df_to_plotly(iscope_df), coloraxis="coloraxis"))

    return cscope_fig, iscope_fig, added_value_fig


def percentage_smpl_producing_cpd(data: DataStorage, cpd_input: list, metadata_filter_input: str):

    cscope_df = data.get_metabolite_production_dataframe()
    iscope_df = data.get_iscope_metabolite_production_dataframe()

    # Check if the iscope dataframe contain all the cpd in cscope dataframe. IF not add them as column filled with 0 value.
    col_diff = cscope_df.columns.difference(iscope_df.columns)

    col_diff_dict = dict.fromkeys(col_diff, 0.0)

    temp_df = pd.DataFrame(col_diff_dict, index=iscope_df.index)

    iscope_df = pd.concat([iscope_df, temp_df], axis=1)

    # Select only the column of intereset.
    cscope_df = cscope_df[["smplID", *cpd_input, metadata_filter_input]].dropna()
    iscope_df = iscope_df[["smplID", *cpd_input, metadata_filter_input]].dropna()

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

    fig1 = px.bar(cscope_df, x = "variable", y = "value", color = "metadata", barmode="group", title="Barplot showing the percentage of sample producing the compounds given in input (at least one bins of the sample).").update_layout(bargap=0.2,xaxis_title="Compounds", yaxis_title="Percent of sample producing the compound")

    fig2 = px.bar(iscope_df, x = "variable", y = "value", color = "metadata", barmode="group", title="Barplot showing the percentage of sample producing the compounds given in input (at least one bins of the sample).").update_layout(bargap=0.2,xaxis_title="Compounds", yaxis_title="Percent of sample producing the compound")

    return fig1, fig2


def col_value_to_percent(col: pd.Series):

    sum_val = col.sum()

    len_val = len(col.values)

    final_val = (sum_val / len_val) * 100

    return final_val
