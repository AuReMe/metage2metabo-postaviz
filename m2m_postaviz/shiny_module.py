import time
import sys

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
    print(df)
    print(metadata)
    df = df.merge(metadata, "inner", "smplID")
    print(df[factor].isin(factor_choice))
    df = df.loc[df[factor].isin(factor_choice)]
    df.sort_index(inplace=True)
    print(df)


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


def metabolites_production_statistical_dataframe(data: DataStorage, metabolites_choices, user_input1, user_input2, multiple_test_correction, correction_method, with_abundance):
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


def run_pcoa(main_df: pd.DataFrame, metadata: pd.DataFrame, distance_method: str = "jaccard"):
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
    main_dataframe = main_df

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
        return px.box(df, y=column_value, title="Total production.")

    elif user_input2 == "None" or user_input1 == user_input2:

        df = df[[column_value,user_input1]]
        df = df.dropna()
        # all_dataframe["global_production_plot_dataframe"] = df

        if du.has_only_unique_value(df, user_input1):

            return px.bar(df, x=user_input1 , y=column_value, color=user_input1, title=f"Total compound production filtered by {user_input1}")

        else:

            fig = px.box(df, x=user_input1 , y=column_value, color=user_input1, title=f"Total compound production filtered by {user_input1}")
            return fig

    else:

        df = df[[column_value,user_input1,user_input2]]
        df = df.dropna()
        # all_dataframe["global_production_plot_dataframe"] = df
        has_unique_value = du.has_only_unique_value(df, user_input1, user_input2)

        return px.bar(df,x=user_input1,y=column_value,color=user_input2) if has_unique_value else px.box(df,x=user_input1,y=column_value,color=user_input2)


def render_reactive_metabolites_production_plot(data: DataStorage, compounds_input, user_input1, color_input, sample_filtering_enabled, sample_filter_button, sample_filter_value = [], with_abundance = None):

    if len(compounds_input) == 0:
        return

    producer_data = data.get_metabolite_production_dataframe()
    # producer_data_iscope = data.get_iscope_metabolite_production_dataframe()

    if sample_filtering_enabled and len(sample_filter_value) != 0:

        if sample_filter_button is "Include":

            producer_data = producer_data.loc[producer_data["smplID"].isin(sample_filter_value)]

        elif sample_filter_button is "Exclude":

            producer_data = producer_data.loc[~producer_data["smplID"].isin(sample_filter_value)]

    # producer_data = producer_data.set_index("smplID")

    if user_input1 == "None":

        if color_input == "None":

            df = producer_data[[*compounds_input]]
            return px.box(df, y=compounds_input)

        else:

            df = producer_data[[*compounds_input, color_input]]
            return px.box(df, y=compounds_input, color=color_input)

    if color_input == "None" or user_input1 == color_input:

        df = producer_data[[*compounds_input,user_input1]]
        df = df.dropna()

        if df[user_input1].dtypes == float:

            df[user_input1] = df[user_input1].astype(str)

        has_unique_value = du.has_only_unique_value(df, user_input1)

        return px.bar(df, x=user_input1, y=compounds_input, color=user_input1) if has_unique_value else px.box(df, x=user_input1, y=compounds_input, color=user_input1)

    df = producer_data[[*compounds_input,user_input1,color_input]]
    df = df.dropna()

    if df[user_input1].dtypes == float:

        df[user_input1] = df[user_input1].astype(str)

    has_unique_value = du.has_only_unique_value(df, user_input1, color_input)

    return px.bar(df, x=user_input1, y=compounds_input,color=color_input) if has_unique_value else px.box(df, x=user_input1, y=compounds_input, color=color_input, boxmode="group")


def metabolites_heatmap(data: DataStorage, user_cpd_list, by, value):

    all_cscope_dataframe = []
    all_iscope_dataframe = []

    for cpd in user_cpd_list:
        
        try:

            tmp_df = data.get_minimal_cpd_dataframe(cpd)
            all_cscope_dataframe.append(tmp_df[0])
            all_iscope_dataframe.append(tmp_df[1])

        except Exception as e:

            print(f'{cpd} not in icscope dataframe')

    dfc = pd.concat(all_cscope_dataframe,axis=0)
    dfi = pd.concat(all_iscope_dataframe,axis=0)

    if by == "taxonomy":

        taxonomy_file = data.get_taxonomic_dataframe()

        mgs_col_taxonomy = taxonomy_file.columns[0]

        dfc["binID"] = dfc["binID"].map(taxonomy_file.set_index(mgs_col_taxonomy)[value])

        dfi["binID"] = dfi["binID"].map(taxonomy_file.set_index(mgs_col_taxonomy)[value])

    if by == "metadata":

        metadata = data.get_metadata()

        dfc["smplID"] = dfc["smplID"].map(metadata.set_index("smplID")[value])

        dfi["smplID"] = dfi["smplID"].map(metadata.set_index("smplID")[value])
    
    dfc.fillna(0,inplace=True)

    dfc.set_index(["binID","smplID"],inplace=True)
    dfc = dfc.groupby(level="binID").sum()

    dfi.fillna(0,inplace=True)

    dfi.set_index(["binID","smplID"],inplace=True)
    dfi = dfi.groupby(level="binID").sum()

    print(sys.getsizeof(dfc) / 1000000000, "Go")

    print(sys.getsizeof(dfi) / 1000000000, "Go")

    fig = make_subplots(1,2,x_title="Compounds",y_title="Bins")

    fig.add_trace(go.Heatmap(df_to_plotly(dfc), coloraxis='coloraxis', ),1,1)

    fig.add_trace(go.Heatmap(df_to_plotly(dfi), coloraxis='coloraxis', ),1,2)

    fig.update_layout(coloraxis=dict(colorscale = 'ylgnbu'))

    return fig


def df_to_plotly(df):
    return {'z': df.values.tolist(),
            'x': df.columns.tolist(),
            'y': df.index.tolist()}
