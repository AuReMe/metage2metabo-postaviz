from m2m_postaviz.data_struct import DataStorage
import time
import pandas as pd
import plotly.express as px

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
    if factor != "None" and len(factor_choice) > 0:
        filter_condition.append((factor, "in", factor_choice))

    df = data.get_bin_dataframe(condition=filter_condition)
    
    unique_sample_in_df = df["smplID"].unique()

    new_serie_production = pd.DataFrame(columns=["smplID", "unique_production_count"])

    for sample in unique_sample_in_df:
        
        tmp_df = df.loc[df["smplID"] == sample][['binID','smplID','Production']]
        all_production = tmp_df["Production"].values

        tmp_production = []

        for prod_list in all_production:
            
            tmp_production += list(prod_list)

        unique_production_count = len(set(tmp_production))
        new_serie_production.loc[len(new_serie_production)] = {"smplID": sample, "unique_production_count": unique_production_count}

    df = df.merge(new_serie_production, how='inner', on="smplID")

    df.sort_index(inplace=True)

    fig1 = px.histogram(df, x="smplID", y="Count_with_abundance" if with_abundance else "unique_production_count", color="smplID" if color =="None" else color, hover_data="binID")

    # If only one bin selected do not make boxplot.
    
    if len(filtered_list_of_bin) > 1: 
        fig3 = px.box(df, x="smplID", y="Count_with_abundance" if with_abundance else "Count", color="smplID" if color =="None" else color, hover_data="binID")
        
    else:
        fig3 = None

    fig2 = px.bar(df, x="smplID", y="Abundance", color="Abundance", hover_data="binID")

    return fig1, fig2, df, time.time() - start_timer, fig3
