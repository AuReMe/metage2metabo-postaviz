import plotly.express as px
import plotly.graph_objects as go
from shiny import App
from shiny import render
from shiny import run_app
from shiny import ui
from shiny import reactive
from shinywidgets import output_widget
from shinywidgets import render_widget
import warnings
import numpy as np

import m2m_postaviz.data_utils as du
from m2m_postaviz.data_struct import DataStorage


def del_list_duplicate(mylist: list):
    return list(dict.fromkeys(mylist))


def run_shiny(data: DataStorage):
    ###
    warnings.filterwarnings("ignore", category=FutureWarning, module="plotly.express")
    list_of_cpd = data.get_cpd_list()

    producer_data = data.get_producer_long_dataframe()
    total_production = data.get_total_production_by_sample()

    if data.HAS_TAXONOMIC_DATA:
        long_taxo_df = data.get_long_taxonomic_data()

    factor_list = data.list_of_factor
    factor_list.insert(0, "None")

    pcoa_dataframe = data.current_pcoadf

    metadata_label = data.get_metadata_label()

    all_dataframe = {"producer_test_df": None, "producer_plot_df": None, "abundance_test_df": None, "abundance_plot_df": None}

    ### ALL CARD OBJECT TO BE ARRANGED ###

    abundance_boxplot =   ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("box_inputx1", "Label for X axis", factor_list),
                    ui.input_select("box_inputx2", "Label for 2nd X axis", factor_list),
                    ui.input_selectize("box_inputy1", "Compound for Y axis", list_of_cpd, multiple=True, selected=list_of_cpd[0]),
                    ui.input_checkbox("ab_norm", "With normalised data"),
                    ui.input_action_button("save_abundance_plot", "Save plot dataframe"),
                    ui.output_text_verbatim("save_ab_plot_txt", True),
                    ui.input_action_button("save_abundance_test", "Save test dataframe"),
                    ui.output_text_verbatim("save_ab_test_txt", True),
                    width=350,
                    gap=30,

            ),
        output_widget("Abundance_boxplot"), 
        ui.output_data_frame("abundance_test_dataframe")
        ),
        full_screen=True
        )

    if data.HAS_TAXONOMIC_DATA:
        taxonomy_boxplot = ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("tax_inpx1", "Label for X axis", factor_list),
                    ui.input_select("tax_inpx2", "Label for 2nd X axis", factor_list),
                    ui.input_selectize("tax_inpy1", "Taxa for Y axis", long_taxo_df["Taxa"].unique().tolist(), multiple=True),
                    ui.input_checkbox("taxo_norm", "With normalised data"),
                    width=350,
                ),
                output_widget("taxonomic_boxplot"),
            ),
            full_screen=True
        )
    else:
        taxonomy_boxplot = ui.output_text_verbatim("no_taxonomy", True),

    ### PRODUCER BOXPLOT CARD
    producer_boxplot = ui.card(
        ui.card_header("Total production of all compound, weighted with the abundance if provided."),
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select("prod_inputx1", "Label for X axis", factor_list),
                ui.input_select("prod_inputx2", "Label for 2nd X axis", factor_list),
                ui.input_checkbox("prod_norm", "Abundance data"),
                ui.input_action_button("save_producer_plot", "Save plot dataframe"),
                ui.output_text_verbatim("save_prod_plot_txt", True),
                ui.input_action_button("save_producer_test", "Save test dataframe"),
                ui.output_text_verbatim("save_prod_test_txt", True),
                width=350,
                gap=30,
            ),
            output_widget("producer_boxplot"),
            ui.output_data_frame("production_test_dataframe")
        
        ),
        full_screen=True
        )

    metadata_table = ui.card(
        ui.row(
                ui.input_select("metadata_factor", "Current column: ", metadata_label, selected=metadata_label[0]),
                ui.input_select("metadata_dtype", "dtype: ", ["category", "str", "int", "float"]),
                ui.input_action_button("dtype_change", "Update")
               ),
        ui.output_text_verbatim("update_metadata_log", True),
        ui.output_data_frame("metadata_table")
        )

    pcoa_plot_dev_table = ui.card(
        ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_select(id="pcoa_color", label="Plot color.", choices=metadata_label, selected=metadata_label[0]),
                        ui.output_ui("pcoa_ui"),
                        width=350,
                        gap=30,
                    ),
        output_widget("pcoa_plot")
        ),
        full_screen=True
    )

    ### APPLICATION TREE ###

    app_ui = ui.page_fillable(
        ui.navset_tab(
            ui.nav("Exploration", producer_boxplot, abundance_boxplot, taxonomy_boxplot),
            ui.nav(
                "Metadata",
                metadata_table,
                # dev_table
            ),
            ui.nav(
                "PCOA",
                    pcoa_plot_dev_table,
                ),
            ),
        )

    def server(input, output, session):

        @render_widget
        def pcoa_plot():

            # Get all parameters.
            selected_col = input.pcoa_color()

            df = pcoa_dataframe

            value = df[selected_col]

            plot_size = None

            min_limit = input.pcoa_float_limit()[0]

            max_limit = input.pcoa_float_limit()[1]

            pcoa_symbol = input.pcoa_symbol()

            pcoa_size = input.pcoa_size()
            # Check column dtype.
            if value.dtype == 'float64' or value.dtype == 'int64':

                df = df.loc[(df[selected_col] <= max_limit) & (df[selected_col] >= min_limit)]
            
            if df[pcoa_size].dtype == 'float64' or df[pcoa_size].dtype == 'int64':

                plot_size = pcoa_size
                

            return px.scatter(df, x="PC1", y="PC2", color=selected_col, symbol=pcoa_symbol, size=plot_size)

        @render.ui
        @reactive.event(input.pcoa_color)
        def pcoa_ui():
            value = pcoa_dataframe[input.pcoa_color()]
            if value.dtype == 'float64' or value.dtype == 'int64':
                return ui.TagList(
                            ui.input_select("pcoa_symbol", "Symbol :",metadata_label,selected=metadata_label[0]),
                            ui.input_select("pcoa_size", "Size :",metadata_label,selected=metadata_label[0]),
                            ui.input_slider("pcoa_float_limit", "Show only :", min=value.min(), max=value.max(), value=[value.min(),value.max()]),
                        )
            else:
                return ui.TagList(
                            ui.input_select("pcoa_symbol", "Symbol :",metadata_label,selected=metadata_label[0]),
                            ui.input_select("pcoa_size", "Size :",metadata_label,selected=metadata_label[0])
                        )
        
        @render_widget
        def taxonomic_boxplot():
            if not data.HAS_TAXONOMIC_DATA:
                return
            df = long_taxo_df
            x1, x2, y1 = input.tax_inpx1(), input.tax_inpx2(), input.tax_inpy1()
            if len(y1) == 0:
                y1 = df["Taxa"].unique()
            if x1 == "None":
                return px.box(df, y="Nb_taxon")
            if x2 == "None":
                fig = px.box(
                    df,
                    x=x1,
                    y="Nb_taxon",
                    color=x1,
                )
                du.add_p_value_annotation(fig, [[0, 1]])
                return fig
            else:
                conditionx2 = df[x2].unique()

                fig = go.Figure()

                for condition in conditionx2:
                    fig.add_trace(
                        go.Box(
                            x=df.loc[df[x2] == condition][x1],
                            y=df["Nb_taxon"],
                            name=str(condition),
                        )
                    )

                fig.update_layout(boxmode="group")
                return fig

        @render.data_frame
        def abundance_test_dataframe():
            # Get input
            with_normalised_data = input.ab_norm()
            if with_normalised_data and data.HAS_ABUNDANCE_DATA:
                df = data.get_melted_norm_ab_dataframe()
                column_value = "Quantity"
            else:
                df = producer_data
                column_value = "Value"

            y1, x1, x2 = input.box_inputy1(), input.box_inputx1(), input.box_inputx2()
            if len(y1) == 0:
                return
            df = df.loc[df["Compound"].isin(y1)]
            # No input selected
            if x1 == "None":
                return 
            # At least first axis selected
            else:
                if x2 == "None":
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        selected_data = df.loc[df[x1] == layer_1][column_value]
                        if len(selected_data) > 1:
                            tested_data[layer_1] = {}
                            tested_data[layer_1]["data"] = selected_data
                            tested_data[layer_1]["n_data"] = len(selected_data)
                    res = du.stat_on_plot(tested_data, 1)
                    all_dataframe["abundance_test_df"] = res
                    return res
                # Both axis have been selected
                else:
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        tested_data[layer_1] = {}
                        for layer_2 in df[x2].unique():
                            selected_data = df.loc[(df[x1] == layer_1) & (df[x2] == layer_2)][column_value]
                            if len(selected_data) > 1:
                                tested_data[layer_1][layer_2] = {}
                                tested_data[layer_1][layer_2]["data"] = selected_data
                                tested_data[layer_1][layer_2]["n_data"] = len(selected_data)
                    res = du.stat_on_plot(tested_data, 2)
                    all_dataframe["abundance_test_df"] = res
                    return res

        @render.data_frame()
        def production_test_dataframe():
            # Get input
            with_normalised_data = input.prod_norm()
            if with_normalised_data and data.HAS_ABUNDANCE_DATA:
                column_value = "Total_abundance_weighted"
            else:
                column_value = "Total_production"
            df = total_production
            x1, x2 = input.prod_inputx1(), input.prod_inputx2()
            # No input selected
            if x1 == "None":
                return 
            # At least first axis selected
            else:
                if x2 == "None":
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        selected_data = df.loc[df[x1] == layer_1][column_value]
                        if len(selected_data) > 1:
                            tested_data[layer_1] = {}
                            tested_data[layer_1]["data"] = selected_data
                            tested_data[layer_1]["n_data"] = len(selected_data)
                    res = du.stat_on_plot(tested_data, 1)
                    all_dataframe["producer_test_df"] = res
                    return res
                # Both axis have been selected
                else:
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        tested_data[layer_1] = {}
                        for layer_2 in df[x2].unique():
                            selected_data = df.loc[(df[x1] == layer_1) & (df[x2] == layer_2)][column_value]
                            if len(selected_data) > 1:
                                tested_data[layer_1][layer_2] = {}
                                tested_data[layer_1][layer_2]["data"] = selected_data.values
                                tested_data[layer_1][layer_2]["n_data"] = len(selected_data)
                    # In case tested_data has no value (barplot)
                    if all(len(tested_data[k]) == 0 for k in tested_data.keys()):
                        print("No values in each keys")
                        return
                    # In case tested_data is not a double layer dict
                    if all(len(tested_data[k]) == 1 for k in tested_data.keys()):
                        print("Only one pair per layer.")
                        return
                    else:
                        res = du.stat_on_plot(tested_data, 2)
                    all_dataframe["producer_test_df"] = res
                    return res
        
        @render.text
        @reactive.event(input.save_producer_test)
        def save_prod_test_txt():
            if bool(all_dataframe):
                if all_dataframe["producer_test_df"] is not None:
                   log = data.save_dataframe(all_dataframe["producer_test_df"], "producer_test_dataframe")
                else:
                    log = "Unable to find any dataframe to save."
            return log
        
        @render.text
        @reactive.event(input.save_abundance_test)
        def save_ab_test_txt():
            if bool(all_dataframe):
                if all_dataframe["abundance_test_df"] is not None:
                   log = data.save_dataframe(all_dataframe["abundance_test_df"], "abundance_test_dataframe")
                else:
                    log = "Unable to find any dataframe to save."
            return log

        @render.text
        @reactive.event(input.save_abundance_plot)
        def save_ab_plot_txt():
            if bool(all_dataframe):
                if all_dataframe["abundance_plot_df"] is not None:
                   log = data.save_dataframe(all_dataframe["abundance_plot_df"], "abundance_plot_dataframe")
                else:
                    log = "Unable to find any dataframe to save."
            return log

        @render.text
        @reactive.event(input.save_producer_plot)
        def save_prod_plot_txt():
            if bool(all_dataframe):
                if all_dataframe["producer_plot_df"] is not None:
                    log = data.save_dataframe(all_dataframe["producer_plot_df"], "producer_plot_dataframe")
                else:
                    log = "Unable to find any dataframe to save."
            return log

        @render_widget
        def Abundance_boxplot():
            # Which type of dataframe
            with_abundance_data = input.ab_norm()
            if with_abundance_data and data.HAS_ABUNDANCE_DATA:
                df = data.get_melted_norm_ab_dataframe()
                column_value = "Quantity"
            else:
                df = producer_data
                column_value = "Value"
                
            all_dataframe["abundance_plot_df"] = df
            # import button input from shiny
            y1, x1, x2 = input.box_inputy1(), input.box_inputx1(), input.box_inputx2()
            df = df.loc[df["Compound"].isin(y1)]

            if len(y1) == 0:
                return
            if x1 == "None":
                return px.box(df, y=df.loc[df["Compound"].isin(y1)][column_value], title=f"Estimated amount of {*y1,} produced by all sample.")
            
            elif x2 == "None":
                # df.sort_values(x1)
                if df.shape[0] == len(df[x1].unique()):
                    return px.bar(df, x=x1 , y=column_value, color=x1, title=f"Estimated amount of {*y1,} produced by {x1}.")
                else:
                    return px.box(df, x=x1 , y=column_value, color=x1, title=f"Estimated amount of {*y1,} produced by {x1}.")

            else:
                # df.sort_values([x1,x2])
                conditionx2 = df[x2].unique()
                conditionx1 = df[x1].unique()

                fig = go.Figure()

                has_unique_value = True

                # If one of the case has multiple y values --> boxplot
                for layer1 in conditionx1:
                    for layer2 in conditionx2:
                        if len(df.loc[(df[x2] == layer2) & (df[x1] == layer1)][column_value]) > 1:
                            has_unique_value = False

                for condition2 in conditionx2:
                    if has_unique_value:
                        fig.add_trace(
                        go.Bar(
                            x=df.loc[df[x2] == condition2][x1],
                            y=df.loc[(df["Compound"].isin(y1)) & (df[x2] == condition2)][column_value],
                            name=str(condition2),
                        )
                        )

                    else:
                        fig.add_trace(
                            go.Box(
                                x=df.loc[df[x2] == condition2][x1],
                                y=df.loc[(df["Compound"].isin(y1)) & (df[x2] == condition2)][column_value],
                                name=str(condition2),
                            )
                        )
                    if has_unique_value:
                        fig.update_layout(
                            barmode="group",
                            hovermode="x",
                            bargroupgap=0.2,
                            title=(f"Estimated quantity of {*y1,} produced by {x1} and {x2}"),
                        )
                    else:
                        fig.update_layout(
                            boxmode="group",
                            hovermode="x",
                            boxgroupgap=0.1,
                            title=(f"Estimated quantity of {*y1,} produced by {x1} and {x2}"),
                        )
                return fig

        @render_widget
        def producer_boxplot():
            with_normalised_data = input.prod_norm()
            if with_normalised_data and data.HAS_ABUNDANCE_DATA:
                column_value = "Total_abundance_weighted"
            else:
                column_value = "Total_production"

            df = total_production
            all_dataframe["producer_plot_df"] = df

            inputx1 , inputx2 = input.prod_inputx1(), input.prod_inputx2()
        
            if inputx1 == "None":
                return px.bar(df, x="smplID", y=column_value, title=f"Total production.", color="smplID")
            elif inputx2 == "None":
                if df.shape[0] == len(df[inputx1].unique()):
                    return px.bar(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
                else:
                    return px.box(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
            else:
                fig = go.Figure()

                has_unique_value = True

                # If one of the case has one y values --> barplot else --> boxplot
                for layer1 in df[inputx1].unique():
                    for layer2 in df[inputx2].unique():   
                        if len(df.loc[(df[inputx2] == layer2) & (df[inputx1] == layer1)][column_value]) > 1:
                            has_unique_value = False

                for x2 in df[inputx2].unique():
                    if has_unique_value:
                        fig.add_trace(
                        go.Bar(x=df.loc[df[inputx2] == x2][inputx1],
                                y=df.loc[df[inputx2] == x2][column_value],
                                name=str(x2),)
                        )
                    else:
                        fig.add_trace(
                        go.Box(x=df.loc[df[inputx2] == x2][inputx1],
                                y=df.loc[df[inputx2] == x2][column_value],
                                name=str(x2),)
                        )
                    if has_unique_value:
                        fig.update_layout(
                                barmode="group",
                                hovermode="x",
                                bargroupgap=0.2,
                                title=(f"Total compound production filtered by {inputx1} and {inputx2}"),
                            )
                    else:
                        fig.update_layout(
                                boxmode="group",
                                hovermode="x",
                                boxgroupgap=0.1,
                                title=(f"Total compound production filtered by {inputx1} and {inputx2}"),
                            )
                return fig

        @render.data_frame
        def dev_table():
            return render.DataGrid(data.get_long_taxonomic_data())

        @render.data_frame
        def metadata_table():
            df = data.get_main_metadata()
            return df
        
        @render.text()
        @reactive.event(input.dtype_change)
        def update_metadata_log():
            text = "No changes applied."
            factor_choice, dtype_choice = input.metadata_factor(), input.metadata_dtype()
            df = data.get_main_metadata()
            df_prod = producer_data
            df_tot = total_production
            df_pcoa = pcoa_dataframe
            try:
                df[factor_choice] = df[factor_choice].astype(dtype_choice)
                df_prod[factor_choice] = df_prod[factor_choice].astype(dtype_choice)
                df_tot[factor_choice] = df_tot[factor_choice].astype(dtype_choice)
                df_pcoa[factor_choice] = df_pcoa[factor_choice].astype(dtype_choice)
                text = f"Column {factor_choice} changed to {dtype_choice}."
                data.set_main_metadata(df)
            except ValueError as e:
                text = f"Cannot perform changes, {e}" 
            return text
        
        @render.text
        def no_taxonomy():
            return "No taxonomic data provided."

    app = App(app_ui, server)
    run_app(app=app, launch_browser=False)
