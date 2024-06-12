import plotly.express as px
import plotly.graph_objects as go
from shiny import App
from shiny import render
from shiny import run_app
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget
import warnings

import m2m_postaviz.data_utils as du
from m2m_postaviz.data_struct import DataStorage


def del_list_duplicate(mylist: list):
    return list(dict.fromkeys(mylist))


def run_shiny(data: DataStorage):
    ###
    warnings.filterwarnings("ignore", category=FutureWarning, module="plotly.express")
    producer_data = data.get_producer_long_dataframe()
    cpd_prod_by_sample = data.get_total_cpd_production_by_sample_and_compound()
    long_taxo_df = data.get_long_taxonomic_data()
    list_of_bin = data.get_bin_list()
    list_of_cpd = data.get_cpd_list()
    total_production = data.get_total_production_by_sample()

    metadata = data.get_main_metadata()
    main_dataframe = data.get_main_dataframe()
    current_dataframe = data.get_main_metadata()

    factor_list = data.list_of_factor
    factor_list.insert(0, "None")

    metadata_label = data.get_metadata_label()

    ### ALL CARD OBJECT TO BE ARRANGED ###

    abundance_boxplot =   ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("box_inputx1", "Label for X axis", factor_list),
                    ui.input_select("box_inputx2", "Label for 2nd X axis", factor_list),
                    ui.input_selectize("box_inputy1", "Compound for Y axis", list_of_cpd, multiple=True, selected=list_of_cpd[0]),
                    ui.input_checkbox("ab_norm", "With normalised data"),
            ),
        output_widget("Abundance_boxplot"), 
        ui.output_data_frame("abundance_test_dataframe")
        ))

    taxonomy_boxplot = ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select("tax_inpx1", "Label for X axis", factor_list),
                ui.input_select("tax_inpx2", "Label for 2nd X axis", factor_list),
                ui.input_selectize("tax_inpy1", "Taxa for Y axis", long_taxo_df["Taxa"].unique().tolist(), multiple=True),
                ui.input_checkbox("taxo_norm", "With normalised data"),
            ),
            output_widget("taxonomic_boxplot"),
        )
    )
    ### PRODUCER BOXPLOT CARD
    producer_boxplot = ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select("prod_inputx1", "Label for X axis", factor_list),
                ui.input_select("prod_inputx2", "Label for 2nd X axis", factor_list),
                ui.input_checkbox("prod_norm", "Normalised data"),
            ),
            output_widget("producer_boxplot"),
            ui.output_data_frame("production_test_dataframe")
        )
    )

    main_table = ui.card(ui.output_data_frame("dev_table"))  # , ui.output_data_frame("main_table"))
    summary_table = ui.card(ui.output_data_frame("summary_table"))

    main_panel_dataframe = ui.card(
        ui.row(
            ui.input_select(id="factor_choice", label="Filter", choices=factor_list, selected=factor_list[0]),
            ui.input_select(id="factor_choice2", label="Value", choices=[]),
        ),
        ui.output_data_frame("main_panel_table"),
    )
    pcoa_plot_dev_table = ui.card(output_widget("pcoa_plot"))

    ### APPLICATION TREE ###

    app_ui = ui.page_fillable(
        ui.navset_tab(
            ui.nav("Exploration", producer_boxplot, abundance_boxplot, taxonomy_boxplot),
            ui.nav(
                "Dev mod",
                main_table,
                # dev_table
            ),
            ui.nav(
                "PCOA",
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_select(id="pcoa_color", label="Plot color.", choices=metadata_label, selected=metadata_label[0]),
                        ui.input_checkbox("pcoa_check", "with abundance data."),
                    ),
                    pcoa_plot_dev_table,
                ),
            ),
        )
    )

    def server(input, output, session):
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
                column_value = "Nb_producers"

            y1, x1, x2 = input.box_inputy1(), input.box_inputx1(), input.box_inputx2()
            if len(y1) == 0:
                return df

            # No input selected
            if x1 == "None":
                return df.loc[df["Compound"].isin(y1)]
            # At least first axis selected
            else:
                if x2 == "None":
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        tested_data[layer_1] = {}
                        tested_data[layer_1]["data"] = df.loc[(df[x1] == layer_1) & (df["Compound"].isin(y1))][column_value]
                        tested_data[layer_1]["n_data"] = len(df.loc[(df[x1] == layer_1) & (df["Compound"].isin(y1))][column_value])
                    return du.stat_on_plot(tested_data, 1)
                # Both axis have been selected
                else:
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        tested_data[layer_1] = {}
                        for layer_2 in df[x2].unique():
                            selected_data = df.loc[(df[x1] == layer_1) & (df[x2] == layer_2) & (df["Compound"].isin(y1))][column_value]
                            if len(selected_data) != 0:
                                tested_data[layer_1][layer_2] = {}
                                tested_data[layer_1][layer_2]["data"] = selected_data
                                tested_data[layer_1][layer_2]["n_data"] = len(selected_data)
                    return du.stat_on_plot(tested_data, 2)

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
                        tested_data[layer_1] = {}
                        tested_data[layer_1]["data"] = df.loc[df[x1] == layer_1][column_value]
                        tested_data[layer_1]["n_data"] = len(df.loc[df[x1] == layer_1][column_value])
                    return du.stat_on_plot(tested_data, 1)
                # Both axis have been selected
                else:
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        tested_data[layer_1] = {}
                        for layer_2 in df[x2].unique():
                            selected_data = df.loc[(df[x1] == layer_1) & (df[x2] == layer_2)][column_value]
                            if len(selected_data) != 0:
                                tested_data[layer_1][layer_2] = {}
                                tested_data[layer_1][layer_2]["data"] = selected_data.values
                                tested_data[layer_1][layer_2]["n_data"] = len(selected_data)
                    if x1 == "Patient" and x2 == "Days":
                        print(tested_data)
                    return du.stat_on_plot(tested_data, 2)
                
        @render_widget
        def Abundance_boxplot():
            # Which type of dataframe
            with_abundance_data = input.ab_norm()
            if with_abundance_data and data.HAS_ABUNDANCE_DATA:
                df = data.get_melted_norm_ab_dataframe()
                column_value = "Quantity"
            else:
                df = producer_data
                column_value = "Nb_producers"
                
            # import button input from shiny
            y1, x1, x2 = input.box_inputy1(), input.box_inputx1(), input.box_inputx2()

            if len(y1) == 0:
                return
            if x1 == "None":
                return px.box(
                    df, y=df.loc[df["Compound"].isin(y1)][column_value], title=f"Estimated amount of {*y1,} produced by all sample."
                )
            else:
                if x2 == "None":
                    new_df = df.loc[df["Compound"].isin(y1)]
                    fig = px.box(
                        new_df,
                        x=x1,
                        y=column_value,
                        color=x1,
                        title=f"Estimated amount of {*y1,} produced by {x1}.",
                    )
                    return fig
                else:
                    conditionx2 = df[x2].unique()
                    conditionx1 = df[x1].unique()

                    fig = go.Figure()

                    has_unique_value = True

                    # If one of the case has multiple y values --> boxplot
                    for layer1 in conditionx1:
                        for layer2 in conditionx2:
                            if not len(df.loc[(df[x2] == layer2) & (df[x1] == layer1) & (df["Compound"].isin(y1))][column_value]) <= 1:
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

            inputx1 , inputx2 = input.prod_inputx1(), input.prod_inputx2()
        
            if inputx1 == "None":
                return px.bar(df, x="smplID", y=column_value, title=f"Total production.")
            elif inputx2 == "None":
                if df.shape[0] == len(df[inputx1].unique()):
                    return px.bar(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
                else:
                    return px.box(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
            else:
                fig = go.Figure()

                has_unique_value = False

                # If one of the case has one y values --> barplot else --> boxplot
                for layer1 in df[inputx1].unique():
                    for layer2 in df[inputx2].unique():   
                        if not len(df.loc[(df[inputx2] == layer2) & (df[inputx1] == layer1)][column_value]) == 1:
                            has_unique_value = True

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

        @render_widget
        def pcoa_plot():
            input_ab = input.pcoa_check()
            input_color = input.pcoa_color()
            if input_ab:
                in_df = data.get_ab_dataframe().drop(["Days", "Group", "Antibiotics", "Patient"], axis=1)
                res = data.run_pcoa(in_df, metadata, "braycurtis")
            else:
                res = data.run_pcoa(main_dataframe, metadata)
            return px.scatter(res.samples, x="PC1", y="PC2", color=input_color)

        @render.data_frame
        def dev_table():
            return render.DataGrid(data.get_long_taxonomic_data())

        @render.data_frame
        def main_table():
            return render.DataGrid(data.abundance_matrix)

    app = App(app_ui, server)
    run_app(app=app, launch_browser=False, reload_dirs="./")
