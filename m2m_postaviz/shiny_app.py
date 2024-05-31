import plotly.express as px
import plotly.graph_objects as go
import scipy.stats as stats
from shiny import App
from shiny import reactive
from shiny import render
from shiny import run_app
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.data_utils as du
from m2m_postaviz.data_struct import DataStorage


def del_list_duplicate(mylist: list):
    return list(dict.fromkeys(mylist))


def run_shiny(data: DataStorage):
    ###
    producer_data = data.get_producer_long_dataframe()
    taxonomic_data = data.get_taxonomic_data()
    long_taxo_df = data.get_long_taxonomic_data()
    list_of_bin = data.get_bin_list()
    list_of_cpd = data.get_cpd_list()

    metadata = data.get_main_metadata()
    main_dataframe = data.get_main_dataframe()
    current_dataframe = data.get_main_metadata()

    factor_list = data.list_of_factor
    factor_list.insert(0, "None")

    metadata_label = data.get_metadata_label()

    ### ALL CARD OBJECT TO BE ARRANGED ###

    abundance_input = ui.row(
        ui.input_select("box_inputx1", "Label for X axis", factor_list),
        ui.input_select("box_inputx2", "Label for 2nd X axis", factor_list),
        ui.input_selectize("box_inputy1", "Compound for Y axis", list_of_cpd, multiple=True, selected=list_of_cpd[0]),
        ui.input_checkbox("ab_norm", "With normalised data")
        # ui.input_action_button("abundance_boxplot_button","Launch")
    )
    abundance_boxplot = ui.card(output_widget("Abundance_boxplot", height="100%", width="100%"),
                                ui.output_data_frame("Abundance_sig_df"))
    ### TAXONOMY BOXPLOT CARD W/I SHINYWIDGET WITH PLOTLY ONLY. (NO HEIGHT OR WIDTH CHANGES POSSIBLE)
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
                ui.input_select("prod_i1", "First X axis", factor_list),
                ui.input_select("prod_i2", "Second X axis", factor_list),
                ui.input_selectize("prod_i3", "Choose compounds", list_of_cpd, multiple=True),
                ui.input_checkbox("with_iscope", "Show iscope data."),
            ),
            output_widget("producer_boxplot"),
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
            ui.nav("Exploration", ui.layout_sidebar(ui.sidebar(abundance_input), abundance_boxplot), taxonomy_boxplot,producer_boxplot),
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
            ui.nav("Taxonomy"),
            ui.nav(
                "Main",
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_select(
                            id="stat_test_choice",
                            label="Statistical test available.",
                            choices=["Choose test", "Student", "Wilcoxon", "ANOVA", "PCOA"],
                        ),
                        ui.panel_conditional(
                            "input.stat_test_choice === 'Wilcoxon' || input.stat_test_choice === 'Student'",
                            ui.input_select(
                                id="wilcoxon_variable_select", label="Select a variable", choices=list(current_dataframe.columns)
                            ),
                            ui.input_select(id="wilcoxon_group1_select", label="Group1 choice", choices=[]),
                            ui.input_select(id="wilcoxon_group2_select", label="Group2 choice", choices=[]),
                            ui.input_select(id="wilcoxon_tested_value_select", label="Tested value choice", choices=[]),
                        ),
                        ui.input_action_button("run_test", "run_test"),
                    ),
                    main_panel_dataframe,
                    ui.output_text(id="stat_test_result"),
                ),
            ),
        )
    )

    def server(input, output, session):
        @render_widget
        def taxonomic_boxplot():
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
        def Abundance_sig_df():
            # Get input
            with_normalised_data = input.ab_norm()
            df = data.get_melted_norm_ab_dataframe()

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
                        tested_data[layer_1] = df.loc[(df[x1] == layer_1) & (df["Compound"].isin(y1))]['Quantity']
                    return du.stat_on_plot(tested_data, 1)
                # Both axis have been selected
                else:
                    tested_data = {}
                    for layer_1 in df[x1].unique():
                        tested_data[layer_1] = {}
                        for layer_2 in df[x2].unique():
                            selected_data = df.loc[(df[x1] == layer_1) & (df[x2] == layer_2) & (df["Compound"].isin(y1))]['Quantity']
                            if len(selected_data) != 0:
                                identifier = str(layer_2)+" (n="+str(len(selected_data))+")"
                                tested_data[layer_1][identifier] = selected_data
                    return du.stat_on_plot(tested_data, 2)


        @render_widget
        def Abundance_boxplot():
            # Which type of dataframe
            with_abundance_data = input.ab_norm()
            if with_abundance_data:
                df = data.get_melted_norm_ab_dataframe()
                variable_of_interest = "Quantity"
            else:
                df = producer_data
                variable_of_interest = "Nb_producers"
            # import button input from shiny
            y1 = input.box_inputy1()
            x1 = input.box_inputx1()
            x2 = input.box_inputx2()
            # If none selected, all compound selected (default)
            if len(y1) == 0:
                return
            if x1 == "None":
                return px.box(df, y=df.loc[df["Compound"].isin(y1)][variable_of_interest],title=f"Estimated amount of {*y1,} produced by all sample.")
            else:
                if x2 == "None":
                    new_df = df.loc[df["Compound"].isin(y1)]
                    fig = px.box(new_df, x=x1, y=variable_of_interest, color=x1, title=f"Estimated amount of {*y1,} produced by {x1}.",)

                    return fig
                else:
                    conditionx2 = df[x2].unique()
                    conditionx1 = df[x1].unique()

                    fig = go.Figure()

                    annotation_list = []
                    for index in range(len(conditionx1)):
                        if not index == 0:
                            annotation_list.append([index - 1, index])

                    for condition2 in conditionx2:

                        fig.add_trace(
                            go.Box(
                                x=df.loc[df[x2] == condition2][x1],
                                y=df.loc[(df["Compound"].isin(y1)) & (df[x2] == condition2)][variable_of_interest],
                                name=str(condition2),
                            )
                        )
                        fig.update_layout(
                            boxmode="group",
                            hovermode="x",
                            boxgroupgap=0.1,
                            title=(f"Estimated quantity of {*y1,} produced by {x1} and {x2}"),
                        )
                    return fig

        @render_widget
        def producer_boxplot():
            df = producer_data
            print(df)
            input1, input2, cpd_choice = input.prod_i1(), input.prod_i2(), input.prod_i3()
            with_iscope = input.with_iscope()
            if len(cpd_choice) == 0:
                cpd_choice = list_of_cpd
            fig = go.Figure()
            if input1 == "None":
                fig.add_trace(go.Box(
                    # x=df["smplID"],
                    y=df["Cscope_cpd_prod"],
                    name="Cscope production",
                    marker_color="#ffc82c",
                    )
                )
                if with_iscope:
                    fig.add_trace(go.Box(
                        # x=df["smplID"],
                        y=df["Iscope_cpd_prod"],
                        name="Iscope production",
                        marker_color="#173a64",
                        )
                    )
                fig.update_layout(barmode="group",bargroupgap=0,bargap=0.5,title="Compound production by community and individual.")
                return fig
            else:
                lenght_of_grid = len(df[input1].unique())
                if input2 == "None":
                    
                    fig = go.Figure()
                    fig.add_trace(go.Box(
                        x=df[input1],
                        y=df["Cscope_cpd_prod"],
                        name="Cscope production",
                        marker_color="#ffc82c",
                        )
                    )
                    if with_iscope:
                        fig.add_trace(go.Box(
                            x=df[input1],
                            y=df["Iscope_cpd_prod"],
                            name="Iscope production",
                            marker_color="#173a64",
                            )
                        )
                    fig.update_layout(boxmode="group",
                                      boxgroupgap=0,boxgap=0.5,
                                      title=f"Compound production by {input1}."
                                      ,xaxis_title=input1
                                      )
                    return fig
                else:
                    fig = go.Figure()
                    # tmp_df = df.loc[df[input2].isin()]
                    for condition in df[input2].unique():

                        fig.add_trace(go.Box(
                            x=df[input1],
                            y=df.loc[df[input2] == condition]["Cscope_cpd_prod"],
                            name=f"Cscope production for {input2} {condition}",
                            # marker_color="#ffc82c",
                            )
                        )
                        if with_iscope:
                            fig.add_trace(go.Box(
                                x=df[input1],
                                y=df.loc[df[input2] == condition]["Iscope_cpd_prod"],
                                name=f"Iscope production for {input2} {condition}",
                                opacity=0.5,
                                marker_color="#173a64",
                                )
                            )
                    fig.update_layout(boxmode="group",boxgroupgap=0,boxgap=0.5,title=f"Compound production by {input1} and {input2}.",xaxis_title=input1,grid={"columns": lenght_of_grid})
                    return fig


            

        @reactive.Effect()
        @reactive.event(input.test_sample1)
        def process():
            current_factor_choice = del_list_duplicate(current_dataframe[input.test_sample1()])
            ui.update_select("test_sample2", choices=current_factor_choice)
            return

        @reactive.Effect()
        @reactive.event(input.factor_choice)
        def process():
            if input.factor_choice() == "None":
                return
            current_factor_choice = del_list_duplicate(current_dataframe[input.factor_choice()])
            # current_factor_choice.insert(0,"None")
            ui.update_select("factor_choice2", choices=current_factor_choice)
            return

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

        @render.data_frame
        def main_panel_table():
            if input.factor_choice() == "None":
                current_dataframe = du.get_metabolic_info(data.sample_data, current_dataframe, taxonomic_data)
                return render.DataGrid(
                    du.get_metabolic_info(data.sample_data, current_dataframe, taxonomic_data), row_selection_mode="single"
                )
            else:
                try:
                    factor_choice2 = int(input.factor_choice2())
                except:
                    return render.DataGrid(
                        current_dataframe.loc[current_dataframe[input.factor_choice()] == input.factor_choice2()],
                        row_selection_mode="single",
                    )
                else:
                    return render.DataGrid(
                        current_dataframe.loc[current_dataframe[input.factor_choice()] == factor_choice2], row_selection_mode="single"
                    )

        @output
        @render.data_frame
        def summary_table(with_filter: bool = False):
            if not with_filter:
                current_dataframe = du.get_metabolic_info(data.sample_data, metadata, taxonomic_data)
                return render.DataGrid(du.get_metabolic_info(data.sample_data, metadata, taxonomic_data), row_selection_mode="single")
            else:
                filter = input.test_sample2()
                dataframe_with_filter = current_dataframe.loc[
                    current_dataframe[input.current_data_choice1()] == input.current_data_choice2()
                ]
                return render.DataGrid(dataframe_with_filter, row_selection_mode="single")

        @output
        @render_widget
        def summary_taxonomy():
            # sample = current_dataframe.iloc[input.summary_table_selected_rows()[0],0]
            # sample_taxonomy = sm.taxonomic_processing(list_of_bin[sample],taxonomic_data)
            # fig = px.bar(sample_taxonomy["Genus"],x=sample_taxonomy["Genus"].index,y="Genus",title="Genus repartition in sample")
            # return fig
            fig = px.box(current_dataframe, x="Name", y="Numbers of models", color="Antibiotics")
            return fig

        @output
        @render.data_frame
        def bin_summary():
            if len(input.summary_table_selected_rows()) != 0:
                return render.DataGrid(
                    taxonomic_data.loc[
                        taxonomic_data["mgs"].isin(list_of_bin[current_dataframe.iloc[input.summary_table_selected_rows()[0], 0]])
                    ]
                )
            else:
                return current_dataframe

        @output
        @render.data_frame
        def bin_group():
            df = taxonomic_data.loc[
                taxonomic_data["mgs"].isin(list_of_bin[current_dataframe.iloc[input.summary_table_selected_rows()[0], 0]])
            ]
            df = df[["mgs", "Genus"]]
            df = df.groupby(["Genus"]).count()
            df = df.reset_index()
            df.columns = ["Genus", "Count"]
            return df

        @output
        @render.text
        def test_result():
            df = current_dataframe
            group1 = df[df["Antibiotics"] == "YES"]
            group2 = df[df["Antibiotics"] == "NO"]
            res = stats.wilcoxon(group1["Numbers of models"], group2["Numbers of models"])
            inshape_res = ("Résultat du test de wilcoxon :", res)
            return inshape_res

    app = App(app_ui, server)
    run_app(app=app, launch_browser=False, reload_dirs="./")
