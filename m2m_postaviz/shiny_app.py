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
    ### Declare PROCESSING VARIABLE

    current_dataframe = data.main_data["metadata"]
    metadata = data.main_data["metadata"]
    taxonomic_data = data.taxonomic_data
    list_of_bin = data.get_bin_list()
    main_dataframe = data.main_data["main_dataframe"].copy()
    factor_list = data.list_of_factor
    factor_list.insert(0, "None")
    print(factor_list)
    metadata_label = data.get_metadata_label()

    ### ALL CARD OBJECT TO BE ARRANGED ###

    abundance_input = ui.row(
        ui.input_select("box_inputx1", "Label for X axis", factor_list),
        ui.input_select("box_inputx2", "Label for 2nd X axis", factor_list),
        ui.input_selectize("box_inputy1", "Compound for Y axis", data.abundance_matrix.columns.tolist(), multiple=True),
        ui.input_checkbox("ab_norm", "With normalised data")
        # ui.input_action_button("abundance_boxplot_button","Launch")
    )
    abundance_boxplot = ui.card(output_widget("Abundance_boxplot", height="100%", width="100%"))

    taxonomy_boxplot = ui.card(
        ui.row(
            ui.input_select("Taxonomic_rank_input", "Choose a taxonomic rank", list(taxonomic_data.columns), selected="Genus"),
            ui.input_selectize("Taxonomic_choice_input", "Multiple choice possible", [], multiple=True),
        ),
        output_widget("Taxonomic_boxplot"),
    )
    main_table = ui.card(ui.output_data_frame("dev_table"), ui.output_data_frame("main_table"))
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
            ui.nav("Abundance", ui.layout_sidebar(ui.sidebar(abundance_input), abundance_boxplot)),
            ui.nav("Taxonomy", output_widget("taxonomy_overview"), taxonomy_boxplot),
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
            ui.nav(
                "Prototype",
                ui.card(ui.output_data_frame("summary_table")),
                ui.accordion(
                    ui.accordion_panel(
                        "Taxonomy panel",
                        ui.layout_column_wrap(ui.output_data_frame("bin_summary"), ui.output_data_frame("bin_group"), width=1 / 2),
                    ),
                    open=False,
                ),
                ui.card(output_widget("summary_taxonomy")),
                ui.card(ui.input_action_button("launch_test", "Run"), ui.output_text(id="test_result")),
            ),
        )
    )

    def server(input, output, session):
        @reactive.Effect()
        @reactive.event(input.Taxonomic_rank_input)
        def update_taxonomy_choice_list():
            updated_list = del_list_duplicate(taxonomic_data[input.Taxonomic_rank_input()])
            ui.update_select("Taxonomic_choice_input", choices=updated_list)
            return

        @render_widget()
        def taxonomy_overview():
            df = du.taxonomic_overview(list_of_bin, taxonomic_data, metadata)
            ### IF NO INPUT -- DEFAULT

            plot = px.box(df, y="nb_taxon")

            # plot = px.box(df, x="Antibiotics", y="Count", color="Days")
            return plot

        @render_widget
        def Abundance_boxplot():
            with_normalised_data = input.ab_norm()
            if with_normalised_data:
                df = data.get_melted_norm_ab_dataframe()
            else:
                df = data.get_melted_ab_dataframe()
            # df.drop(df.loc[df["Quantity"] == 0].index, inplace=True)
            y1 = input.box_inputy1()
            x1 = input.box_inputx1()
            x2 = input.box_inputx2()
            if len(y1) == 0:
                y1 = df["Compound"].unique()
            if x1 == "None":
                return px.box(df, y=df.loc[df["Quantity"] != 0 & df["Compound"].isin(y1)]["Quantity"])
            else:
                if x2 == "None":
                    return px.box(
                        df,
                        x=df.loc[(df["Compound"].isin(y1)) & (df["Quantity"] != 0)][x1],
                        y=df.loc[(df["Compound"].isin(y1)) & (df["Quantity"] != 0)]["Quantity"],
                    )
                else:
                    conditionx2 = df[x2]
                    conditionx2 = conditionx2.unique()

                    fig = go.Figure()

                    for condition in conditionx2:
                        fig.add_trace(
                            go.Box(
                                x=df.loc[df[x2] == condition][x1],
                                y=df.loc[(df[x2] == condition) & (df["Quantity"] != 0) & (df["Compound"].isin(y1))]["Quantity"],
                                name=str(condition),
                                # marker_color='green'
                            )
                        )

                    fig.update_layout(boxmode="group", hovermode="y")
                    return fig

        @output
        @render_widget
        def Taxonomic_boxplot():
            try:
                df = du.taxonomic_dataframe_from_input(
                    input.Taxonomic_rank_input(), list_of_bin, input.Taxonomic_choice_input(), taxonomic_data, metadata
                )
                plot = px.box(df, x="Antibiotics", y="Count", color="Antibiotics")
                return plot
            except:
                return

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
            return render.DataGrid(data.abundance_dataframe)

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
