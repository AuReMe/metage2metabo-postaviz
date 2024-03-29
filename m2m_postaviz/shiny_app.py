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
import m2m_postaviz.shiny_module as sm
from m2m_postaviz.data_struct import DataStorage


def del_list_duplicate(mylist: list):
    return list(dict.fromkeys(mylist))

def run_shiny(data: DataStorage):
    ### Declare PROCESSING VARIABLE

    current_dataframe = data.main_data["metadata"]
    metadata = data.main_data["metadata"]
    taxonomic_data = data.taxonomic_data
    list_of_bin = data.get_bin_list()
    main_dataframe = data.main_data["main_dataframe"]
    factor_list = data.list_of_factor
    factor_list.insert(0, 'None')
    print(factor_list)


    ### ALL CARD OBJECT TO BE ARRANGED ###

    abundance_input = ui.row(
            ui.input_select("box_inputx1", "Label for X axis", factor_list),
            ui.input_select("box_inputx2", "Label for 2nd X axis", factor_list),
            ui.input_selectize("box_inputy1", "Compound for Y axis", data.abundance_matrix.columns.tolist(),multiple=True)
        )
    abundance_boxplot = ui.card(
        output_widget("Abundance_boxplot")
    )

    taxonomy_boxplot = ui.card(
        ui.row(
            ui.input_select("Taxonomic_rank_input", "Choose a taxonomic rank", list(taxonomic_data.columns), selected="Genus"),
            ui.input_selectize("Taxonomic_choice_input", "Multiple choice possible", [], multiple=True),
        ),
        output_widget("Taxonomic_boxplot"),
    )
    main_table = ui.card(ui.output_data_frame("dev_table"),output_widget("main_table"))
    # dev_table = ui.card(ui.output_text("dev_text"))
    summary_table = ui.card(ui.output_data_frame("summary_table"))

    main_panel_dataframe = ui.card(
        ui.row(
            ui.input_select(id="factor_choice", label="Filter", choices=factor_list, selected=factor_list[0]),
            ui.input_select(id="factor_choice2", label="Value", choices=[]),
        ),
        ui.output_data_frame("main_panel_table"),
    )

    ### APPLICATION TREE ###

    app_ui = ui.page_fluid(
        ui.navset_tab(
            ui.nav("Dev mod",
                main_table,
                # dev_table
            ),
            ui.nav("Abundance",
                   abundance_input,
                   abundance_boxplot
                   ),
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

        @output
        @render_widget()
        def taxonomy_overview():
            df = du.taxonomic_overview(list_of_bin, taxonomic_data, metadata)
            ### IF NO INPUT -- DEFAULT

            plot = px.box(df, y="nb_taxon")

            # plot = px.box(df, x="Antibiotics", y="Count", color="Days")
            return plot


        @output
        @render_widget
        def Abundance_boxplot():
            df = data.produce_long_abundance_dataframe()
            df.drop(df.loc[df["Quantity"] == 0].index, inplace=True)
            print(len(input.box_inputy1()))
            y1 = input.box_inputy1()
            if len(y1) == 0:
                y1 = 'Quantity'
            if not input.box_inputx1() == 'None':
                if input.box_inputx2() == 'None':
                    return px.box(df,x=input.box_inputx1(), y='Quantity')
                else:
                    conditionx2 = df[input.box_inputx2()]
                    conditionx2 = conditionx2.unique()

                    fig = go.Figure()
                    
                    for condition in conditionx2:

                        fig.add_trace(go.Box(
                            x=df.loc[df[input.box_inputx2()] == condition][input.box_inputx1()] ,
                            y=df.loc[(df[input.box_inputx2()] == condition) & (df["Quantity"] != 0)]['Quantity'],
                            name=str(condition),
                            # marker_color='green'
                        ))
                    # fig.add_trace(go.Box(
                    #     x=df[input.box_inputx1()],
                    #     y=df.loc[(df["Group"] == 'Treatment') & (df["Quantity"] != 0)]['Quantity'],
                    #     name='Treatment',
                    #     marker_color='blue'
                    # ))
                    fig.update_layout(boxmode='group')
                    return fig
            ##### MAKE PLOT
            return px.box(df, y='Quantity')

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


        @render.data_frame
        def dev_table():
            return render.DataGrid(data.melted_abundance_dataframe)


        @render_widget
        def main_table():
            df = data.produce_long_abundance_dataframe()
            fig = px.box(df, x='Days', y='Quantity', color='Group')
            return fig


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
    run_app(app=app, launch_browser=True,reload_dirs="./")
