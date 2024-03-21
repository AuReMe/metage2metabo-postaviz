import plotly.express as px
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


def custom_ontology():
    df_ontology = du.open_tsv("/home/lbrindel/m2m-postaviz/tests/compounds_26_5_level1.tsv")
    df_ontology_lvl2 = du.open_tsv("/home/lbrindel/Downloads/compounds_26_5_level2.tsv")
    du.add_row(df_ontology, ["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology, ["PSICOSE", "PSICOSE", "Others"])
    du.add_row(df_ontology_lvl2, ["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology_lvl2, ["PSICOSE", "PSICOSE", "Others"])
    df_ontology = du.deal_with_duplicated_row(df_ontology, "compound_name")
    df_ontology_lvl2 = du.deal_with_duplicated_row(df_ontology_lvl2, "compound_name")
    return df_ontology, df_ontology_lvl2


def run_shiny(data: DataStorage):
    ### Declare PROCESSING VARIABLE

    current_dataframe = data.main_data["metadata"]
    metadata = data.main_data["metadata"]
    taxonomic_data = data.taxonomic_data
    list_of_bin = data.get_bin_list()
    main_dataframe = data.main_data["main_dataframe"]
    factor_list = data.get_factor_list()
    print(factor_list)

    ### ALL GLOBAL OBJECT, TO BE REMOVED AT SOME POINT ###
    # main_df = du.merge_metadata_with_df(global_data["main_dataframe"],global_data["metadata"])
    # df_ontology, df_ontology_lvl2 = custom_ontology()

    ### Pathway ontology
    individual_pathway_ontology = sm.pathway_data_processing()
    ### ALL CARD OBJECT TO BE ARRANGED ###

    abundance_boxplot = ui.card(
        ui.row(
            ui.input_select("abplot_input1", "Label for X axis", factor_list),
            ui.input_selectize("abplot_input2", "Compound for Y axis", list(data.abundance_data.columns.values[data.get_factor_len():]))
        ),
        output_widget("Abundance_boxplot")
    )

    taxonomy_boxplot = ui.card(
        ui.row(
            ui.input_select("Taxonomic_rank_input", "Choose a taxonomic rank", list(taxonomic_data.columns), selected="Genus"),
            ui.input_selectize("Taxonomic_choice_input", "Multiple choice possible", [], multiple=True),
        ),
        output_widget("Taxonomic_boxplot"),
    )
    main_table = ui.card(ui.output_data_frame("main_table"))
    summary_table = ui.card(ui.output_data_frame("summary_table"))

    iscope_tab_card1 = ui.card(
        ui.row(
            ui.input_select("iscope_sample", label="Sample", choices=[key for key in list_of_bin.keys()]),
            ui.input_select("iscope_bin", label="Bin", choices=[key for key in individual_pathway_ontology[0].keys()]),
        )
    )

    iscope_taxonomic_card = ui.card(
        ui.layout_column_wrap(
            output_widget("taxonomic_domain"),
            output_widget("taxonomic_kingdom"),
            output_widget("taxonomic_phylum"),
            output_widget("taxonomic_order"),
            output_widget("taxonomic_family"),
            output_widget("taxonomic_genus"),
            width=1 / 2,
        )
    )

    iscope_pathway_card = ui.card(
        ui.layout_column_wrap(
            output_widget("pathwaylvl1_widget"),
            output_widget("pathwaylvl2_widget"),
            width=1 / 2,
        )
    )

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
                main_table
            ),
            ui.nav("Abundance",
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
                "Exploration",
                iscope_tab_card1,
                iscope_taxonomic_card,
                iscope_pathway_card,
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
        @render_widget
        def taxonomy_overview():
            df = du.taxonomic_overview(list_of_bin, taxonomic_data, metadata)
            plot = px.box(df, x="Antibiotics", y="Count", color="Days")
            return plot


        @output
        @render_widget
        def Abundance_boxplot():
            try:
                return px.box(data.abundance_data, x=input.abplot_input1(),y=input.abplot_input2())
            except:
                print("No valid input.")

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

        @render.text
        @reactive.event(input.run_test)
        def stat_test_result():
            if input.stat_test_choice() == "Student":
                return du.student_test(
                    current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group1_select()][
                        input.wilcoxon_tested_value_select()
                    ],
                    current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group2_select()][
                        input.wilcoxon_tested_value_select()
                    ],
                )
            if input.stat_test_choice() == "Wilcoxon":
                return du.wilcoxon_test(
                    current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group1_select()][
                        input.wilcoxon_tested_value_select()
                    ],
                    current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group2_select()][
                        input.wilcoxon_tested_value_select()
                    ],
                )

        @reactive.Effect()
        @reactive.event(input.wilcoxon_variable_select)
        def process_wilcoxon():
            current_factor_choice = del_list_duplicate(current_dataframe[input.wilcoxon_variable_select()])
            ui.update_select("wilcoxon_group1_select", choices=current_factor_choice)
            ui.update_select("wilcoxon_group2_select", choices=current_factor_choice)
            ui.update_select("wilcoxon_tested_value_select", choices=list(current_dataframe.columns))
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
        def main_table():
            return render.DataGrid(main_dataframe)


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

        @output
        @render.data_frame
        def iscope_table():
            # data = du.open_tsv("~/Downloads/mapping_mgs_genus.txt")
            df = data.sample_data[input.sample()]["iscope"]
            df = du.merge_df(df, taxonomic_data)
            return render.DataGrid(df, row_selection_mode="single")

        @output
        @render_widget
        def taxonomic_domain():
            current_sample = sm.taxonomic_processing(list_of_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Domain"], x=current_sample["Domain"].index, y="Domain", title="Domain repartition in sample")
            return fig

        @output
        @render_widget
        def taxonomic_kingdom():
            current_sample = sm.taxonomic_processing(list_of_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Kingdom"], x=current_sample["Kingdom"].index, y="Kingdom", title="Kingdom repartition in sample")
            return fig

        @output
        @render_widget
        def taxonomic_phylum():
            current_sample = sm.taxonomic_processing(list_of_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Phylum"], x=current_sample["Phylum"].index, y="Phylum", title="Phylum repartition in sample")
            return fig

        @output
        @render_widget
        def taxonomic_order():
            current_sample = sm.taxonomic_processing(list_of_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Order"], x=current_sample["Order"].index, y="Order", title="Order repartition in sample")
            return fig

        @output
        @render_widget
        def taxonomic_family():
            current_sample = sm.taxonomic_processing(list_of_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Family"], x=current_sample["Family"].index, y="Family", title="Family repartition in sample")
            return fig

        @output
        @render_widget
        def taxonomic_genus():
            current_sample = sm.taxonomic_processing(list_of_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Genus"], x=current_sample["Genus"].index, y="Genus", title="Genus repartition in sample")
            return fig

        @output
        @render_widget
        def pathwaylvl1_widget():
            fig = px.bar(
                individual_pathway_ontology[0][input.iscope_bin()],
                x="category",
                y="count",
                title="Ontology lvl1 for " + str(input.iscope_bin()),
            )
            return fig

        @output
        @render_widget
        def pathwaylvl2_widget():
            fig = px.bar(
                individual_pathway_ontology[1][input.iscope_bin()],
                x="category",
                y="count",
                title="Ontology lvl2 for " + str(input.iscope_bin()),
            )
            return fig

    app = App(app_ui, server)
    run_app(app=app, launch_browser=True)
