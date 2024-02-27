import plotly.express as px
from shiny import App, render, run_app, ui, reactive
import scipy.stats as stats
from m2m_postaviz.time_decorator import timeit

from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.data_utils as du
import m2m_postaviz.rpy2_utils as ru
import m2m_postaviz.shiny_module as sm


def clean_list_duplicate(mylist: list):
    return list(dict.fromkeys(mylist))


def custom_ontology():
    df_ontology = du.open_tsv("/home/lbrindel/m2m-postaviz/tests/compounds_26_5_level1.tsv")
    df_ontology_lvl2 = du.open_tsv("/home/lbrindel/Downloads/compounds_26_5_level2.tsv")
    du.add_row(df_ontology,["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology,["PSICOSE", "PSICOSE", "Others"])
    du.add_row(df_ontology_lvl2,["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology_lvl2,["PSICOSE", "PSICOSE", "Others"])
    df_ontology = du.deal_with_duplicated_row(df_ontology, "compound_name")
    df_ontology_lvl2 = du.deal_with_duplicated_row(df_ontology_lvl2, "compound_name")
    return df_ontology, df_ontology_lvl2


def run_shiny(global_data, sample_data):

    ### Declare PROCESSING VARIABLE

    global_data["current_dataframe"] = global_data["metadata"]
    current_dataframe = global_data["current_dataframe"]
    factor_list = list(global_data["current_dataframe"].columns.values)
    factor_list.insert(0,"None")
    print(factor_list, type(factor_list))

    ### ALL GLOBAL OBJECT, TO BE REMOVED AT SOME POINT ###
    main_df = du.merge_metadata_with_df(global_data["main_dataframe"],global_data["metadata"])
    converter = ru.Rconverter
    pcoa_results = converter.pcoa(converter, main_df, du.get_column_size(global_data["metadata"]))
    df_ontology, df_ontology_lvl2 = custom_ontology()

    ### Variable from shiny module ###

    ### Taxonomy
    taxonomic_data = du.open_tsv("/home/lbrindel/m2m-postaviz/tests/taxonomic_database.tsv")
    list_of_individual_bin = sm.read_rev_iscope("rev_iscope.tsv", "/home/lbrindel/cscope_metadata/", "Unnamed: 0")
    ### Pathway ontology 
    individual_pathway_ontology = sm.pathway_data_processing()
    ### ALL CARD OBJECT TO BE ARRANGED ###

    taxonomy_boxplot = ui.card(ui.row(ui.input_select("Taxonomic_rank_input", "Choose a taxonomic rank", list(taxonomic_data.columns), selected="Genus"),
                                      ui.input_selectize("Taxonomic_choice_input", "Multiple choice possible", [], multiple=True)
                                      ),
                                output_widget("Taxonomic_boxplot")
                               )

    summary_table = ui.card(ui.output_data_frame("summary_table"))

    ontology_card = ui.card(ui.card_header("Ontology level 1"),ui.row(ui.input_select("scope_ontology", label="Sample", choices=list(global_data["main_dataframe"]["Name"]))
                    ),
                    output_widget("ontology_barplot"),
                    )
    
    ontology_card2 = ui.card(ui.card_header("Ontology level 2"),ui.row(ui.input_select("scope_ontology2", label="Sample", choices=list(global_data["main_dataframe"]["Name"]))
                    ),
                    output_widget("ontology_barplot2"),
                    )

    iscope_tab_card1 = ui.card(
                            ui.row(
                                    ui.input_select("iscope_sample", label="Sample", choices=[key for key in list_of_individual_bin.keys()]),
                                    ui.input_select("iscope_bin", label="Bin", choices=[key for key in individual_pathway_ontology[0].keys()])
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
                                        width=1/2,
                                    )
                                )
                                
    iscope_pathway_card = ui.card(
                                    ui.layout_column_wrap(
                                        output_widget("pathwaylvl1_widget"),
                                        output_widget("pathwaylvl2_widget"),
                                        width=1/2,
                                    )
                                )

    main_panel_dataframe = ui.card(
                            ui.row(
                                ui.input_select(id="factor_choice",label="Filter",choices=factor_list,selected=factor_list[0]),
                                ui.input_select(id="factor_choice2",label="Value",choices=[]),
                            ),
        ui.output_data_frame("main_panel_table")
        )

    ### APPLICATION TREE ###

    app_ui = ui.page_fluid(
        ui.navset_tab(
            ui.nav("Taxonomy",
                   output_widget("taxonomy_overview"),
                    taxonomy_boxplot
                   ),
            ui.nav("Main",
                   ui.layout_sidebar(
                       ui.sidebar(ui.input_select(id="stat_test_choice",label="Statistical test available.",choices=["Choose test","Student","Wilcoxon","ANOVA","PCOA"]),
                                  ui.panel_conditional(
                                    "input.stat_test_choice === 'Wilcoxon' || input.stat_test_choice === 'Student'",
                                        ui.input_select(id="wilcoxon_variable_select",label="Select a variable",choices=list(current_dataframe.columns)),
                                        ui.input_select(id="wilcoxon_group1_select",label="Group1 choice",choices=[]),
                                        ui.input_select(id="wilcoxon_group2_select",label="Group2 choice",choices=[]),
                                        ui.input_select(id="wilcoxon_tested_value_select",label="Tested value choice",choices=[]),
                                  ),
                                  ui.input_action_button("run_test", "run_test"),
                                  
                                  ),
                       main_panel_dataframe,
                       ui.output_text(id="stat_test_result")
                   )
                   ),
            ui.nav("Exploration",
                   iscope_tab_card1,
                   iscope_taxonomic_card,
                    iscope_pathway_card,
                    ontology_card,
                    ontology_card2
            ),
            ui.nav("Prototype",
                   ui.card( 
                            ui.output_data_frame("summary_table")
                        ),
                        
                        ui.accordion(
                        ui.accordion_panel("Taxonomy panel",
                        
                            ui.layout_column_wrap(
                                      ui.output_data_frame("bin_summary"),
                                      ui.output_data_frame("bin_group"),
                                      width=1/2)
                        ),open=False),
                        ui.card(output_widget("summary_taxonomy")),

                        ui.card(
                            ui.input_action_button("launch_test", "Run"),
                            ui.output_text(id="test_result")
                              ),
                            )
        )
    )
    def server(input, output, session):


        @reactive.Effect()
        @reactive.event(input.Taxonomic_rank_input)
        def updapte_taxonomy_choice_list():
            updated_list = clean_list_duplicate(taxonomic_data[input.Taxonomic_rank_input()])
            ui.update_select("Taxonomic_choice_input", choices=updated_list)
            return


        @output
        @render_widget
        def taxonomy_overview():
            df = du.taxonomic_overview(list_of_individual_bin, taxonomic_data, global_data["metadata"])
            plot = px.box(df, x="Antibiotics",y='Count',color='Days')
            return plot


        @output
        @render_widget
        def Taxonomic_boxplot():
            df = du.taxonomic_dataframe_from_input(input.Taxonomic_rank_input(), list_of_individual_bin, input.Taxonomic_choice_input(), taxonomic_data, global_data["metadata"])
            plot = px.box(df, x="Antibiotics",y='Count',color='Antibiotics')
            return plot


        @render.text
        @reactive.event(input.run_test)
        def stat_test_result():
            if input.stat_test_choice() == "Student":
                return du.student_test(
                current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group1_select()][input.wilcoxon_tested_value_select()],
                current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group2_select()][input.wilcoxon_tested_value_select()]
                )
            if input.stat_test_choice() == "Wilcoxon":
                return du.wilcoxon_test(
                current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group1_select()][input.wilcoxon_tested_value_select()],
                current_dataframe.loc[current_dataframe[input.wilcoxon_variable_select()] == input.wilcoxon_group2_select()][input.wilcoxon_tested_value_select()]
                )


        @reactive.Effect()
        @reactive.event(input.wilcoxon_variable_select)
        def process_wilcoxon():
            current_factor_choice = clean_list_duplicate(current_dataframe[input.wilcoxon_variable_select()])
            ui.update_select("wilcoxon_group1_select", choices=current_factor_choice)
            ui.update_select("wilcoxon_group2_select", choices=current_factor_choice)
            ui.update_select("wilcoxon_tested_value_select", choices=list(current_dataframe.columns))
            return


        @reactive.Effect()
        @reactive.event(input.test_sample1)
        def process():
            current_factor_choice = clean_list_duplicate(current_dataframe[input.test_sample1()])
            ui.update_select("test_sample2", choices=current_factor_choice)
            return
        

        @reactive.Effect()
        @reactive.event(input.factor_choice)
        def process():
            if input.factor_choice() == "None":
                return
            current_factor_choice = clean_list_duplicate(current_dataframe[input.factor_choice()])
            # current_factor_choice.insert(0,"None")
            ui.update_select("factor_choice2", choices=current_factor_choice)
            return


        @output
        @render_widget
        def ontology_barplot():
            current_df = global_data["main_dataframe"].set_index('Name')
            df = df_ontology.loc[df_ontology.index.isin(current_df.columns.values)]
            count_matrix = current_df@df
            percent_matrix = du.get_percent_value(count_matrix,True)
            results = du.get_threshold_value(percent_matrix, 1)
            fig = px.bar(results[input.scope_ontology()],x=results[input.scope_ontology()].index, y=results[input.scope_ontology()])
            return fig


        @output
        @render_widget
        def ontology_barplot2():
            current_df = global_data["main_dataframe"].set_index('Name')
            df = df_ontology_lvl2.loc[df_ontology_lvl2.index.isin(current_df.columns.values)]
            count_matrix = current_df@df
            percent_matrix = du.get_percent_value(count_matrix,True)
            results = du.get_threshold_value(percent_matrix, 1)
            fig = px.bar(results[input.scope_ontology2()],x=results[input.scope_ontology2()].index, y=results[input.scope_ontology2()])
            return fig


        @render.data_frame
        def main_panel_table():
            if input.factor_choice() == factor_list[0]:
                global_data["current_dataframe"] = du.get_metabolic_info(sample_data, global_data["current_dataframe"], taxonomic_data)
                return render.DataGrid(du.get_metabolic_info(sample_data, global_data["current_dataframe"], taxonomic_data), row_selection_mode="single")
            else:
                try:
                    factor_choice2 = int(input.factor_choice2())
                except:
                    return render.DataGrid(global_data["current_dataframe"].loc[global_data["current_dataframe"][input.factor_choice()] == input.factor_choice2()], row_selection_mode="single")
                else:
                    return render.DataGrid(global_data["current_dataframe"].loc[global_data["current_dataframe"][input.factor_choice()] == factor_choice2], row_selection_mode="single")


        @output
        @render.data_frame
        def summary_table(with_filter: bool = False):
            if not with_filter:
                global_data["current_dataframe"] = du.get_metabolic_info(sample_data, global_data["metadata"], taxonomic_data)
                return render.DataGrid(du.get_metabolic_info(sample_data, global_data["metadata"], taxonomic_data), row_selection_mode="single")
            else:
                filter = input.test_sample2()
                dataframe_with_filter = global_data["current_dataframe"].loc[global_data["current_dataframe"][input.current_data_choice1()] == input.current_data_choice2()]
                return render.DataGrid(dataframe_with_filter,row_selection_mode="single")


        @output
        @render_widget
        def summary_taxonomy():
            # sample = current_dataframe.iloc[input.summary_table_selected_rows()[0],0]
            # sample_taxonomy = sm.taxonomic_processing(list_of_individual_bin[sample],taxonomic_data)
            # fig = px.bar(sample_taxonomy["Genus"],x=sample_taxonomy["Genus"].index,y="Genus",title="Genus repartition in sample")
            # return fig
            fig = px.box(global_data["current_dataframe"], x="Name", y="Numbers of models", color="Antibiotics")
            return fig
        

        @output
        @render.data_frame
        def bin_summary():
            if len(input.summary_table_selected_rows()) != 0:
                return render.DataGrid(taxonomic_data.loc[taxonomic_data["mgs"].isin(list_of_individual_bin[global_data["current_dataframe"].iloc[input.summary_table_selected_rows()[0],0]])])
            else: 
                return global_data["current_dataframe"]


        @output
        @render.data_frame
        def bin_group():
            df = taxonomic_data.loc[taxonomic_data["mgs"].isin(list_of_individual_bin[global_data["current_dataframe"].iloc[input.summary_table_selected_rows()[0],0]])]
            df = df[["mgs",'Genus']]
            df = df.groupby(["Genus"]).count()
            df = df.reset_index()
            df.columns = ["Genus", "Count"]
            return df


        @output
        @render.text
        def test_result():
            df = global_data["current_dataframe"]
            group1 = df[df["Antibiotics"]=="YES"]
            group2 = df[df["Antibiotics"]=="NO"]
            res = stats.wilcoxon(group1["Numbers of models"],group2["Numbers of models"])
            inshape_res = ("Résultat du test de wilcoxon :", res)
            return inshape_res


        @output
        @render.data_frame
        def iscope_table():
            # data = du.open_tsv("~/Downloads/mapping_mgs_genus.txt")
            df = sample_data[input.sample()]["iscope"]
            df = du.merge_df(df,taxonomic_data)
            return render.DataGrid(df, row_selection_mode="single")
        

        @output
        @render_widget
        def taxonomic_domain():
            current_sample = sm.taxonomic_processing(list_of_individual_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Domain"],x=current_sample["Domain"].index,y="Domain",title="Domain repartition in sample")
            return fig
        

        @output
        @render_widget
        def taxonomic_kingdom():
            current_sample = sm.taxonomic_processing(list_of_individual_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Kingdom"],x=current_sample["Kingdom"].index,y="Kingdom",title="Kingdom repartition in sample")
            return fig
        

        @output
        @render_widget
        def taxonomic_phylum():
            current_sample = sm.taxonomic_processing(list_of_individual_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Phylum"],x=current_sample["Phylum"].index,y="Phylum",title="Phylum repartition in sample")
            return fig
        

        @output
        @render_widget
        def taxonomic_order():
            current_sample = sm.taxonomic_processing(list_of_individual_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Order"],x=current_sample["Order"].index,y="Order",title="Order repartition in sample")
            return fig
        

        @output
        @render_widget
        def taxonomic_family():
            current_sample = sm.taxonomic_processing(list_of_individual_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Family"],x=current_sample["Family"].index,y="Family",title="Family repartition in sample")
            return fig
        

        @output
        @render_widget
        def taxonomic_genus():
            current_sample = sm.taxonomic_processing(list_of_individual_bin[input.iscope_sample()], taxonomic_data)
            fig = px.bar(current_sample["Genus"],x=current_sample["Genus"].index,y="Genus",title="Genus repartition in sample")
            return fig
        

        @output
        @render_widget
        def pathwaylvl1_widget():
            fig = px.bar(individual_pathway_ontology[0][input.iscope_bin()],x="category",y="count",title="Ontology lvl1 for "+str(input.iscope_bin()))
            return fig
        

        @output
        @render_widget
        def pathwaylvl2_widget():
            fig = px.bar(individual_pathway_ontology[1][input.iscope_bin()],x="category",y="count",title="Ontology lvl2 for "+str(input.iscope_bin()))
            return fig
        

    app = App(app_ui, server)
    run_app(app=app, launch_browser=True)
