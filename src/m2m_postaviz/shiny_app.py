import plotly.express as px
from shiny import App, render, run_app, ui, reactive

from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.data_utils as du
import m2m_postaviz.rpy2_utils as ru
import m2m_postaviz.shiny_module as sm

def run_shiny(global_data, sample_data):

    ### ALL GLOBAL OBJECT, TO BE REMOVED AT SOME POINT ###
    main_df = du.merge_metadata_with_df(global_data["main_dataframe"],global_data["metadata"])
    converter = ru.Rconverter
    pcoa_results = converter.pcoa(converter, main_df, du.get_column_size(global_data["metadata"]))
    df_ontology = du.open_tsv("/home/lbrindel/m2m-postaviz/tests/compounds_26_5_level1.tsv")
    df_ontology_lvl2 = du.open_tsv("/home/lbrindel/Downloads/compounds_26_5_level2.tsv")
    du.add_row(df_ontology,["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology,["PSICOSE", "PSICOSE", "Others"])
    du.add_row(df_ontology_lvl2,["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology_lvl2,["PSICOSE", "PSICOSE", "Others"])
    df_ontology = du.deal_with_duplicated_row(df_ontology, "compound_name")
    df_ontology_lvl2 = du.deal_with_duplicated_row(df_ontology_lvl2, "compound_name")

    ### Variable from shiny module ###

    ### Taxonomy
    taxonomic_data = du.open_tsv("/home/lbrindel/m2m-postaviz/tests/taxonomic_database.tsv")
    list_of_individual_bin = sm.read_rev_iscope("rev_iscope.tsv", "/home/lbrindel/cscope_metadata/", "Unnamed: 0")
    ### Pathway ontology 
    individual_pathway_ontology = sm.pathway_data_processing()
    ### ALL CARD OBJECT TO BE ARRANGED ###
    meta_card = ui.card(ui.input_selectize(
                id="search_metadata", label="meta_searchbar", choices=list(global_data["metadata"].columns.values), width="400px", multiple=True
                ),
                ui.output_data_frame("metadata_table"),
                )
    
    data_card = ui.card(ui.input_selectize(
                id="search_data", label="data_searchbar", choices=list(global_data["main_dataframe"].columns.values[du.get_column_size(global_data["metadata"]) + 1 :]), width="500px", multiple=True,
                ),
                ui.output_data_frame("data_table"),
                )
    
    output_pcoa = ui.card(ui.row(ui.input_select("x", label="Symbol", choices=list(global_data["metadata"].columns.values[1:])),
                  ui.input_select("color", label="Color", choices=list(global_data["metadata"].columns.values[1:])),
                  ),
                  output_widget("my_widget", width=600),
                  )

    sample_choice =ui.card(ui.row(ui.input_select("sample", label="Sample", choices=list(key for key in sample_data.keys())),ui.input_select("col", label="Column", choices=["antibio", "days"])
                    ),
                    ui.layout_column_wrap(ui.output_data_frame("iscope_table"),ui.output_text(id="row_been_selected")),
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

    ### APPLICATION TREE ###

    app_ui = ui.page_fluid(
        ui.navset_tab(
            ui.nav("cscope",
                ui.layout_column_wrap(
                meta_card,
                data_card,
                ontology_card2,
                ontology_card,
                output_pcoa,
                width= 1/2
                )
                   ),
            ui.nav("iscope",
                   sample_choice
                   ),
            ui.nav("Exploration",
                   iscope_tab_card1,
                   iscope_taxonomic_card,
                    iscope_pathway_card
            ),
            ui.nav("Prototype",
                   summary_table,
                   
            )
        )

    )
    def server(input, output, session):
        @output
        @render_widget
        def my_widget():
            fig = px.scatter(
                pcoa_results,
                x="Dim1",
                y="Dim2",
                color=input.color(),
                symbol=input.x(),
            )
            return fig
        
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


        @output
        @render.data_frame
        def summary_table():
            return du.get_metabolic_info(sample_data)


        @output
        @render.data_frame
        def metadata_table():
            column = du.get_columns_index(main_df, input.search_metadata())
            df = main_df.iloc[:, column]
            return df


        @output
        @render.data_frame
        def data_table():
            column = du.get_columns_index(global_data["main_dataframe"], input.search_data())
            df = global_data["main_dataframe"].iloc[:, column]
            return df


        @output
        @render.data_frame
        def iscope_table():
            # data = du.open_tsv("~/Downloads/mapping_mgs_genus.txt")
            df = sample_data[input.sample()]["iscope"]
            df = du.merge_df(df,taxonomic_data)
            return render.DataGrid(df, row_selection_mode="single")
        

        @output
        @render.text
        def row_been_selected():
            selected = input.iscope_table_selected_rows()
            thing = ""
            for f in selected:
                thing += str(f)
            return thing


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
