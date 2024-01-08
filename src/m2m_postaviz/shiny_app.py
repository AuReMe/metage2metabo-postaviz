import plotly.express as px
from shiny import App
from shiny import render
from shiny import run_app
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.data_utils as du
import m2m_postaviz.rpy2_utils as ru


def run_shiny(all_data, data):
    main_df = du.merge_metadata_with_df(all_data["main_dataframe"],data)
    converter = ru.Rconverter
    pcoa_results = converter.pcoa(converter, main_df, du.get_column_size(data))
    df_ontology = du.open_tsv("/home/lbrindel/m2m-postaviz/tests/compounds_26_5_level1.tsv")
    df_ontology_lvl2 = du.open_tsv("/home/lbrindel/Downloads/compounds_26_5_level2.tsv")
    du.add_row(df_ontology,["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology,["PSICOSE", "PSICOSE", "Others"])
    du.add_row(df_ontology_lvl2,["ALTROSE", "ALTROSE", "Others"])
    du.add_row(df_ontology_lvl2,["PSICOSE", "PSICOSE", "Others"])
    df_ontology = du.deal_with_duplicated_row(df_ontology, "compound_name")
    df_ontology_lvl2 = du.deal_with_duplicated_row(df_ontology_lvl2, "compound_name")

    meta_card = ui.card(ui.input_selectize(
                id="search_metadata", label="meta_searchbar", choices=list(data.columns.values), width="400px", multiple=True
                ),
                ui.output_table("metadata_table"),
                )
    
    data_card = ui.card(ui.input_selectize(
                id="search_data", label="data_searchbar", choices=list(all_data["main_dataframe"].columns.values[du.get_column_size(data) + 1 :]), width="500px", multiple=True,
                ),
                ui.output_table("data_table"),
                )
    
    output_pcoa = ui.card(ui.row(ui.input_select("x", label="Symbol", choices=list(data.columns.values[1:])),
                  ui.input_select("color", label="Color", choices=list(data.columns.values[1:])),
                  ),
                  output_widget("my_widget", width=600),
                  )

    sample_choice =ui.card(ui.row(ui.input_select("sample", label="Sample", choices=list(all_data["iscope"])),ui.input_select("col", label="Column", choices=["antibio", "days"])
                    ),
                    ui.output_table("iscope_table"),
                    )

    ontology_card = ui.card(ui.card_header("Ontology level 1"),ui.row(ui.input_select("scope_ontology", label="Sample", choices=list(all_data["main_dataframe"]["Name"]))
                    ),
                    output_widget("ontology_barplot"),
                    )
    
    ontology_card2 = ui.card(ui.card_header("Ontology level 2"),ui.row(ui.input_select("scope_ontology2", label="Sample", choices=list(all_data["main_dataframe"]["Name"]))
                    ),
                    output_widget("ontology_barplot2"),
                    )


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
            current_df = all_data["main_dataframe"].set_index('Name')
            df = df_ontology.loc[df_ontology.index.isin(current_df.columns.values)]
            count_matrix = current_df@df
            percent_matrix = du.get_percent_value(count_matrix,True)
            results = du.get_threshold_value(percent_matrix, 1)
            fig = px.bar(results[input.scope_ontology()],x=results[input.scope_ontology()].index, y=results[input.scope_ontology()])
            return fig

        @output
        @render_widget
        def ontology_barplot2():
            current_df = all_data["main_dataframe"].set_index('Name')
            df = df_ontology_lvl2.loc[df_ontology_lvl2.index.isin(current_df.columns.values)]
            count_matrix = current_df@df
            percent_matrix = du.get_percent_value(count_matrix,True)
            results = du.get_threshold_value(percent_matrix, 1)
            fig = px.bar(results[input.scope_ontology2()],x=results[input.scope_ontology2()].index, y=results[input.scope_ontology2()])
            return fig

        @output
        @render.table
        def metadata_table():
            column = du.get_columns_index(main_df, input.search_metadata())
            df = main_df.iloc[:, column]
            return df

        @output
        @render.table
        def data_table():
            column = du.get_columns_index(all_data["main_dataframe"], input.search_data())
            df = all_data["main_dataframe"].iloc[:, column]
            return df

        @output
        @render.table
        def iscope_table():
            data = du.open_tsv("~/Downloads/mapping_mgs_genus.txt")
            df = all_data["iscope"][input.sample()]
            df = du.merge_df(df,data)
            df = df.groupby("genus")
            return df
        
    app = App(app_ui, server)
    run_app(app=app, launch_browser=True)
