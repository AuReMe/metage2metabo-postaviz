import plotly.express as px
from shiny import App
from shiny import render
from shiny import run_app
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.data_utils as du
import m2m_postaviz.rpy2_utils as ru


def run_shiny(main_df, data, all_iscopes):
    converter = ru.Rconverter
    pcoa_results = converter.pcoa(converter, main_df, du.get_column_size(data))


    meta_card = ui.card(ui.input_selectize(
                id="search_metadata", label="meta_searchbar", choices=list(data.columns.values), width="400px", multiple=True
                ),
                ui.output_table("metadata_table"),
                )
    
    data_card = ui.card(ui.input_selectize(
                id="search_data", label="data_searchbar", choices=list(main_df.columns.values[du.get_column_size(data) + 1 :]), width="500px", multiple=True,
                ),
                ui.output_table("data_table"),
                )
    
    output_pcoa = ui.card(ui.row(ui.input_select("x", label="Symbol", choices=list(data.columns.values[1:])),
            ui.input_select("color", label="Color", choices=list(data.columns.values[1:])),
            ),
            output_widget("my_widget", width=600),
    )

    sample_choice =ui.card(ui.row(ui.input_select("sample", label="Sample", choices=list(all_iscopes)),ui.input_select("col", label="Column", choices=["antibio", "days"])
                    ),
                    ui.output_table("iscope_table"),
                    )



    app_ui = ui.page_fluid(
        ui.navset_tab(
            ui.nav("cscope",
                ui.layout_column_wrap(
                meta_card,
                data_card,
                output_pcoa,
                width= 1/1
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
        @render.table
        def metadata_table():
            column = du.get_columns_index(main_df, input.search_metadata())
            df = main_df.iloc[:, column]
            return df

        @output
        @render.table
        def data_table():
            column = du.get_columns_index(main_df, input.search_data())
            df = main_df.iloc[:, column]
            return df

        @output
        @render.table
        def iscope_table():
            df = all_iscopes[input.sample()]
            df = df.iloc[:,:5]
            return df
        
    app = App(app_ui, server)
    run_app(app=app, launch_browser=True)
