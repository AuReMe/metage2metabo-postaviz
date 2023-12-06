import plotly.express as px
from shiny import App
from shiny import render
from shiny import run_app
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.data_utils as du
import m2m_postaviz.rpy2_utils as ru


def run_shiny(main_df, data):
    converter = ru.Rconverter
    pcoa_results = converter.pcoa(converter, main_df, du.get_column_size(data))

    app_ui = ui.page_fluid(
        ui.div(
            ui.input_select("x", label="Symbol", choices=list(data.columns.values[1:])),
            ui.input_select("color", label="Color", choices=list(data.columns.values[1:])),
            ui.input_selectize(id="search", label="searchbar", choices=list(main_df.columns.values), width="600px"),
            class_="d-flex gap-3",
        ),
        output_widget("my_widget", width=600),
        ui.output_table("table_search"),
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
        def table_search():
            column = main_df.columns.get_loc(input.search())
            df = main_df.iloc[:, [0, column]]
            return df

    app = App(app_ui, server)
    run_app(app=app, launch_browser=True)
