import plotly.express as px
from shiny import module
from shiny import reactive
from shiny import render
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.data_utils as du
import m2m_postaviz.shiny_module as sm
from m2m_postaviz.data_struct import DataStorage


@module.ui
def pcoa_module_ui(Data: DataStorage):

    metadata_label = Data.get_factors()
    metadata_label.remove("smplID")

    pcoa_card = ui.card(ui.card_header("Principal Coordinates Analysis made with all samples. Change color input to see different clusters."),
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select(id="pcoa_color", label="Plot color.", choices=metadata_label, selected=metadata_label[0]),
                ui.output_ui("pcoa_ui"),
                ui.help_text(ui.output_text("display_warning_pcoa")),
                width=300,
                gap=30,
                bg="lightgrey",
                
            ),
        output_widget("pcoa_plot")
        ),
        full_screen=True,
        min_height="600px"
        )


    custom_pcoa_card = ui.card(ui.card_header("Customize the Principal Coordinates Analysis by filtering samples with their metadata value."),
        ui.layout_sidebar(
            ui.sidebar(

                ui.input_select(id="custom_pcoa_selection", label="Choose column", choices=metadata_label, selected=metadata_label[0]),
                ui.output_ui("pcoa_custom_choice"),
                ui.input_select(id="pcoa_custom_color", label="Color.", choices=metadata_label, selected=metadata_label[0]),
                ui.input_checkbox("pcoa_custom_abundance_check", "Use abundance data."),

                ui.input_task_button("run_custom_pcoa","Go"),
                width=300,
                gap=35,
                bg="lightgrey"
            ),
        output_widget("pcoa_custom_plot"),
        ),
        full_screen=True,
        min_height="600px"
    )

    return pcoa_card, custom_pcoa_card

@module.server
def pcoa_module_server(input, output, session, Data: DataStorage):

    @render.text
    def display_warning_pcoa():
        return "Warning: this is just for display, Pcoa's dataframe is not recalculated."

    @reactive.effect
    @reactive.event(input.run_custom_pcoa, ignore_none=False)
    def handle_click():
        make_custom_pcoa(input.custom_pcoa_selection(), input.pcoa_custom_choice(), input.pcoa_custom_abundance_check(), input.pcoa_custom_color())

    @ui.bind_task_button(button_id="run_custom_pcoa")
    @reactive.extended_task
    async def make_custom_pcoa(column, choices, abundance, color):

        return sm.make_pcoa(Data, column, choices, abundance, color)

    @render_widget
    def pcoa_custom_plot():
        return make_custom_pcoa.result()

    @render.ui
    @reactive.event(input.custom_pcoa_selection)
    def pcoa_custom_choice():

        df = Data.get_metadata()
        value = df[input.custom_pcoa_selection()]

        if not du.serie_is_float(value):

            return ui.TagList(
                    ui.card(" ",
                        ui.input_selectize("pcoa_custom_choice", "Get only in column:", value.unique().tolist(),
                                            selected=value.unique().tolist(),
                                            multiple=True,
                                            options={"plugins": ["clear_button"]})
                    ,max_height="400px"
                    ))
        else:

            return ui.TagList(
                        ui.input_slider("pcoa_custom_choice", "Set min/max filter:", min=value.min(), max=value.max(), value=[value.min(),value.max()]),)

    @render_widget
    def pcoa_plot():
        # Get all parameters.
        selected_col = input.pcoa_color()

        df = Data.get_pcoa_dataframe()

        # Check column dtype.
        if du.serie_is_float(df[selected_col]):

            min_limit = input.pcoa_slider()[0]

            max_limit = input.pcoa_slider()[1]

            df = df.loc[(df[selected_col] <= max_limit) & (df[selected_col] >= min_limit)]

        else:

            show_only = input.pcoa_selectize()

            df = df.loc[df[selected_col].isin(show_only)]

        fig = px.scatter(df, x="PC1", y="PC2", color=selected_col)

        return fig

    @render.ui
    @reactive.event(input.pcoa_color)
    def pcoa_ui():

        df = Data.get_pcoa_dataframe()
        value = df[input.pcoa_color()]

        if not du.serie_is_float(value):

            return ui.TagList(
                    ui.card(" ",
                        ui.input_selectize("pcoa_selectize", "Show only:", value.unique().tolist(),
                                            selected=value.unique().tolist(),
                                            multiple=True,
                                            options={"plugins": ["clear_button"]}),
                        max_height="400px"
                        ))
        else:

            return ui.TagList(
                        ui.input_slider("pcoa_slider", "Set min/max filter:", min=value.min(), max=value.max(), value=[value.min(),value.max()]),)
