import warnings

from shiny import App
from shiny import reactive
from shiny import render
from shiny import run_app
from shiny import ui

from m2m_postaviz.bin_exploration_module import bin_exp_server
from m2m_postaviz.bin_exploration_module import bin_exp_ui
from m2m_postaviz.cpd_exploration_module import cpd_tab_server
from m2m_postaviz.cpd_exploration_module import cpd_tab_ui
from m2m_postaviz.data_struct import DataStorage
from m2m_postaviz.overview_module import overview_module_server
from m2m_postaviz.overview_module import overview_module_ui


def run_shiny(data: DataStorage):

    warnings.filterwarnings("ignore", category=FutureWarning, module="plotly.express")

    metadata_label = data.get_factors(with_dtype=True)

    if data.HAS_TAXONOMIC_DATA:
        taxonomic_rank = data.get_taxonomy_rank()
        taxonomic_rank.insert(0, "all")

    metadata_table = ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select("metadata_factor", "Current column: ", metadata_label, selected=metadata_label[0]),
                ui.input_select("metadata_dtype", "dtype: ", ["category", "str", "int", "float"]),
                ui.input_action_button("dtype_change", "Update", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                bg="lightgrey"
               ),
        ui.output_text_verbatim("update_metadata_log", True),
        ui.output_data_frame("metadata_table")
        ),full_screen=True, min_height="800px")

    ### APPLICATION TREE ###

    app_ui = ui.page_fillable(
        ui.navset_tab(
            ui.nav_panel("Overview",
                overview_module_ui("module_pcoa", data)
            ),

            ui.nav_panel("Taxonomy-based exploration",
                bin_exp_ui("module_bin_exp", data)
            ),

            ui.nav_panel("Metabolite-based exploration",
                cpd_tab_ui("module_cpd_exp", data)
            ),

            ui.nav_panel("Metadata management",
                metadata_table
                ),

            ),
        )

    def server(input, output, session):

        cpd_tab_server("module_cpd_exp", Data=data)

        bin_exp_server("module_bin_exp", data)

        overview_module_server("module_pcoa", data)

        @render.data_frame
        def metadata_table():
            df = data.get_metadata()
            return df

        @render.text()
        @reactive.event(input.dtype_change)
        def update_metadata_log():

            factor_choice, dtype_choice = input.metadata_factor().split("/")[0], input.metadata_dtype()

            df = data.get_metadata().to_pandas()

            try:
                df[factor_choice] = df[factor_choice].astype(dtype_choice)

                text = f"Column {factor_choice} changed to {dtype_choice}."
                data.set_metadata(df)

            except ValueError as e:

                text = f"Cannot perform changes, {e}"

            new_metadata_label = data.get_factors(with_dtype=True)
            ui.update_select("metadata_factor", choices=new_metadata_label)

            return text

        @render.text
        def no_taxonomy():
            return "No taxonomic data provided."

    app = App(app_ui, server)
    run_app(app=app, launch_browser=True)
