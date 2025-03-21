import warnings

from shiny import App
from shiny import reactive
from shiny import render
from shiny import run_app
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.shiny_module as sm
from m2m_postaviz.bin_exploration_module import bin_exp_server
from m2m_postaviz.bin_exploration_module import bin_exp_ui
from m2m_postaviz.cpd_exploration_module import cpd_tab_server
from m2m_postaviz.cpd_exploration_module import cpd_tab_ui
from m2m_postaviz.data_struct import DataStorage
from m2m_postaviz.pcoa_module import pcoa_module_server
from m2m_postaviz.pcoa_module import pcoa_module_ui


def run_shiny(data: DataStorage):

    warnings.filterwarnings("ignore", category=FutureWarning, module="plotly.express")

    factor_list = data.get_factors()
    factor_list.insert(0, "None")

    factor_list_no_smpl = data.get_factors()
    factor_list_no_smpl.remove("smplID")
    factor_list_no_smpl.insert(0, "None")

    metadata_label = data.get_factors(with_dtype=True)

    if data.HAS_TAXONOMIC_DATA:
        taxonomic_rank = data.get_taxonomy_rank()
        taxonomic_rank.insert(0, "all")

    all_dataframe = {"global_production_test_dataframe": None, "global_production_plot_dataframe": None, "metabolites_production_test_dataframe": None, "metabolites_production_plot_dataframe": None}

    ### ALL CARD OBJECT TO BE ARRANGED ###

    total_production_plot = ui.card(
        ui.card_header("Total production of all compound, weighted with the abundance if provided."),
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select("prod_inputx1", "Label for X axis", factor_list),
                ui.input_select("prod_inputx2", "Label for 2nd X axis", factor_list),

                ui.input_checkbox("prod_norm", "Abundance data"),
                ui.input_checkbox("multiple_correction_global_plot", "Multiple test correction"),
                ui.panel_conditional("input.multiple_correction_global_plot",ui.input_select("multiple_test_method_global","Method",
                                                                                                     ["bonferroni","sidak","holm-sidak","holm","simes-hochberg","hommel","fdr_bh","fdr_by","fdr_tsbh","fdr_tsbky"],
                                                                                                     selected="bonferroni",)),

                ui.input_action_button("export_global_production_plot_dataframe_button", "Save plot dataframe"),
                ui.output_text_verbatim("export_global_production_plot_dataframe_txt", True),
                ui.input_action_button("export_global_production_test_button", "Export stats dataframe"),
                ui.output_text_verbatim("export_global_production_test_dataframe", True),
                width=350,
                gap=30,
            ),

            ui.card(ui.p(output_widget("total_production_plot")),full_screen=True),

            ui.card(ui.output_data_frame("production_test_dataframe"),full_screen=True)

        ),
        full_screen=True
        )

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
            ui.nav_panel("Exploration",
                ui.layout_column_wrap(

                    ui.value_box(
                        "Numbers of unique metabolic network:",
                        ui.output_text("unique_total_bin_count"),
                        theme="bg-gradient-indigo-purple",
                    ),

                    ui.value_box(
                        "Numbers of samples:",
                        ui.output_text("total_samples_count"),
                        theme="cyan",
                    ),

                    ui.value_box(
                        "Numbers of unique compounds produced:",
                        ui.output_text("total_unique_cpd"),
                        theme="bg-gradient-blue-purple",
                    ),

                fill=False,
                ),
                total_production_plot
                ),

            ui.nav_panel("Metadata",
                metadata_table
                ),

            ui.nav_panel("PCOA",
                pcoa_module_ui("module_pcoa", data)
            ),

            ui.nav_panel("Bins",
                bin_exp_ui("module_bin_exp", data)
            ),

            ui.nav_panel("Cpd",
                cpd_tab_ui("module_cpd_exp", data)
            )

            ),
        )

    def server(input, output, session):

        cpd_tab_server("module_cpd_exp", Data=data)

        bin_exp_server("module_bin_exp", data)

        pcoa_module_server("module_pcoa", data)

        @render.text
        def unique_total_bin_count():
            return data.get_total_unique_bins_count()

        @render.text
        def total_samples_count():
            return str(data.get_main_dataframe().shape[0])

        @render.text
        def total_unique_cpd():
            return str(data.get_main_dataframe().shape[1])

        @render.data_frame
        def production_test_dataframe():

            return sm.global_production_statistical_dataframe(data, input.prod_inputx1(),
                                                              input.prod_inputx2(),
                                                              input.multiple_correction_global_plot(),
                                                              input.multiple_test_method_global(),
                                                              input.prod_norm())

        @render.text
        @reactive.event(input.export_global_production_test_button)
        def export_global_production_test_dataframe():

            if bool(all_dataframe):
                if all_dataframe["global_production_test_dataframe"] is not None:
                   log = data.save_dataframe(all_dataframe["global_production_test_dataframe"], "global_production_test_dataframe")
                else:
                    log = "Unable to find any dataframe to save."

            return log

        @render.text
        @reactive.event(input.export_metabolites_test_button)
        def export_metabolites_test_dataframe_txt():

            if bool(all_dataframe):
                if all_dataframe["metabolites_production_test_dataframe"] is not None:
                   log = data.save_dataframe(all_dataframe["metabolites_production_test_dataframe"], "metabolites_production_test_dataframe")
                else:
                    log = "Unable to find any dataframe to save."

            return log

        @render.text
        @reactive.event(input.export_metabolites_plot_button)
        def export_metabolites_plot_dataframe_txt():

            if bool(all_dataframe):
                if all_dataframe["metabolites_production_plot_dataframe"] is not None:
                   log = data.save_dataframe(all_dataframe["metabolites_production_plot_dataframe"], "metabolites_production_plot_dataframe")
                else:
                    log = "Unable to find any dataframe to save."

            return log

        @render.text
        @reactive.event(input.export_global_production_plot_dataframe_button)
        def export_global_production_plot_dataframe_txt():

            if bool(all_dataframe):
                if all_dataframe["global_production_plot_dataframe"] is not None:
                    log = data.save_dataframe(all_dataframe["global_production_plot_dataframe"], "producer_plot_dataframe")
                else:
                    log = "Unable to find any dataframe to save."

            return log

        @render_widget
        def total_production_plot():

            return sm.render_reactive_total_production_plot(data, input.prod_inputx1(), input.prod_inputx2(), input.prod_norm())

        @render.data_frame
        def metadata_table():
            df = data.get_metadata()
            return df

        @render.text()
        @reactive.event(input.dtype_change)
        def update_metadata_log():
            
            factor_choice, dtype_choice = input.metadata_factor().split("/")[0], input.metadata_dtype()

            df = data.get_metadata()
            print(df[factor_choice].dtype)
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
