from shiny import module
from shiny import reactive
from shiny import render
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget

import m2m_postaviz.shiny_module as sm
from m2m_postaviz.data_struct import DataStorage


@module.ui
def cpd_tab_ui(Data: DataStorage):

    if Data.USE_METACYC_PADMET:

        cpd_exploration_all_category = Data.get_metacyc_category_list()

        cpd_exploration_all_category.insert(0,"None")

    else:

        cpd_exploration_all_category = ["None"]

    compounds_exploration_card = ui.card(
    ui.card_header("Metabolites exploration."),
    ui.card_body(
        ui.layout_sidebar(
            ui.sidebar(

                ui.input_select("first_category_input","Select any compounds category from Metacyc_database.\n(Compounds in data / Compounds in metacyc database)",choices=cpd_exploration_all_category),
                ui.input_selectize("category_choice_input", "Select a sub category", choices=cpd_exploration_all_category, selected=cpd_exploration_all_category[0], multiple=False, width="400px"),

                ui.card(" ",
                    ui.input_selectize("compounds_choice_input", "Select compounds", choices=Data.get_compound_list(without_compartment=True), multiple=True, remove_button=True),
                    max_height="400px",
                    min_height="250px",
                    ),
                ui.input_select("metadata_filter_input", "Filter by metadata:", Data.get_factors(remove_smpl_col=True, insert_none=True)),

                ui.input_selectize("color_filter_input", "Add color by metadata: ", Data.get_factors(remove_smpl_col=False, insert_none=True), multiple=False, width="400px"),

                    ui.input_radio_buttons("sample_filter_choice_input", "Include or exclude samples from the plots.", ["All", "Include", "Exclude"]),
                    ui.input_selectize("sample_filter_selection_input", "Selection of samples to filter.", choices=Data.get_sample_list(), multiple=True, remove_button=True),
                    ui.input_action_button("reset_sample_filter_button", "Reset",width="50%"),

                ui.row(
                    ui.input_checkbox("row_cluster", "Add rows clustering on Heatmap."),
                    ui.input_checkbox("col_cluster", "Add columns clustering on Heatmap."),
                ),

                ui.input_task_button("run_plot_generation","Generate plots"),

                ##### stat ui

                ui.input_checkbox("exp_cpd_generate_stat_dataframe", "Generate statistical test dataframe."),

                ui.panel_conditional("input.exp_cpd_generate_stat_dataframe",

                    ui.input_checkbox("multiple_correction", "Multiple test correction"),
                    ui.panel_conditional("input.multiple_correction",
                                            ui.input_select("multiple_correction_method","Method",Data.get_list_of_tests(),selected=Data.get_list_of_tests()[0],)),
                    ),

            width=400,
            gap=35,
            bg="lightgrey"
        ),

        ui.navset_card_tab(
            ui.nav_panel("Cscope", ui.card(ui.output_plot("heatmap_cscope"), full_screen=True)),
            ui.nav_panel("Iscope", ui.card(ui.output_plot("heatmap_iscope"), full_screen=True)),
            ui.nav_panel("Added value", ui.card(ui.output_plot("heatmap_added_value"), full_screen=True)),
            title= "Heatmap of the number of bins producers of the compounds for each samples (or selected samples). Metadata filtering and hierarchical clustering on both samples and compounds is available."),

        ui.navset_card_tab(
            ui.nav_panel("Cscope", ui.card(output_widget("sample_percentage_production_cscope"), full_screen=True)),
            ui.nav_panel("Iscope", ui.card(output_widget("sample_percentage_production_iscope"), full_screen=True)),
            title= "Barplot showing the percentage of sample producing the compounds (at least one genomes producers). Both in Cscope and Iscope."),

        ui.card(ui.card_header("Boxplot of the numbers of genomes producers (Y-axis) for each compounds (X-axis) in input. Can be filtered by metadata and grouped by color input."),
            output_widget("cpd_exp_producers_plot"),full_screen=True),

        ui.card(ui.output_data_frame("cpd_exp_stat_dataframe"),full_screen=True),

    ),

    min_height="1500px"
    ),

    full_screen=True,)

    return compounds_exploration_card


@module.server
def cpd_tab_server(input, output, session, Data: DataStorage):

        @reactive.effect
        @reactive.event(input.reset_sample_filter_button, ignore_none=True)
        def _reset_sample_filter_choice():
            return ui.update_selectize("sample_filter_selection_input", choices=Data.get_sample_list(), selected=None)

        @render.plot
        def heatmap_cscope():
            return cpd_plot_generation.result()[3][0]

        @render.plot
        def heatmap_iscope():
            return cpd_plot_generation.result()[3][1]

        @render.plot
        def heatmap_added_value():
            return cpd_plot_generation.result()[3][2]

        @render.data_frame
        def cpd_exp_stat_dataframe():
            return cpd_plot_generation.result()[2]

        @render_widget
        def sample_percentage_production_cscope():
            try:
                plot = cpd_plot_generation.result()[1][0]
            except TypeError as e:
                ui.notification_show(
                "Sample Percentage production plot needs a metadata filter input.",
                type="warning",
                duration=6,)
                plot = None

            return plot

        @render_widget
        def sample_percentage_production_iscope():
            return cpd_plot_generation.result()[1][1]

        @render_widget
        def cpd_exp_producers_plot():
            return cpd_plot_generation.result()[0]

        @ui.bind_task_button(button_id="run_plot_generation")
        @reactive.extended_task
        async def cpd_plot_generation(selected_compounds, user_input1, user_color_input, sample_filter_mode, sample_filter_value, with_statistic, with_multiple_correction, multiple_correction_method, row_cluster, col_cluster):

            cpd_filtered_list = []
            for cpd in Data.get_compound_list():
                if cpd[:-3] in selected_compounds:
                    cpd_filtered_list.append(cpd)

            if len(selected_compounds) == 0:
                return

            nb_producers_boxplot = sm.render_reactive_metabolites_production_plot(Data, cpd_filtered_list, user_input1, user_color_input, sample_filter_mode, sample_filter_value) ###

            if user_input1 != "None":

                percent_barplot = sm.percentage_smpl_producing_cpd(Data, cpd_filtered_list, user_input1,)

            else:
                percent_barplot = None

            if with_statistic:

                try:
                    stat_dataframe = sm.metabolites_production_statistical_dataframe(Data, cpd_filtered_list, user_input1, "None", with_multiple_correction, multiple_correction_method)
                except:  # noqa: E722
                    stat_dataframe = None
            else:

                stat_dataframe = None

            cscope_heatmap, iscope_heatmap, added_value_heatmap = sm.sns_clustermap(Data, cpd_filtered_list, user_input1, row_cluster, col_cluster, sample_filter_mode, sample_filter_value) ###

            return nb_producers_boxplot, percent_barplot, stat_dataframe, (cscope_heatmap, iscope_heatmap, added_value_heatmap)


        @reactive.effect
        @reactive.event(input.run_plot_generation, ignore_none=True)
        def handle_click_cpd_exploration():
            cpd_plot_generation(input.compounds_choice_input(), input.metadata_filter_input(),
                                input.color_filter_input(),
                                input.sample_filter_choice_input(), input.sample_filter_selection_input(),
                                input.exp_cpd_generate_stat_dataframe(), input.multiple_correction(),
                                input.multiple_correction_method(), input.row_cluster(), input.col_cluster())


        @reactive.effect
        def _update_sub_category_choices():

            if not Data.USE_METACYC_PADMET:

                return

            metacyc_category_first_input = input.first_category_input().split(" ")[0]

            if metacyc_category_first_input == "None" or metacyc_category_first_input == "":

                return ui.update_selectize("category_choice_input", choices=[])

            category_node = []
            print(metacyc_category_first_input)
            Data.get_sub_tree_recursive(Data.get_cpd_category_tree(), metacyc_category_first_input, category_node)

            category_node = category_node[0]

            new_sub_category_list = Data.get_metacyc_category_list(category_node)

            new_sub_category_list.insert(0, input.first_category_input())

            return ui.update_selectize("category_choice_input", choices=new_sub_category_list)

        @reactive.effect
        def _update_compounds_choices():

            if not Data.USE_METACYC_PADMET:

                return

            category_level = input.category_choice_input().split(" ")[0]

            if category_level == "":

                return ui.update_selectize("compounds_choice_input", choices=[])

            if category_level == "None":

                return ui.update_selectize("compounds_choice_input", choices=Data.get_compound_list(without_compartment=True))

            category_node = []

            Data.get_sub_tree_recursive(Data.get_cpd_category_tree(), category_level, category_node)

            category_node = category_node[0]

            cpds_found = []

            Data.find_compounds_from_category(category_node, cpds_found)

            cpd_list = Data.get_compound_list(True)

            final_cpd_list = [cpd for cpd in cpd_list if cpd in cpds_found]

            return ui.update_selectize("compounds_choice_input", choices=final_cpd_list, selected=final_cpd_list)

