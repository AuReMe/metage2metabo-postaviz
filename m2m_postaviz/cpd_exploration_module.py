import polars as pl
from shiny import module
from shiny import reactive
from shiny import render
from shiny import ui
from shinywidgets import output_widget
from shinywidgets import render_widget
from faicons import icon_svg as icon

import m2m_postaviz.shiny_module as sm
from m2m_postaviz.data_struct import DataStorage


@module.ui
def cpd_tab_ui(Data: DataStorage):

    if Data.USE_METACYC_PADMET:

        cpd_exploration_all_category = Data.get_metacyc_category_list()

        cpd_exploration_all_category.insert(0, Data.get_outsider_cpd()[1])

        cpd_exploration_all_category.insert(0,"None")

    else:

        cpd_exploration_all_category = ["None"]

    welcome_card = ui.card(ui.output_ui("Starting_message"))

    compounds_exploration_card = ui.card(
    ui.card_header("Metabolites exploration.",
        ui.tooltip(
            ui.input_action_button("_cpd_tab", " ", icon=icon("circle-question")), "If you provided Metacyc database information as input, you can select metabolite families to be explored. If not, select one or several metabolites to check their producibility. A metadata variable can be used to group samples samples, and another one can be used to additionally color the barplot. You can also filter samples (include, exclude) based on a metadata variable.")),  
    ui.card_body(
        ui.layout_sidebar(
            ui.sidebar(

                ui.input_selectize("first_category_input","Select any compound category from the Metacyc database.\n(Compounds in data / Compounds in metacyc database)",choices=cpd_exploration_all_category, multiple=False),
                ui.input_selectize("category_choice_input", "Select a sub category of the Metacyc database", choices=cpd_exploration_all_category, selected=cpd_exploration_all_category[0], multiple=False, width="400px"),

                ui.card(" ",
                    ui.input_selectize("compounds_choice_input", "Select compounds", choices=Data.get_compound_list(without_compartment=True), multiple=True, remove_button=True),
                    max_height="400px",
                    min_height="250px",
                    ),
                ui.input_select("metadata_filter_input", "Group by metadata", Data.get_factors(remove_smpl_col=True, insert_none=True)),

                ui.input_selectize("color_filter_input", "Second grouping for the barplot ", Data.get_factors(remove_smpl_col=False, insert_none=True), multiple=False, width="400px"),

                ui.card(" ",
                    ui.card_header("Sample Filtering (Optional)"),
                    ui.input_radio_buttons("sample_filter_choice_input", "Include or exclude samples from the analysis.", ["All", "Include", "Exclude"]),
                    ui.input_select("sample_filter_metadata1"," ",choices=Data.get_factors(True)),
                    ui.input_selectize("sample_filter_metadata2"," ",choices=[],multiple=True),
                    ui.input_selectize("sample_filter_selection_input", "Selection of samples to filter.", choices=Data.get_sample_list(), multiple=True, remove_button=True),
                    ui.card_footer(ui.input_action_button("reset_sample_filter_button", "Reset",width="50%")),

                    max_height="400px"),
                ui.row(
                    ui.input_checkbox("row_cluster", "Add rows clustering on Heatmap."),
                    ui.input_checkbox("col_cluster", "Add columns clustering on Heatmap."),
                ),
                ui.input_checkbox("render_cpd_abundance", "Use the abundance dataframe to weight the production of bins by their respective abundance in their sample."),

                ui.input_checkbox("exp_cpd_generate_stat_dataframe", "Generate statistical test dataframe."),

                ui.panel_conditional("input.exp_cpd_generate_stat_dataframe",

                    ui.input_checkbox("multiple_correction", "Multiple test correction"),
                    ui.panel_conditional("input.multiple_correction",
                                            ui.input_select("multiple_correction_method","Method",Data.get_list_of_tests(),selected=Data.get_list_of_tests()[0],)),
                    ),

                ui.input_task_button("run_plot_generation","Generate plots"),

            width=400,
            gap=35,
            bg="lightgrey"
        ),

        ui.navset_card_tab(
            ui.nav_panel("Community metabolic potential", ui.card(ui.output_plot("heatmap_cscope"), ui.card_footer(ui.input_action_button("save_cscope_heatmap", "save plot."), ui.output_text_verbatim("log_cscope_save")), full_screen=True)),
            ui.nav_panel("Individual metabolic potential", ui.card(ui.output_plot("heatmap_iscope"), ui.card_footer(ui.input_action_button("save_iscope_heatmap", "save plot."), ui.output_text_verbatim("log_iscope_save")), full_screen=True)),
            ui.nav_panel("Metabolite producible through cooperation only", ui.card(ui.output_plot("heatmap_added_value"), ui.card_footer(ui.input_action_button("save_advalue_heatmap", "save plot."), ui.output_text_verbatim("log_advalue_save")), full_screen=True)),
            title= "Heatmap depicting the number of metabolic networks able to produce the selected metabolites in samples"),

        ui.navset_card_tab(
            ui.nav_panel("Community metabolic potential", ui.card(output_widget("sample_percentage_production_cscope"), full_screen=True)),
            ui.nav_panel("Individual metabolic potential", ui.card(output_widget("sample_percentage_production_iscope"), full_screen=True)),
            title= "Barplot showing the percentage of samples having at least one metabolic network able to produce the metabolites, either individually or considering interactions across populations."),

        ui.card(ui.card_header("Boxplot showing the number of metabolic network producers for selected metabolites"),
            output_widget("cpd_exp_producers_plot"),full_screen=True),

        ui.card(ui.output_data_frame("cpd_exp_stat_dataframe"),full_screen=True),

    ),

    min_height="1500px"
    ),

    full_screen=True,)

    return welcome_card, compounds_exploration_card


@module.server
def cpd_tab_server(input, output, session, Data: DataStorage):

    @render.ui
    def Starting_message():
        msg = (
            '<div style="white-space: normal;">'
            "This is the <i><b>Metabolite-based exploration</b></i> tab.<br>"
            "Here, you can explore the producibility of certain metabolites of interest across samples of the dataset.<br>"
            "You can select metabolites and check whether they can be produced by metabolic networks individually, or collectively through metabolic interactions occurring across populations.<br><br>"
            "<i> A few tips: </i><br>"
            "<ul>"
                "<li>You can weigh the producibility values of metabolites by the relative abundance of the producer instead of using {0,1} values.</li>"
                "<li>If your data refers to the Metacyc database and if you have provided a padmet file (see documentation), you can also select a category of interest and explore the compounds in that category rather than picking metabolites one by one.</li>"
            "</ul>"
            "<i>Several plots are generated: </i><br>"
            "<ul>"
                "<li>The first ones are heatmaps displaying the number of metabolite producers in each sample, both at the community level (Cscope) and at the individual population level (Iscope), differences between both being the 'added-value of cooperation'. You can navigate between these three plots in the tabs below.</li>"
                "<li>A barplot shows the percentage of samples having at least one metabolic network producing the metabolites, either individually, or taking into account interactions across populations.<br> </li>"
                "<li>The last plot is a boxplot showing the number of producers for each metabolite, which can be filtered by metadata and colored by a variable of interest.</li>"
            "</ul>"
            'If you have any questions, please refer to the online '
            '<a href="https://metage2metabo-postaviz.readthedocs.io/en/latest/reference/m2m_postaviz.html" target="_blank">documentation</a> '
            'or raise an issue on <a href="https://github.com/AuReMe/metage2metabo-postaviz/tree/main" target="_blank" > GitHub</a>.<br>'
            '</div>'
        )
        return ui.HTML(msg)

    @render.text
    @reactive.event(input.save_cscope_heatmap, ignore_none=True, ignore_init=True)
    def log_cscope_save():
        obj_to_save = Data.get_working_dataframe("cscope_heatmap")
        if obj_to_save is None:
            return "Obj to save is None."
        else:
            log_msg = Data.save_seaborn_plot(obj_to_save, "cscope_heatmap.pdf")
            return log_msg

    @render.text
    @reactive.event(input.save_iscope_heatmap, ignore_none=True, ignore_init=True)
    def log_iscope_save():
        obj_to_save = Data.get_working_dataframe("iscope_heatmap")
        if obj_to_save is None:
            return "Obj to save is None."
        else:
            log_msg = Data.save_seaborn_plot(obj_to_save, "iscope_heatmap.pdf")
            return log_msg

    @render.text
    @reactive.event(input.save_advalue_heatmap, ignore_none=True, ignore_init=True)
    def log_advalue_save():
        obj_to_save = Data.get_working_dataframe("advalue_heatmap")
        if obj_to_save is None:
            return
        else:
            log_msg = Data.save_seaborn_plot(obj_to_save, "advalue_heatmap.pdf")
            return log_msg

    @reactive.effect
    @reactive.event(input.reset_sample_filter_button, ignore_none=True)
    def _reset_sample_filter_choice():
        return ui.update_selectize("sample_filter_selection_input", choices=Data.get_sample_list(), selected=None)

    @render.plot
    def heatmap_cscope():
        plot_object = cpd_plot_generation.result()[3][0]
        Data.keep_working_dataframe("cscope_heatmap", plot_object)
        return plot_object

    @render.plot
    def heatmap_iscope():
        plot_object = cpd_plot_generation.result()[3][1]
        Data.keep_working_dataframe("iscope_heatmap", plot_object)
        return plot_object

    @render.plot
    def heatmap_added_value():
        plot_object = cpd_plot_generation.result()[3][2]
        Data.keep_working_dataframe("advalue_heatmap", plot_object)
        return plot_object

    @render.data_frame
    def cpd_exp_stat_dataframe():
        return cpd_plot_generation.result()[2]

    @render_widget
    def sample_percentage_production_cscope():
        try:
            plot = cpd_plot_generation.result()[1][0]
        except TypeError:
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
    async def cpd_plot_generation(selected_compounds, user_input1, user_color_input, sample_filter_mode, sample_filter_value, with_statistic, with_multiple_correction, multiple_correction_method, row_cluster, col_cluster, render_cpd_abundance):

        cpd_filtered_list = []
        for cpd in Data.get_compound_list():
            if cpd[:-3] in selected_compounds:
                cpd_filtered_list.append(cpd)

        if len(selected_compounds) == 0:
            return

        nb_producers_boxplot = sm.render_reactive_metabolites_production_plot(Data, cpd_filtered_list, user_input1, user_color_input, sample_filter_mode, sample_filter_value, render_cpd_abundance) ###

        if user_input1 != "None":

            percent_barplot = sm.percentage_smpl_producing_cpd(Data, cpd_filtered_list, user_input1, sample_filter_mode, sample_filter_value)

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
                            input.multiple_correction_method(), input.row_cluster(), input.col_cluster(), input.render_cpd_abundance())


    @reactive.effect
    def _update_sub_category_choices():

        if not Data.USE_METACYC_PADMET:

            return

        metacyc_category_first_input = input.first_category_input().split(" ")[0]

        if metacyc_category_first_input == "Others":

            return ui.update_selectize("category_choice_input", choices=["Others"])

        if metacyc_category_first_input == "None" or metacyc_category_first_input == "":

            return ui.update_selectize("category_choice_input", choices=[])

        category_node = []

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

        if category_level == "Others":

            return ui.update_selectize("compounds_choice_input", choices=Data.get_outsider_cpd()[0])

        category_node = []

        Data.get_sub_tree_recursive(Data.get_cpd_category_tree(), category_level, category_node)

        category_node = category_node[0]

        cpds_found = []

        Data.get_compounds_from_category(category_node, cpds_found)

        cpd_list = Data.get_compound_list(True)

        final_cpd_list = [cpd for cpd in cpd_list if cpd in cpds_found]

        return ui.update_selectize("compounds_choice_input", choices=final_cpd_list, selected=final_cpd_list)

    @reactive.effect
    def _update_metadata_sample_filter():

        metadata_input1 = input.sample_filter_metadata1()

        try:
            metadata = Data.get_metadata().get_column(metadata_input1).unique().to_list()
        except Exception as e:
            ui.notification_show(
            f"Sample filter metadata error, {e}",
            type="warning",
            duration=6,)
            return

        return ui.update_selectize("sample_filter_metadata2", choices=metadata)

    @reactive.effect
    def _update_sample_filter_selection_input():

        sample_metadata_filter1 = input.sample_filter_metadata1()
        sample_metadata_filter2 = input.sample_filter_metadata2()
        metadata = Data.get_metadata()

        try:
            metadata = metadata.filter(pl.col(sample_metadata_filter1).is_in(sample_metadata_filter2))
        except:
            metadata = metadata.with_columns(pl.col(sample_metadata_filter1)).cast(pl.String)
            metadata = metadata.filter(pl.col(sample_metadata_filter1).is_in(sample_metadata_filter2))

        metadata = metadata.get_column("smplID").to_list()

        return ui.update_selectize("sample_filter_selection_input", choices=metadata, selected = metadata)
