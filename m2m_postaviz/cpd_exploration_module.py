from shiny import module, ui
from shinywidgets import output_widget
from shinywidgets import render_widget
from m2m_postaviz.data_struct import DataStorage
from shiny import reactive
from shiny import render
import m2m_postaviz.shiny_module as sm


@module.ui
def cpd_tab_ui(Data: DataStorage):
        
    cpd_exploration_all_category = ["None", "Ontology level 1", "Ontology level 2", "Pathway level 1", "Pathway level 2"]

    compounds_exploration_card = ui.card(
    ui.card_header("Metabolites exploration."),
    ui.card_body(
        ui.layout_sidebar(
            ui.sidebar(

                ui.input_select("category_level_input", "Search by compounds, ontology or pathway", choices= cpd_exploration_all_category, selected= cpd_exploration_all_category[0]),

                ui.input_selectize("category_choice_input", "Select one category", choices=Data.get_compounds_category_list()[0], selected=Data.get_compounds_category_list()[0][0], multiple=False, width="400px"),
                
                ui.input_selectize("compounds_choice_input", "Select compounds", choices=Data.get_compound_list(without_compartment=True), multiple=True, remove_button=True),
            
                ui.input_select("metadata_filter_input", "Filter by metadata:", Data.get_factors(remove_smpl_col=True, insert_none=True)),

                ui.input_selectize("color_filter_input", "Add color by metadata: ", Data.get_factors(remove_smpl_col=False, insert_none=True), multiple=False, width="400px"),

                ui.input_checkbox("sample_filter_enable_input","Filter plot by sample."),

                ui.panel_conditional("input.sample_filter_enable_input",

                    ui.input_radio_buttons("sample_filter_choice_input", "Exclude or include samples into the plot.", ["Include", "Exclude"]),
                    ui.input_selectize("sample_filter_selection_input", "Selection of samples to filter.", choices=Data.get_sample_list(), multiple=True, remove_button=True),
                    ui.input_action_button("reset_sample_filter_button", "Reset")
                    ),

                ui.input_task_button("run_plot_generation","Generate plots"),

                ##### stat ui

                ui.input_checkbox("exp_cpd_generate_stat_dataframe", "Generate statistical test dataframe."),

                ui.panel_conditional("input.exp_cpd_generate_stat_dataframe",
                                        
                    ui.input_checkbox("exp_cpd_multiple_correction", "Multiple test correction"),
                    ui.panel_conditional("input.exp_cpd_multiple_correction",
                                            ui.input_select("exp_cpd_multiple_correction_method","Method",Data.get_list_of_tests(),selected=Data.get_list_of_tests()[0],)),
                    ),

            width=400,
            gap=35,
            bg="lightgrey"
        ),
        ui.card(
            ui.input_radio_buttons("heatmap_radio_button", " ", ["Cscope", "Iscope", "Added value"], selected="Cscope"),
            output_widget("cpd_exp_diff_heatmap"),full_screen=True
            ),

        ui.card(
            ui.input_radio_buttons("percent_barplot_radio_button", " ", ["Cscope", "Iscope"], selected="Cscope"),
            output_widget("cpd_exp_percent"),full_screen=True
            ),

        ui.card(output_widget("cpd_exp_producers_plot"),full_screen=True),

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
                
        @render_widget
        def cpd_exp_diff_heatmap():

            if input.heatmap_radio_button() == "Cscope":

                return cpd_plot_generation.result()[3][0]
            
            if input.heatmap_radio_button() == "Iscope":

                return cpd_plot_generation.result()[3][1]
    
            if input.heatmap_radio_button() == "Added value":

                return cpd_plot_generation.result()[3][2]

        @render.data_frame
        def cpd_exp_stat_dataframe():
            return cpd_plot_generation.result()[2]

        @render_widget
        def cpd_exp_percent():
            if input.percent_barplot_radio_button() == "Cscope":
                
                return cpd_plot_generation.result()[1][0]

            if input.percent_barplot_radio_button() == "Iscope":
                
                return cpd_plot_generation.result()[1][1]
                

        @render_widget
        def cpd_exp_producers_plot():
            return cpd_plot_generation.result()[0]

        @ui.bind_task_button(button_id="run_plot_generation")
        @reactive.extended_task
        async def cpd_plot_generation(selected_compounds, user_input1, user_color_input, sample_filtering_enabled, sample_filter_mode, sample_filter_value, with_statistic, with_multiple_correction, multiple_correction_method):

            cpd_filtered_list = []
            for cpd in Data.get_compound_list():
                if cpd[:-3] in selected_compounds:
                    cpd_filtered_list.append(cpd)

            if len(selected_compounds) == 0:
                return

            nb_producers_boxplot = sm.render_reactive_metabolites_production_plot(Data, cpd_filtered_list, user_input1, user_color_input, sample_filtering_enabled, sample_filter_mode, sample_filter_value)

            if user_input1 != "None":

                percent_barplot = sm.percentage_smpl_producing_cpd(Data, cpd_filtered_list, user_input1,)

            else:
                percent_barplot = None

            if with_statistic:

                try:
                    stat_dataframe = sm.metabolites_production_statistical_dataframe(Data, cpd_filtered_list, user_input1, "None", with_multiple_correction, multiple_correction_method)
                except:
                    stat_dataframe = None
            else:

                stat_dataframe = None

            cscope_heatmap, iscope_heatmap, added_value_heatmap = sm.added_value_heatmap(Data, cpd_filtered_list, sample_filtering_enabled, sample_filter_mode, sample_filter_value)

            return nb_producers_boxplot, percent_barplot, stat_dataframe, (cscope_heatmap, iscope_heatmap, added_value_heatmap)


        @reactive.effect
        @reactive.event(input.run_plot_generation, ignore_none=True)
        def handle_click_cpd_exploration():
            cpd_plot_generation(input.compounds_choice_input(), input.metadata_filter_input(), 
                                input.color_filter_input(), input.sample_filter_enable_input(),
                                input.sample_filter_choice_input(), input.sample_filter_selection_input(), 
                                input.exp_cpd_generate_stat_dataframe(), input.exp_cpd_multiple_correction(),
                                input.exp_cpd_multiple_correction_method())


        @reactive.effect
        def _update_category_choices():

            category_level = input.category_level_input()

            if category_level == "None":

                return ui.update_selectize("category_choice_input", choices=["None"], selected="None")

            if category_level == "Ontology level 1":

                return ui.update_selectize("category_choice_input", choices=Data.get_compounds_category_list()[0], selected=Data.get_compounds_category_list()[0][0])
                
            if category_level == "Ontology level 2":

                return ui.update_selectize("category_choice_input", choices=Data.get_compounds_category_list()[1], selected=Data.get_compounds_category_list()[1][0])

            if category_level == "Pathway level 1":

                return ui.update_selectize("category_choice_input", choices=[])

            if category_level == "Pathway level 1":

                return ui.update_selectize("category_choice_input", choices=[])
            

        @reactive.effect
        def _update_compounds_choices():

            category_level = input.category_level_input()
            category_choice = input.category_choice_input()
            cpd_list_no_compart = Data.get_compound_list(without_compartment=True)

            if category_level == "None":

                return ui.update_selectize("compounds_choice_input", choices=Data.get_compound_list(True))

            if category_level == "Ontology level 1":

                df = Data.get_compounds_category_dataframe()[0]

            if category_level == "Ontology level 2":

                df = Data.get_compounds_category_dataframe()[1]

            if category_level == "Pathway level 1":

                return ui.update_selectize("compounds_choice_input", choices=[])

            if category_level == "Pathway level 2":

                return ui.update_selectize("compounds_choice_input", choices=[])

            df = df[df['category'] == category_choice]
            input_choices = df.loc[df["compound_id"].isin(cpd_list_no_compart)]["compound_id"].tolist()

            return ui.update_selectize("compounds_choice_input", choices=input_choices, selected=input_choices)
        