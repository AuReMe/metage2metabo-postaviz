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
def bin_exp_ui(Data: DataStorage):

    factor_list = Data.get_factors()
    factor_list.insert(0, "None")

    if Data.HAS_TAXONOMIC_DATA:
        taxonomic_rank = Data.get_taxonomy_rank()
        taxonomic_rank.insert(0, "all")

        welcome_card = ui.card(ui.output_text("Starting_message"))

        bins_exploration_card = ui.card(
            ui.card_header("Bins exploration",
                ui.tooltip(
                    ui.input_action_button("_bin_tab", " ", icon=icon("circle-question")), "blablabla tooltips !")),  #FC BLABLA HERE !!
            ui.card_body(
                ui.layout_sidebar(
                    ui.sidebar(

                        ui.input_selectize("rank_choice", "Choose a taxonomic rank.", taxonomic_rank, selected=taxonomic_rank[0], multiple=False),
                        ui.output_ui("rank_unique_choice"),

                        ui.input_selectize("bin_factor", "Filter", factor_list, selected=factor_list[0], multiple=False),
                        ui.output_ui("bin_factor_unique"),

                        ui.input_selectize("bin_color", "Color", factor_list, selected=factor_list[0], multiple=False),

                        ui.input_checkbox("with_bin_abundance", "Weigh the producibility value by the relative abundance of the producer instead of using {0,1} values."),

                        ui.input_task_button("run_bin_exploration","Go"),

                        ui.output_text("bin_size_text"),
                        ui.output_text("Timer_info"),

                    width=400,
                    gap=35,
                    bg="lightgrey"
                ),
                ui.card(ui.card_header("Maximum production of unique compounds by taxa/metagenomes in community."),output_widget("bin_unique_count_cscope_histplot"),full_screen=True),
                ui.card(ui.card_header("Maximum production of unique compounds by taxa/metagenomes as individual."),output_widget("bin_unique_count_iscope_histplot"),full_screen=True),
                # ui.card(ui.card_header(),output_widget("bin_boxplot_count"),full_screen=True),
                ui.card(ui.card_header("Respective abundance of taxa/metagenomes in their sample."),output_widget("bin_abundance_plot"),full_screen=True),
            )
        ),full_screen=True)

    else:

        bins_exploration_card = ui.output_text_verbatim("no_taxonomy_provided")

    return welcome_card, bins_exploration_card


@module.server
def bin_exp_server(input, output, session, Data: DataStorage):

    @render.text
    def Starting_message():
        msg = "blablabla"
        return msg

    @render.text
    def no_taxonomy_provided():
        return "No taxonomy file has been provided. Bins exploration disabled.\n You can use the -t option to provide a taxonomy file in txt or tsv/csv format. Metagenomes column must be the first column."

    list_of_bins = Data.get_bins_list()

    if Data.HAS_TAXONOMIC_DATA:
        taxonomic_rank = Data.get_taxonomy_rank()
        taxonomic_rank.insert(0, "all")
        converted_bin_list = Data.associate_bin_taxonomy(list_of_bins)

    @render.ui
    def bin_factor_unique():

        factor_choice = input.bin_factor()

        if factor_choice == "None":
            return ui.TagList(
            ui.input_selectize("bin_factor_unique", "Select", choices=[], multiple=True, remove_button=True)
            )

        df = Data.get_metadata()

        choices = df.get_column(factor_choice).unique().to_list()

        return ui.TagList(
            ui.input_selectize("bin_factor_unique", "Select", choices=choices, multiple=True, remove_button=True)
            )

    @render.ui
    def rank_unique_choice():

        rank_choice = input.rank_choice()

        if rank_choice == "all":

            return ui.TagList(
            ui.input_selectize("rank_unique_choice", "Select", choices=[], multiple=False,)
            )

        if rank_choice == "mgs":

            return ui.TagList(
            ui.input_selectize("rank_unique_choice", "Select", choices=converted_bin_list, multiple=False,)
            )

        df = Data.get_taxonomic_dataframe()

        choices = df.get_column(rank_choice).unique().to_list()

        return ui.TagList(
            ui.input_selectize("rank_unique_choice", "Select", choices=choices, multiple=False,)
            )

    @render.text
    def bin_size_text():
        """Display the number of bins within the selection of user's input.
        If only one bin is in selection, then display in how many samples this bin is present.

        Returns:
            str: str to display inside an output text UI.
        """
        rank_choice, rank_unique_choice = input.rank_choice(), input.rank_unique_choice()

        if rank_choice == "mgs":

            rank_unique_choice = rank_unique_choice.split(" ")[0]

        if rank_choice == "all":

            return f"All {len(list_of_bins)} bins selected"

        list_of_bin_in_rank = Data.get_bin_list_from_taxonomic_rank(rank_choice, rank_unique_choice)

        filtered_list_of_bin = []

        for x in list_of_bin_in_rank:
            if x in list_of_bins:

                filtered_list_of_bin.append(x)

        if len(filtered_list_of_bin) == 0:
            return "No bin in selection."

        if len(filtered_list_of_bin) == 1:
            return "Bin: "+str(filtered_list_of_bin[0])+" Found in "+str(Data.get_bins_count()[filtered_list_of_bin[0]])+" sample(s)"

        return f"{len(filtered_list_of_bin)} bins found in selection."

    @render.text
    def Timer_info():
        timer = run_exploration.result()[3]
        return f"Took {timer} seconds to run."

    @render_widget
    def bin_boxplot_count():
        return run_exploration.result()[4]

    @render_widget
    def bin_abundance_plot():
        return run_exploration.result()[1]

    @render_widget
    def bin_unique_count_cscope_histplot():
        return run_exploration.result()[0][0]

    @render_widget
    def bin_unique_count_iscope_histplot():
        return run_exploration.result()[0][1]

    # @ui.bind_task_button(button_id="run_custom_pcoa")
    @reactive.extended_task
    async def run_exploration(factor, factor_choice, rank, rank_choice, with_abundance, color):

        return sm.bin_exploration_processing(Data, factor, factor_choice, rank, rank_choice, with_abundance, color)

    @reactive.effect
    @reactive.event(input.run_bin_exploration, ignore_none=True)
    def handle_click_bin_exploration():
        run_exploration(input.bin_factor(), input.bin_factor_unique(), input.rank_choice(), input.rank_unique_choice(), input.with_bin_abundance(), input.bin_color())
