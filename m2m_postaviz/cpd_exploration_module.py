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

    # Build first-level options: top-level categories first, then all others (prefixed)
    def _first_selector_options():
        tree = Data.get_cpd_category_tree()
        # Unwrap artificial root
        unwrapped = tree
        if len(unwrapped) == 1 and all(isinstance(v, dict) for v in unwrapped.values()):
            (_, unwrapped) = next(iter(unwrapped.items()))

        # Top-level options (direct children)
        data_cpd_list = Data.get_compound_list(without_compartment=True)
        top_level = []
        for key, child in unwrapped.items():
            if not child:
                continue
            cpd_in_category = []
            Data.get_compounds_from_category(child, cpd_in_category)
            final_cpd_list = [cpd for cpd in data_cpd_list if cpd in cpd_in_category]
            top_level.append(f"{key} ({len(final_cpd_list)}/{len(cpd_in_category)})")

        # All categories (entire tree)
        all_with_counts = Data.get_metacyc_category_list()

        # Merge with top-level first, then the rest (excluding duplicates and the single artificial root if present)
        seen = set()
        merged = []
        for item in top_level:
            merged.append(item)
            seen.add(item.split(" (")[0])
        for item in all_with_counts:
            name = item.split(" (")[0]
            if name in ("All_metabolites",):
                continue
            if name not in seen:
                merged.append(f"-- {item}")  # prefix to visually separate deeper nodes
                seen.add(name)
        return merged

    if Data.USE_METACYC_PADMET:

        cpd_exploration_all_category = _first_selector_options()

        # Add global All/Others with counts for the first selector
        data_cpd_list = Data.get_compound_list(without_compartment=True)
        all_cpd_list = []
        Data.get_compounds_from_category(Data.get_cpd_category_tree(), all_cpd_list)
        outsiders_full = Data.get_outsider_cpd()[0]
        all_cpd_list = list(set(all_cpd_list + outsiders_full))
        all_in_data = len([cpd for cpd in all_cpd_list if cpd in data_cpd_list])
        all_label = f"All ({all_in_data}/{len(all_cpd_list)})"

        outsiders = outsiders_full
        outsiders_in_data = len([cpd for cpd in outsiders if cpd in data_cpd_list])
        outsiders_label = f"Others ({outsiders_in_data}/{len(all_cpd_list)})"

        cpd_exploration_all_category.insert(0, outsiders_label)
        cpd_exploration_all_category.insert(0, all_label)
        cpd_exploration_all_category.insert(0, "None")

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

                ui.input_selectize("first_category_input","Select any compound category from the Metacyc database (optional).\n(Compounds in data / Compounds in metacyc database)",choices=cpd_exploration_all_category, multiple=False),
                ui.markdown("**Available top-level categories:** " + ", ".join([c.rsplit(" (", 1)[0] for c in cpd_exploration_all_category if c not in ("None",) and not c.startswith("Others") and not c.startswith("-- ")])),
                ui.input_selectize("category_choice_input", "Select a sub category of the Metacyc database (optional)", choices=cpd_exploration_all_category, selected=cpd_exploration_all_category[0], multiple=False, width="400px"),
                ui.input_selectize("category_choice_input_lvl3", "Select a third-level category (optional)", choices=[], selected=None, multiple=False, width="400px"),

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
                ui.input_checkbox("save_raw_data", "Save dataframe used to generate plots."),
                ui.output_text_verbatim("save_raw_data_logs"),

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

        ui.navset_card_tab(
            ui.nav_panel("Community metabolic potential", ui.card(output_widget("cpd_exp_producers_plot"), full_screen=True)),
            ui.nav_panel("Individual metabolic potential", ui.card(output_widget("cpd_exp_producers_plot2"), full_screen=True)),
            title= "Boxplot showing the number of metabolic network producers for selected metabolites."),

        ui.card(ui.output_data_frame("cpd_exp_stat_dataframe"),full_screen=True),

    ),

    min_height="1500px"
    ),

    full_screen=True,)

    return welcome_card, compounds_exploration_card


@module.server
def cpd_tab_server(input, output, session, Data: DataStorage):

    # Helper to strip trailing counts added to labels, e.g., "Alkaloids (10/20)"
    def _strip_label(label: str) -> str:
        if label is None:
            return ""
        # Remove visual prefixes like '-- '
        while label.startswith('-') or label.startswith(' '):
            label = label.lstrip('-').lstrip()
        if "(" in label:
            return label.rsplit(" (", 1)[0]
        # Handle labels with trailing counts without parentheses, e.g., "Others 3/10"
        parts = label.split()
        if len(parts) > 1 and "/" in parts[-1]:
            return " ".join(parts[:-1])
        return label

    # Helper to get ALL available categories with optional parent filter
    def _get_all_categories_with_counts(parent_key: str | None = None):
        """Get all categories at the next level, or all categories if parent_key is None."""
        tree = Data.get_cpd_category_tree()
        
        # Navigate to the parent subtree if specified
        if parent_key:
            subtree_list = []
            Data.get_sub_tree_recursive(tree, parent_key, subtree_list)
            if not subtree_list:
                return []
            subtree = subtree_list[0]
        else:
            subtree = tree
            # If the subtree is wrapped in a single artificial root, unwrap it
            if len(subtree) == 1 and all(isinstance(v, dict) for v in subtree.values()):
                (_, subtree) = next(iter(subtree.items()))

        data_cpd_list = Data.get_compound_list(without_compartment=True)
        options = []
        for key, child in subtree.items():
            if not child:
                continue  # skip leaves (compounds)
            cpd_in_category = []
            Data.get_compounds_from_category(child, cpd_in_category)
            final_cpd_list = [cpd for cpd in data_cpd_list if cpd in cpd_in_category]
            options.append(f"{key} ({len(final_cpd_list)}/{len(cpd_in_category)})")
        return options

    def _get_children_with_counts_within(subtree: dict | None):
        """Get direct children with counts within a specific subtree (no cross-branch lookups)."""
        if subtree is None:
            return []
        data_cpd_list = Data.get_compound_list(without_compartment=True)
        options = []
        for key, child in subtree.items():
            if not child:
                continue
            cpd_in_category = []
            Data.get_compounds_from_category(child, cpd_in_category)
            final_cpd_list = [cpd for cpd in data_cpd_list if cpd in cpd_in_category]
            options.append(f"{key} ({len(final_cpd_list)}/{len(cpd_in_category)})")
        return options

    # Helpers for All/Others counts on a subtree
    def _subtree_for_key(key: str | None):
        tree = Data.get_cpd_category_tree()
        if key is None or key == "":
            return tree
        lst = []
        Data.get_sub_tree_recursive(tree, key, lst)
        return lst[0] if lst else None

    def _subtree_for_key_within(root: dict | None, key: str | None):
        """Search for a key within a given subtree only (avoids cross-branch collisions)."""
        if root is None or key is None:
            return None
        if key in root:
            return root.get(key)
        for _ck, cv in root.items():
            if isinstance(cv, dict):
                found = _subtree_for_key_within(cv, key)
                if found is not None:
                    return found
        return None

    def _compounds_recursive(subtree: dict):
        cpds = []
        if subtree is None:
            return cpds
        Data.get_compounds_from_category(subtree, cpds)
        return list(set(cpds))

    def _compounds_direct(subtree: dict):
        cpds = []
        if subtree is None:
            return cpds
        for key, child in subtree.items():
            if not child:
                cpds.append(key)
        return list(set(cpds))

    def _format_counts_label(name: str, cpds: list[str]):
        unique_cpds = list(set(cpds))
        data_cpd_list = set(Data.get_compound_list(True))
        in_data = len([c for c in unique_cpds if c in data_cpd_list])
        total = len(unique_cpds)
        return f"{name} ({in_data}/{total})"

    def _direct_child_keys(subtree: dict):
        return [k for k, v in subtree.items() if isinstance(v, dict) and v]

    def _aggregate_grandchildren_with_counts(first_key: str):
        first_sub = _subtree_for_key(first_key)
        if first_sub is None:
            return []
        # Collect unique grandchildren keys
        grandchildren = set()
        for child_k, child_v in first_sub.items():
            if isinstance(child_v, dict) and child_v:
                for gc_k, gc_v in child_v.items():
                    if isinstance(gc_v, dict) and gc_v:
                        grandchildren.add(gc_k)
        # Build labels with counts
        labels = []
        for gc in sorted(grandchildren):
            sub = _subtree_for_key(gc)
            cpds = _compounds_recursive(sub)
            labels.append(_format_counts_label(gc, cpds))
        return labels
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
    def save_raw_data_logs():
        return f"Data will be saved in {Data.raw_data_path}."

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
        Data.keep_working_dataframe("cscope_heatmap", plot_object, True)
        return plot_object

    @render.plot
    def heatmap_iscope():
        plot_object = cpd_plot_generation.result()[3][1]
        Data.keep_working_dataframe("iscope_heatmap", plot_object, True)
        return plot_object

    @render.plot
    def heatmap_added_value():
        plot_object = cpd_plot_generation.result()[3][2]
        Data.keep_working_dataframe("advalue_heatmap", plot_object, True)
        return plot_object

    @render.data_frame
    def cpd_exp_stat_dataframe():
        return cpd_plot_generation.result()[2]

    @render_widget
    def sample_percentage_production_cscope():

        return cpd_plot_generation.result()[1][0]

    @render_widget
    def sample_percentage_production_iscope():
        return cpd_plot_generation.result()[1][1]

    @render_widget
    def cpd_exp_producers_plot():
        return cpd_plot_generation.result()[0][0]

    @render_widget
    def cpd_exp_producers_plot2():
        return cpd_plot_generation.result()[0][1]

    @ui.bind_task_button(button_id="run_plot_generation")
    @reactive.extended_task
    async def cpd_plot_generation(selected_compounds, user_input1, user_color_input, sample_filter_mode, sample_filter_value, with_statistic, with_multiple_correction, multiple_correction_method, row_cluster, col_cluster, render_cpd_abundance, save_raw_data):

        if len(selected_compounds) == 0:
            return

        cpd_filtered_list = []
        for cpd in Data.get_compound_list():
            if cpd[:-3] in selected_compounds:
                cpd_filtered_list.append(cpd)

        if len(selected_compounds) == 1 and col_cluster == True:
            # Columns clustering on one compound(column) will throw an error. EmptyMatrix
            col_cluster = False

        try:
            nb_producers_boxplot = sm.render_reactive_metabolites_production_plot(Data, cpd_filtered_list, user_input1, user_color_input, sample_filter_mode, sample_filter_value, render_cpd_abundance, save_raw_data) ###
        except Exception as e:
            print(e)
            nb_producers_boxplot = [None, None]

        if user_input1 != "None":
            percent_barplot = sm.percentage_smpl_producing_cpd(Data, cpd_filtered_list, user_input1, sample_filter_mode, sample_filter_value, save_raw_data)

        else:
            ui.notification_show(
            "Sample Percentage production plot needs a metadata filter input.",
            type="warning",
            duration=6,)
            percent_barplot = [None, None]

        if with_statistic:
            try:
                stat_dataframe = sm.metabolites_production_statistical_dataframe(Data, cpd_filtered_list, user_input1, "None", with_multiple_correction, multiple_correction_method, save_raw_data)
            except Exception as e:
                print(e)
                stat_dataframe = None
        else:
            stat_dataframe = None

        cscope_heatmap, iscope_heatmap, added_value_heatmap = sm.sns_clustermap(Data, cpd_filtered_list, user_input1, row_cluster, col_cluster, sample_filter_mode, sample_filter_value, save_raw_data) ###

        return nb_producers_boxplot, percent_barplot, stat_dataframe, (cscope_heatmap, iscope_heatmap, added_value_heatmap)


    @reactive.effect
    @reactive.event(input.run_plot_generation, ignore_none=True)
    def handle_click_cpd_exploration():
        cpd_plot_generation(input.compounds_choice_input(), input.metadata_filter_input(),
                            input.color_filter_input(),
                            input.sample_filter_choice_input(), input.sample_filter_selection_input(),
                            input.exp_cpd_generate_stat_dataframe(), input.multiple_correction(),
                            input.multiple_correction_method(), input.row_cluster(), input.col_cluster(),
                            input.render_cpd_abundance(), input.save_raw_data())


    @reactive.effect
    def _update_sub_category_choices():

        if not Data.USE_METACYC_PADMET:

            return

        first_label = input.first_category_input()
        first_key = _strip_label(first_label)

        if first_key in ("None", ""):
            # No first category: clear second/third
            ui.update_selectize("category_choice_input", choices=["None"], selected="None")
            return ui.update_selectize("category_choice_input_lvl3", choices=["None"], selected="None")

        if first_key == "Others":
            outsiders = Data.get_outsider_cpd()[0]
            all_label = _format_counts_label("All", outsiders)
            ui.update_selectize("category_choice_input", choices=[all_label, "None"], selected="None")
            return ui.update_selectize("category_choice_input_lvl3", choices=["None"], selected="None")

        if first_key == "All":
            all_cpds_global = Data.get_compound_list(without_compartment=True)
            all_label = _format_counts_label("All", all_cpds_global)
            ui.update_selectize("category_choice_input", choices=["None"], selected="None")
            return ui.update_selectize("category_choice_input_lvl3", choices=["None"], selected="None")

        # Get child categories of the selected first category and add All/Others
        first_sub = _subtree_for_key(first_key)
        # All compounds under first
        all_cpds_first = _compounds_recursive(first_sub)
        all_label = _format_counts_label("All", all_cpds_first)
        # Others: direct compounds under first
        others_cpds_first = _compounds_direct(first_sub)
        others_label = _format_counts_label("Others", others_cpds_first)

        second_level_choices = _get_all_categories_with_counts(first_key)
        final_choices = [all_label, others_label] + second_level_choices

        ui.update_selectize("category_choice_input", choices=final_choices, selected=all_label)
        return ui.update_selectize("category_choice_input_lvl3", choices=[], selected=None)

    @reactive.effect
    def _update_third_category_choices():

        if not Data.USE_METACYC_PADMET:

            return

        second_label = input.category_choice_input()
        second_key = _strip_label(second_label)

        # If none selected, clear third-level choices
        if second_key in ("None", "", None):
            return ui.update_selectize("category_choice_input_lvl3", choices=["None"], selected="None")

        first_key = _strip_label(input.first_category_input())

        # If first is Others or missing, no third level applies
        if first_key in ("Others", "", None):
            return ui.update_selectize("category_choice_input_lvl3", choices=["None"], selected="None")

        # If second is Others, no third level applies
        if second_key == "Others":
            return ui.update_selectize("category_choice_input_lvl3", choices=["None"], selected="None")

        # Build choices depending on second selection
        if second_key == "All":
            # Aggregate direct children across all second-level categories
            third_level_choices = _aggregate_grandchildren_with_counts(first_key)
            # All/Others at third relative to first (aggregate)
            first_sub = _subtree_for_key(first_key)
            if first_sub is None:
                return ui.update_selectize("category_choice_input_lvl3", choices=[], selected=None)
            all_cpds_first = _compounds_recursive(first_sub)
            others_cpds_agg = []
            for child_k, child_v in first_sub.items():
                if isinstance(child_v, dict) and child_v:
                    others_cpds_agg.extend(_compounds_direct(child_v))
            others_cpds_agg = list(set(others_cpds_agg))
            all_label = _format_counts_label("All", all_cpds_first)
            others_label = _format_counts_label("Others", others_cpds_agg)
            final_choices = [all_label, others_label] + third_level_choices
            return ui.update_selectize("category_choice_input_lvl3", choices=final_choices, selected=all_label)
        else:
            # Third-level under the specific second category
            first_sub = _subtree_for_key(first_key)
            second_sub = _subtree_for_key_within(first_sub, second_key)
            if second_sub is None:
                return ui.update_selectize("category_choice_input_lvl3", choices=[], selected=None)

            third_level_choices = _get_children_with_counts_within(second_sub)
            all_cpds_second = _compounds_recursive(second_sub)
            others_cpds_second = _compounds_direct(second_sub)
            all_label = _format_counts_label("All", all_cpds_second)
            others_label = _format_counts_label("Others", others_cpds_second)
            final_choices = [all_label, others_label] + third_level_choices
            return ui.update_selectize("category_choice_input_lvl3", choices=final_choices, selected=all_label)

    @reactive.effect
    def _update_compounds_choices():

        if not Data.USE_METACYC_PADMET:

            return

        # Prefer the deepest selection available: third level > second level > first level
        third_label = input.category_choice_input_lvl3()
        second_label = input.category_choice_input()
        first_label = input.first_category_input()

        # Prefer deepest non-None/non-"None"; otherwise fall back to first
        target_label = None
        if third_label not in (None, "None"):
            target_label = third_label
        elif second_label not in (None, "None"):
            target_label = second_label
        else:
            target_label = first_label

        target_key = _strip_label(target_label)

        if target_key in ("", "None", None):
            return ui.update_selectize("compounds_choice_input", choices=[])

        # Determine compounds based on All/Others logic
        first_key = _strip_label(first_label)
        second_key = _strip_label(second_label)
        third_key = _strip_label(third_label) if third_label else None

        data_cpd_list = set(Data.get_compound_list(True))

        def _cpds_in_data(cpds):
            return sorted(list(set([c for c in cpds if c in data_cpd_list])))

        # Special handling when first category is None/empty
        if first_key in ("None", ""):
            return ui.update_selectize("compounds_choice_input", choices=[])

        # Special handling when first category is All
        if first_key == "All":
            final_cpd_list = _cpds_in_data(Data.get_compound_list(without_compartment=True))
            return ui.update_selectize("compounds_choice_input", choices=final_cpd_list, selected=final_cpd_list)

        # Special handling when first category is Others
        if first_key == "Others":
            outsiders = Data.get_outsider_cpd()[0]
            final_cpd_list = _cpds_in_data(outsiders)
            return ui.update_selectize("compounds_choice_input", choices=final_cpd_list, selected=final_cpd_list)

        # If target is Others but subtree missing, return empty
        if target_key == "Others" and _subtree_for_key(target_key) is None:
            return ui.update_selectize("compounds_choice_input", choices=[])

        first_sub = _subtree_for_key(first_key)

        # Compute target subtree and compound list
        if third_key and third_key not in ("None", ""):
            # Third-level selected
            if third_key == "All":
                # Compounds under second selection (or first if second is All)
                if second_key == "All":
                    sub = first_sub
                else:
                    sub = _subtree_for_key_within(first_sub, second_key)
                cpds_found = _compounds_recursive(sub)
            elif third_key == "Others":
                # Direct compounds under second level (or aggregated across all seconds if second is All)
                if second_key == "All":
                    cpds_found = []
                    if first_sub:
                        for ck, cv in first_sub.items():
                            if isinstance(cv, dict) and cv:
                                cpds_found.extend(_compounds_direct(cv))
                    cpds_found = list(set(cpds_found))
                else:
                    sub = _subtree_for_key_within(first_sub, second_key)
                    cpds_found = _compounds_direct(sub)
            else:
                # Specific third-level node within the selected second branch (or all children when second is All)
                if second_key == "All":
                    sub = _subtree_for_key_within(first_sub, third_key)
                else:
                    second_sub = _subtree_for_key_within(first_sub, second_key)
                    sub = _subtree_for_key_within(second_sub, third_key)
                cpds_found = _compounds_recursive(sub)
        elif second_key and second_key not in ("", "None"):
            # Second-level selected
            if second_key == "All":
                sub = first_sub
                cpds_found = _compounds_recursive(sub)
            elif second_key == "Others":
                sub = first_sub
                cpds_found = _compounds_direct(sub)
            else:
                sub = _subtree_for_key_within(first_sub, second_key)
                cpds_found = _compounds_recursive(sub)
        else:
            # Only first-level selected
            sub = _subtree_for_key(first_key)
            cpds_found = _compounds_recursive(sub)

        final_cpd_list = _cpds_in_data(cpds_found)

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
