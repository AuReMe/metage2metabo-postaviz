import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from shiny import App
from shiny import render
from shiny import run_app
from shiny import ui
from shiny import reactive
from shinywidgets import output_widget
from shinywidgets import render_widget
import warnings
import pandas as pd
from ontosunburst.ontosunburst import ontosunburst as ontos
import time

import m2m_postaviz.data_utils as du
from m2m_postaviz.data_struct import DataStorage

# from profilehooks import profile
# from pympler.asizeof import asizeof

def run_shiny(data: DataStorage):
    
    warnings.filterwarnings("ignore", category=FutureWarning, module="plotly.express")

    list_of_cpd = data.get_compound_list()

    factor_list = data.get_factors()
    factor_list.insert(0, "None")

    metadata_label = data.get_factors()
    metadata_label.remove("smplID")

    taxonomic_rank = data.get_taxonomy_rank()


    list_of_bins = data.get_bins_list()
    if data.HAS_TAXONOMIC_DATA:
        list_of_bins = du.associate_bin_taxonomy(list_of_bins, data.get_taxonomic_dataframe())


    bins_count = data.get_bins_count()

    all_dataframe = {"global_production_test_dataframe": None, "global_production_plot_dataframe": None, "metabolites_production_test_dataframe": None, "metabolites_production_plot_dataframe": None}

    current_bin_dataframe = None

    ### ALL CARD OBJECT TO BE ARRANGED ###

    producer_plot =   ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("box_inputx1", "Label for X axis", factor_list),
                    ui.input_select("box_inputx2", "Label for 2nd X axis", factor_list),
                    ui.input_selectize("box_inputy1", "Compound for Y axis", list_of_cpd, multiple=True, selected=list_of_cpd[0]),

                    ui.input_checkbox("ab_norm", "With normalised data"),
                    ui.input_checkbox("multiple_correction_metabo_plot", "Multiple test correction"),
                    ui.panel_conditional("input.multiple_correction_metabo_plot",ui.input_select("multiple_test_method_metabo","Method",
                                                                                                ["bonferroni","sidak","holm-sidak","holm","simes-hochberg","hommel","fdr_bh","fdr_by","fdr_tsbh","fdr_tsbky"],
                                                                                                selected="bonferroni",)),


                    ui.input_action_button("export_metabolites_plot_button", "Export plot dataframe"),
                    ui.output_text_verbatim("export_metabolites_plot_dataframe_txt", True),
                    ui.input_action_button("export_metabolites_test_button", "Export stats dataframe"),
                    ui.output_text_verbatim("export_metabolites_test_dataframe_txt", True),
                    width=350,
                    gap=30,

            ),
            ui.layout_column_wrap(

                ui.card(output_widget("producer_plot"),full_screen=True),
         
                ui.card(ui.output_data_frame("producer_test_dataframe"),full_screen=True)
                
            ),
        ),
        full_screen=True
        )

    if data.HAS_TAXONOMIC_DATA:
        taxonomy_boxplot = ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("tax_inpx1", "Label for X axis", factor_list),
                    ui.input_select("tax_inpx2", "Label for 2nd X axis", factor_list),
                    # ui.input_selectize("tax_inpy1", "Taxa for Y axis", long_taxo_df["Taxa"].unique().tolist(), multiple=True),
                    ui.input_checkbox("taxo_norm", "With normalised data"),
                    width=350,
                ),
                output_widget("taxonomic_boxplot"),
            ),
            full_screen=True
        )

    else:
        taxonomy_boxplot = ui.output_text_verbatim("no_taxonomy", True),

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
        ui.row(
                ui.input_select("metadata_factor", "Current column: ", metadata_label, selected=metadata_label[0]),
                ui.input_select("metadata_dtype", "dtype: ", ["category", "str", "int", "float"]),
                ui.input_action_button("dtype_change", "Update")
               ),
        ui.output_text_verbatim("update_metadata_log", True),
        ui.output_data_frame("metadata_table")
        )

    pcoa_card = ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select(id="pcoa_color", label="Plot color.", choices=metadata_label, selected=metadata_label[0]),
                ui.output_ui("pcoa_ui"),
                ui.help_text(ui.output_text("display_warning_pcoa")),
                width=300,
                gap=30,
            ),
        output_widget("pcoa_plot")
        ),
        full_screen=True
    )

    custom_pcoa_card = ui.card(ui.card_header("Make a custom pcoa filtered by columns values."),
        ui.layout_sidebar(
            ui.sidebar(

                ui.input_checkbox("pcoa_custom_abundance_check", "Use abundance data."),
                ui.input_select(id="custom_pcoa_selection", label="Choose column", choices=metadata_label, selected=metadata_label[0]),
                ui.output_ui("pcoa_custom_choice"),
                ui.input_select(id="pcoa_custom_color", label="Color.", choices=metadata_label, selected=metadata_label[0]),


                ui.input_task_button("run_custom_pcoa","Go"),
                width=300,
                gap=35,
            ),
        output_widget("pcoa_custom_plot"),
        ),
        full_screen=True
    )

    bins_exploration_card = ui.card(
        ui.card_header("Bins exploration"),
        ui.card_body(
            ui.layout_sidebar(
                ui.sidebar(

                    ui.input_selectize("rank_choice", "Choose a taxonomic rank.", taxonomic_rank, selected=taxonomic_rank[0], multiple=False, width='400px'),
                    ui.output_ui("rank_unique_choice"),

                    ui.input_selectize("bin_choice", "Choose", list_of_bins, selected=list_of_bins[0], multiple=False, width='400px'),
                    ui.output_text("bin_size_text"),
                    ui.input_selectize("bin_factor", "Choose", factor_list, selected=factor_list[0], multiple=False, width='400px'),
                    ui.input_task_button("run_bin_exploration","Go"),
                    ui.output_text("iscope_info"),
                    
                    ui.output_ui("bin_sample_select"),
                width=350,
                gap=35,
                bg='lightgrey'
            ),
            ui.card(output_widget("bin_count_plot"),full_screen=True),
            ui.card(output_widget("bin_abundance_plot"),full_screen=True),
            ui.card(output_widget("ontosunburst"),full_screen=True),
        )
    ),full_screen=True)
    
    ### APPLICATION TREE ###

    app_ui = ui.page_fillable(
        ui.navset_tab(
            # ui.nav_panel("Exploration", total_production_plot, producer_plot, taxonomy_boxplot),
            ui.nav_panel(
                "Metadata",
                metadata_table,
                # dev_table
            ),
            ui.nav_panel(
                "PCOA",
                    pcoa_card,
                    custom_pcoa_card
            ),
            ui.nav_panel("Bins",
                         bins_exploration_card
            )
            ),
        )

    def server(input, output, session):


        @render.ui
        def rank_unique_choice():

            rank_choice = input.rank_choice()
            
            df = data.get_taxonomic_dataframe()

            choices = df[rank_choice].unique().tolist()

            return ui.TagList(
                ui.input_selectize("rank_unique_choice", "Select", choices=choices, multiple=False,)
                )

        @render.ui
        def bin_sample_select():
            df = run_exploration.result()[2]
            choice = df["smplID"].tolist()

            return ui.TagList(
                ui.input_selectize("bin_sample_select_input", "Choose sample", choices=choice, multiple=False,)
                )

        @render.text
        def bin_size_text():
            bin_input = input.bin_choice().split("/")[0]

            return "Found in "+str(bins_count[bin_input])+" sample(s)"

        @render.text
        def iscope_info():
            # bin_input = input.bin_choice().split("/")[0]
            # iscope_prod = data.get_iscope_production(bin_input)
            # df = run_exploration.result()[2]
            # cscope_prod = df.loc[df["Count"] == df["Count"].max()]["Production"].tolist()

            # res = []
            # for cpd in cscope_prod:
            #     if cpd not in iscope_prod:
            #         res.append(cpd)
            timer = run_exploration.result()[3]
            return f"Took {timer} seconds to run."

        @render_widget
        def ontosunburst():

            smpl_choice = input.bin_sample_select_input()
            bin_input = input.bin_choice().split("/")[0]

            df = run_exploration.result()[2]
            
            if not du.is_indexed_by_id(df):
                df.set_index("smplID", inplace=True)

            cscope_prod = df.at[smpl_choice,"Production"]
            iscope_prod = data.get_iscope_production(bin_input)

            res = []
            reference_set = []
            for cpd in cscope_prod:
                if cpd not in iscope_prod:
                    res.append(cpd)
                reference_set.append(cpd)

            for i, cpd in enumerate(res):
                res[i] = cpd[:-3]

            for i, cpd in enumerate(cscope_prod):
                reference_set[i] = cpd[:-3]

            interest_set = pd.Series( (cpd for cpd in res) )


            fig = ontos(interest_set=interest_set, reference_set=reference_set, ontology='metacyc')
            return fig

        @render_widget
        def bin_abundance_plot():
            return run_exploration.result()[1]

        @render_widget
        def bin_count_plot():
            return run_exploration.result()[0]

        @ui.bind_task_button(button_id="run_custom_pcoa")
        @reactive.extended_task
        async def run_exploration(bin_choice, factor, rank, rank_choice):
            start_timer = time.time()
            list_of_bin_in_rank = data.get_bin_list_from_taxonomic_rank(rank, rank_choice)

            #### Taxonomic dataframe can contain MORE information and MORE bin than the data who can be un subset of the whole data. Filtering is needed.
            set_bin = set(list_of_bins)
            [x for x in list_of_bin_in_rank if x in set_bin]

            res = []
            for bin in [x for x in list_of_bin_in_rank if x in set_bin]:
                df = data.from_bin_get_dataframe(bin, factor)
                df.insert(0 , "binID", bin)
                res.append(df)
            
            res = pd.concat(res)
            res.fillna(0)
            print(res)
            if factor == "None":
                res.sort_index(inplace=True)
            else:
                res.sort_values(by=factor,inplace=True)

            max_count_range = res["Count"].max() + res["Count"].max() * 0.01
            min_count_range = res["Count"].min() - res["Count"].min() * 0.02

            fig1 = px.bar(res, x="smplID", y="Count", text="Count", color="smplID" if factor =="None" else factor)

            fig2 = px.bar(res, x="smplID", y="Abundance", color="Abundance")


            # # Create figure with secondary y-axis
            # fig = make_subplots(specs=[[{"secondary_y": True}]])

            # # Add traces
            # fig.add_trace(
            #     go.Bar(x=res.index, y=res["Count"], name="count subplot", text=res["Count"], offsetgroup=1),
            #     secondary_y=False
            # )

            # fig.add_trace(
            #     go.Bar(x=res.index, y=res["Abundance"], name="abundance subplot", text=res["Abundance"], offsetgroup=2),
            #     secondary_y=True
            # )

            # # Add figure title
            # fig.update_layout(
            #     title_text=f"Number of compounds {bin_choice} produce by sample(s) with their respective abundance.",
            #     barmode="group"
            # )

            # # Set x-axis title
            # fig.update_xaxes(title_text="Sample ID")

            # # Set y-axes titles
            # fig.update_yaxes(title_text="<b>Compounds produced</b>", secondary_y=False)
            # fig.update_yaxes(title_text="<b>Abundance in %</b>", secondary_y=True)

            return fig1, fig2, res, time.time() - start_timer

        @reactive.effect
        @reactive.event(input.run_bin_exploration, ignore_none=True)
        def handle_click_bin_exploration():
            run_exploration(input.bin_choice().split("/")[0], input.bin_factor(), input.rank_choice(), input.rank_unique_choice())

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
            
            if abundance:
                df = data.get_normalised_abundance_dataframe()
            else:
                df = data.get_main_dataframe()

            metadata = data.get_metadata()

            if du.is_indexed_by_id(df):
                df.reset_index(inplace=True)

            if du.is_indexed_by_id(metadata):
                metadata.reset_index(inplace=True)

            if du.serie_is_float(metadata[column]):

                selected_sample = metadata.loc[(metadata[column] >= choices[0]) & (metadata[column] <= choices[1])]["smplID"]
                df = df.loc[df["smplID"].isin(selected_sample)]
                metadata = metadata.loc[metadata["smplID"].isin(selected_sample)]

            else:

                selected_sample = metadata.loc[metadata[column].isin(choices)]["smplID"]
                df = df.loc[df["smplID"].isin(selected_sample)]
                metadata = metadata.loc[metadata["smplID"].isin(selected_sample)]

            plot_df = du.pcoa_alternative_method(df, metadata)

            fig = px.scatter(plot_df, x="PC1", y="PC2",
                             color= color  
                             )

            return fig

        @render_widget
        def pcoa_custom_plot():
            return make_custom_pcoa.result()

        @render.ui
        @reactive.event(input.custom_pcoa_selection)
        def pcoa_custom_choice():

            df = data.get_metadata()
            value = df[input.custom_pcoa_selection()]

            if not du.serie_is_float(value):

                return ui.TagList(
                            ui.input_selectize("pcoa_custom_choice", "Get only in column:", value.unique().tolist(),
                                               selected=value.unique().tolist(),
                                               multiple=True,
                                               options={"plugins": ['clear_button']}),)
            else:

                return ui.TagList(
                            ui.input_slider("pcoa_custom_choice", "Set min/max filter:", min=value.min(), max=value.max(), value=[value.min(),value.max()]),)
        
        @render_widget
        def pcoa_plot():
            # Get all parameters.
            selected_col = input.pcoa_color()

            df = data.get_pcoa_dataframe()

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

            df = data.get_pcoa_dataframe()
            value = df[input.pcoa_color()]

            if not du.serie_is_float(value):

                return ui.TagList(
                            ui.input_selectize("pcoa_selectize", "Show only:", value.unique().tolist(),
                                               selected=value.unique().tolist(),
                                               multiple=True,
                                               options={"plugins": ['clear_button']}),)
            else:

                return ui.TagList(
                            ui.input_slider("pcoa_slider", "Set min/max filter:", min=value.min(), max=value.max(), value=[value.min(),value.max()]),)
        
        @render_widget
        def taxonomic_boxplot():

            if not data.HAS_TAXONOMIC_DATA:
                return
            
            df = data.get_taxonomic_dataframe()
            x1, x2, y1 = input.tax_inpx1(), input.tax_inpx2(), input.tax_inpy1()

            if len(y1) == 0:
                y1 = df["Taxa"].unique()

            if x1 == "None":
                return px.box(df, y="Nb_taxon")
            
            if x2 == "None":
                fig = px.box(
                    df,
                    x=x1,
                    y="Nb_taxon",
                    color=x1,
                )

                return fig
            
            else:

                conditionx2 = df[x2].unique()

                fig = go.Figure()

                for condition in conditionx2:
                    fig.add_trace(
                        go.Box(
                            x=df.loc[df[x2] == condition][x1],
                            y=df["Nb_taxon"],
                            name=str(condition),
                        )
                    )

                fig.update_layout(boxmode="group")

                return fig

        @render.data_frame
        def producer_test_dataframe():

            y1, x1, x2 = list(input.box_inputy1()), input.box_inputx1(), input.box_inputx2()

            if len(y1) == 0:
                return
            
            if x1 == "None":
                return 

            multipletests_correction = input.multiple_correction_metabo_plot()

            if multipletests_correction:
                multipletests_method = input.multiple_test_method_metabo()
            else:
                multipletests_method = "hs"

            if x2 == "None":

                df = data.get_metabolite_production_dataframe()[[*y1,x1]]
                df = df.dropna()

                if du.serie_is_float(df[x1]):

                    correlation_results = []

                    for y_value in y1:

                        correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), x1))

                    return pd.concat(correlation_results)

               
                res = du.preprocessing_for_statistical_tests(df, y1, x1, multipletests = multipletests_correction, multipletests_method = multipletests_method)
                all_dataframe["metabolites_production_test_dataframe"] = res

                return res
            
            if x1 != x2:

                df = data.get_metabolite_production_dataframe()[[*y1,x1,x2]]
                df = df.dropna()

                if du.serie_is_float(df[x1]): # First input is Float type

                    if du.serie_is_float(df[x2]): # Second input is Float type

                        correlation_results = []

                        for y_value in y1:

                            correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), str(x1+" "+y_value)))
                            correlation_results.append(du.correlation_test(df[y_value].to_numpy(), df[x1].to_numpy(), str(x2+" "+y_value)))

                        return pd.concat(correlation_results)

                    else: # Second input is not Float type

                        correlation_results = []

                        for y_value in y1:

                            for x2_unique_value in df[x2].unique():
                          
                                factor_array = df.loc[df[x2] == x2_unique_value][x1]
                                value_array = df.loc[df[x2] == x2_unique_value][y_value]

                                correlation_results.append(du.correlation_test(value_array.to_numpy(), factor_array.to_numpy(), str(x2_unique_value)+" "+y_value))

                        return pd.concat(correlation_results)
                    
                else:   

                    res = du.preprocessing_for_statistical_tests(df, y1, x1, x2, multipletests = multipletests_correction, multipletests_method = multipletests_method)
                    all_dataframe["metabolites_production_test_dataframe"] = res

                return res
                
        @render_widget
        def producer_plot():

            compound_input, first_input, second_input = list(input.box_inputy1()), input.box_inputx1(), input.box_inputx2()

            if len(compound_input) == 0:
                return

            producer_data = data.get_metabolite_production_dataframe()
            producer_data = producer_data.set_index("smplID")

            if first_input == "None":
                df = producer_data[[*compound_input]]
                all_dataframe["metabolites_production_plot_dataframe"] = df
                return px.box(df,y=compound_input)

            if second_input == "None" or first_input == second_input:
                df = producer_data[[*compound_input,first_input]]
                df = df.dropna()
                all_dataframe["metabolites_production_plot_dataframe"] = df
                has_unique_value = du.has_only_unique_value(df, first_input)


                return px.bar(df,x=first_input,y=compound_input,color=first_input) if has_unique_value else px.box(df,x=first_input,y=compound_input,color=first_input)

            df = producer_data[[*compound_input,first_input,second_input]]
            df = df.dropna()
            all_dataframe["metabolites_production_plot_dataframe"] = df
            has_unique_value = du.has_only_unique_value(df, first_input, second_input)

            return px.bar(df, x=first_input, y=compound_input,color=second_input) if has_unique_value else px.box(df, x=first_input, y=compound_input,color=second_input, boxmode="group")

        @render.data_frame
        def production_test_dataframe():

            x1, x2 = input.prod_inputx1(), input.prod_inputx2()
            
            # No input selected
            if x1 == "None":
                return 
            
            multipletests_correction = input.multiple_correction_global_plot()

            if multipletests_correction:
                multipletests_method = input.multiple_test_method_global()
            else:
                multipletests_method = "hs"

            with_normalised_data = input.prod_norm()

            if with_normalised_data and data.HAS_ABUNDANCE_DATA:
                column_value = "Total_abundance_weighted"
            else:
                column_value = "Total_production"

            df = data.get_global_production_dataframe()
            
            # At least first axis selected
            if x2 == "None":
                df = df[[column_value,x1]]
                df = df.dropna()

                if du.serie_is_float(df[x1]):

                    res = du.correlation_test(df[column_value].to_numpy(), df[x1].to_numpy(), x1)

                    return res
                    

                res = du.preprocessing_for_statistical_tests(df, [column_value], x1, multipletests=multipletests_correction, multipletests_method=multipletests_method)
                all_dataframe["global_production_test_dataframe"] = res

                return res
            
            # Both axis have been selected and are not equal.
            if x1 != x2:

                df = df[[column_value,x1,x2]]
                df = df.dropna()

                if du.serie_is_float(df[x1]):

                    if du.serie_is_float(df[x2]):

                        # Double cor
                        res1 = du.correlation_test(df[column_value].to_numpy(), df[x1].to_numpy(), x1)
                        res2 = du.correlation_test(df[column_value].to_numpy(), df[x2].to_numpy(), x2)
                        return pd.concat([res1,res2])
                    
                    else:

                        # cor filtered by second categorical factor .loc
                        all_results = []
                        for unique_x2_value in df[x2].unique():

                            value_array = df.loc[df[x2] == unique_x2_value][column_value]
                            factor_array = df.loc[df[x2] == unique_x2_value][x1]

                            all_results.append(du.correlation_test(value_array, factor_array, unique_x2_value))

                        return pd.concat(all_results)

                res = du.preprocessing_for_statistical_tests(df, [column_value], x1, x2, multipletests=multipletests_correction, multipletests_method=multipletests_method)
                all_dataframe["global_production_test_dataframe"] = res

                return res
            
            return
        
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

            with_normalised_data = input.prod_norm()

            if with_normalised_data and data.HAS_ABUNDANCE_DATA:
                column_value = "Total_abundance_weighted"
            else:
                column_value = "Total_production"

            df = data.get_global_production_dataframe()
            

            inputx1 , inputx2 = input.prod_inputx1(), input.prod_inputx2()
        
            if inputx1 == "None":

                all_dataframe["global_production_plot_dataframe"] = df
                return px.box(df, y=column_value, title=f"Total production.")
            
            elif inputx2 == "None" or inputx1 == inputx2:

                df = df[[column_value,inputx1]]
                drop_count = len(df)
                df = df.dropna()

                if drop_count - len(df) > 0:
                    ui.notification_show(f'{drop_count - len(df)} lines dropped(NaN values).',duration=8,close_button=True,type="warning")

                all_dataframe["global_production_plot_dataframe"] = df

                if du.has_only_unique_value(df, inputx1):

                    return px.bar(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
                
                else:
                    
                    fig = px.box(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
                    return fig
                
            else:

                df = df[[column_value,inputx1,inputx2]]
                df = df.dropna()
                all_dataframe["global_production_plot_dataframe"] = df
                has_unique_value = du.has_only_unique_value(df, inputx1, inputx2)

                return px.bar(df,x=inputx1,y=column_value,color=inputx2) if has_unique_value else px.box(df,x=inputx1,y=column_value,color=inputx2)

        @render.data_frame
        def metadata_table():
            df = data.get_metadata()
            return df
        
        @render.text()
        @reactive.event(input.dtype_change)
        def update_metadata_log():
            return
            text = "No changes applied."

            factor_choice, dtype_choice = input.metadata_factor(), input.metadata_dtype()

            df = data.get_metadata()

            df_prod = data.get_metabolite_production_dataframe

            df_tot = data.get_global_production_dataframe()

            df_pcoa = data.get_pcoa_dataframe()

            try:
                df[factor_choice] = df[factor_choice].astype(dtype_choice)
                df_prod[factor_choice] = df_prod[factor_choice].astype(dtype_choice)
                df_tot[factor_choice] = df_tot[factor_choice].astype(dtype_choice)
                df_pcoa[factor_choice] = df_pcoa[factor_choice].astype(dtype_choice)
                text = f"Column {factor_choice} changed to {dtype_choice}."
                data.set_main_metadata(df)

            except ValueError as e:
                text = f"Cannot perform changes, {e}"

            return text
        
        @render.text
        def no_taxonomy():
            return "No taxonomic data provided."

    app = App(app_ui, server)
    run_app(app=app)
