import plotly.express as px
import plotly.graph_objects as go
from shiny import App
from shiny import render
from shiny import run_app
from shiny import ui
from shiny import reactive
from shinywidgets import output_widget
from shinywidgets import render_widget
import warnings

import m2m_postaviz.data_utils as du
from m2m_postaviz.data_struct import DataStorage

from profilehooks import profile
from pympler.asizeof import asizeof

def run_shiny(data: DataStorage):
    ###
    warnings.filterwarnings("ignore", category=FutureWarning, module="plotly.express")

    list_of_cpd = data.get_compound_list()

    factor_list = data.get_factors()
    factor_list.insert(0, "None")

    metadata_label = data.get_factors()
    metadata_label.remove("smplID")

    all_dataframe = {"global_production_test_dataframe": None, "global_production_plot_dataframe": None, "metabolites_production_test_dataframe": None, "metabolites_production_plot_dataframe": None}

    ### ALL CARD OBJECT TO BE ARRANGED ###

    producer_plot =   ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("box_inputx1", "Label for X axis", factor_list),
                    ui.input_select("box_inputx2", "Label for 2nd X axis", factor_list),
                    ui.input_selectize("box_inputy1", "Compound for Y axis", list_of_cpd, multiple=True, selected=list_of_cpd[0]),
                    ui.input_checkbox("ab_norm", "With normalised data"),
                    ui.input_action_button("export_metabolites_plot_button", "Export plot dataframe"),
                    ui.output_text_verbatim("export_metabolites_plot_dataframe_txt", True),
                    ui.input_action_button("export_metabolites_test_button", "Export stats dataframe"),
                    ui.output_text_verbatim("export_metabolites_test_dataframe_txt", True),
                    width=350,
                    gap=30,

            ),
        output_widget("producer_plot"), 
        ui.output_data_frame("producer_test_dataframe")
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
                ui.input_action_button("export_global_production_plot_dataframe_button", "Save plot dataframe"),
                ui.output_text_verbatim("export_global_production_plot_dataframe_txt", True),
                ui.input_action_button("export_global_production_test_button", "Export stats dataframe"),
                ui.output_text_verbatim("export_global_production_test_dataframe", True),
                width=350,
                gap=30,
            ),
            output_widget("total_production_plot"),
            ui.output_data_frame("production_test_dataframe")
        
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

    pcoa_plot_dev_table = ui.card(
        ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_select(id="pcoa_color", label="Plot color.", choices=metadata_label, selected=metadata_label[0]),
                        ui.output_ui("pcoa_ui"),
                        width=350,
                        gap=30,
                    ),
        output_widget("pcoa_plot")
        ),
        full_screen=True
    )


    ### APPLICATION TREE ###

    app_ui = ui.page_fillable(
        ui.navset_tab(
            ui.nav("Exploration", total_production_plot, producer_plot, taxonomy_boxplot),
            ui.nav(
                "Metadata",
                metadata_table,
                # dev_table
            ),
            ui.nav(
                "PCOA",
                    pcoa_plot_dev_table,
                ),
            ),
        )

    def server(input, output, session):

        @render_widget
        def pcoa_plot():

            # Get all parameters.
            selected_col = input.pcoa_color()

            df = data.get_pcoa_dataframe()

            # value = df[selected_col]

            # min_limit = input.pcoa_float_limit()[0]

            # max_limit = input.pcoa_float_limit()[1]

            # Check column dtype.
            # if value.dtype == 'float64' or value.dtype == 'int64':

            #     df = df.loc[(df[selected_col] <= max_limit) & (df[selected_col] >= min_limit)]

            fig = px.scatter(df, x="PC1", y="PC2", color=selected_col)

            return fig

        # Value.dtype does not seem to be reliable with Nan in columns.
        # @render.ui
        # @reactive.event(input.pcoa_color)
        # def pcoa_ui():
        #     value = pcoa_dataframe[input.pcoa_color()]
        #     if value.dtype == 'float64' or value.dtype == 'int64':
        #         return ui.TagList(
        #                     ui.input_slider("pcoa_float_limit", "Show only :", min=value.min(), max=value.max(), value=[value.min(),value.max()]),
        #                 )
        #     else:
        #         return 
        
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
            # Get input

            y1, x1, x2 = list(input.box_inputy1()), input.box_inputx1(), input.box_inputx2()

            if len(y1) == 0:
                return
            # No input selected
            if x1 == "None":
                return 

            # At least first axis selected
            if x2 == "None":

                df = data.get_metabolite_production_dataframe()[[*y1,x1]]
                df = df.dropna()
                if df[x1].dtype == 'float64' and len(df[x1]) > 100:

                    print(len(df[x1].unique()), "unique value")
                    df[x1] = df[x1].round(1)
                    print(len(df[x1].unique()), "unique value after round")

                tested_data = {}
                for layer_1 in df[x1].unique():
                    selected_data = df.loc[df[x1] == layer_1][y1]

                    if len(selected_data) > 1:
                        tested_data[layer_1] = {}
                        tested_data[layer_1]["data"] = selected_data
                        tested_data[layer_1]["n_data"] = len(selected_data)

                res = du.stat_on_plot(tested_data, 1)
                all_dataframe["metabolites_production_test_dataframe"] = res

                return res
            # Both axis have been selected
            else:

                df = data.get_metabolite_production_dataframe()[[*y1,x1,x2]]
                df = df.dropna()
                if df[x1].dtype == 'float64' and len(df[x1]) > 100:

                    print(len(df[x1].unique()), "unique value")
                    df[x1] = df[x1].round(1)
                    print(len(df[x1].unique()), "unique value after round")

                if df[x2].dtype == 'float64' and len(df[x2]) > 100:

                    print(len(df[x2].unique()), "unique value")
                    df[x2] = df[x2].round(1)
                    print(len(df[x2].unique()), "unique value after round")

                tested_data = {}

                for layer_1 in df[x1].unique():
                    tested_data[layer_1] = {}

                    for layer_2 in df[x2].unique():
                        selected_data = df.loc[(df[x1] == layer_1) & (df[x2] == layer_2)][y1]

                        if len(selected_data) > 1:
                            tested_data[layer_1][layer_2] = {}
                            tested_data[layer_1][layer_2]["data"] = selected_data
                            tested_data[layer_1][layer_2]["n_data"] = len(selected_data)

                res = du.stat_on_plot(tested_data, 2)
                all_dataframe["metabolites_production_test_dataframe"] = res

                return res
                
        @render_widget
        def producer_plot():

            compound_input, first_input, second_input = list(input.box_inputy1()), input.box_inputx1(), input.box_inputx2()

            producer_data = data.get_metabolite_production_dataframe()
            producer_data = producer_data.set_index("smplID")

            if len(compound_input) == 0:
                return

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

        @render.data_frame()
        def production_test_dataframe():

            x1, x2 = input.prod_inputx1(), input.prod_inputx2()

            # No input selected
            if x1 == "None":
                return 
            
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
                if df[x1].dtype == 'float64' and len(df[x1]) > 100:

                    print(df[x1].dtype)
                    print(len(df[x1].unique()), "unique value")
                    df[x1] = df[x1].round(1)
                    print(len(df[x1].unique()), "unique value after round")

                tested_data = {}

                for layer_1 in df[x1].unique():
                    selected_data = df.loc[df[x1] == layer_1][column_value]

                    if len(selected_data) > 1:
                        tested_data[layer_1] = {}
                        tested_data[layer_1]["data"] = selected_data
                        tested_data[layer_1]["n_data"] = len(selected_data)

                res = du.stat_on_plot(tested_data, 1)
                all_dataframe["global_production_test_dataframe"] = res

                return res
            
            # Both axis have been selected and are not equal.
            if x1 != x2:

                df = df[[column_value,x1,x2]]
                df = df.dropna()
                for factor_columns in df[[x1,x2]].columns:

                    # If Series type is Float64 AND lenght of uniques values > 200 -> Round all numbers to the 10th. PERFORMANCE ISSUE
                    if df[factor_columns].dtype == 'float64' and len(df[factor_columns]) > 100:

                        print(df[x2].dtype)
                        print(len(df[x2].unique()), "unique value")
                        df[x2] = df[x2].round(1)
                        print(len(df[x2].unique()), "unique value after round")

                tested_data = {}

                for layer_1 in df[x1].unique():
                    tested_data[layer_1] = {}

                    for layer_2 in df[x2].unique():
                        selected_data = df.loc[(df[x1] == layer_1) & (df[x2] == layer_2)][column_value]

                        if len(selected_data) > 1:
                            tested_data[layer_1][layer_2] = {}
                            tested_data[layer_1][layer_2]["data"] = selected_data.values
                            tested_data[layer_1][layer_2]["n_data"] = len(selected_data)

                # In case tested_data has no value (barplot)
                if all(len(tested_data[k]) == 0 for k in tested_data.keys()):
                    print("No values in each keys")
                    return
                
                # In case tested_data is not a double layer dict
                if all(len(tested_data[k]) == 1 for k in tested_data.keys()):
                    print("Only one pair per layer.")
                    return
                
                else:
                    res = du.stat_on_plot(tested_data, 2)
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
        def export_metabolites_test_dataframe():

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
                df = df.dropna()
                all_dataframe["global_production_plot_dataframe"] = df

                if du.has_only_unique_value(df, inputx1):
                    return px.bar(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
                else:
                    return px.box(df, x=inputx1 , y=column_value, color=inputx1, title=f"Total compound production filtered by {inputx1}")
                
            else:

                df = df[[column_value,inputx1,inputx2]]
                df = df.dropna()
                all_dataframe["global_production_plot_dataframe"] = df
                has_unique_value = du.has_only_unique_value(df, inputx1, inputx2)

                return px.bar(df,x=inputx2,y=column_value,color=inputx1) if has_unique_value else px.box(df,x=inputx2,y=column_value,color=inputx1)

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
    run_app(app=app, launch_browser=False)
