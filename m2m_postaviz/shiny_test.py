import plotly.express as px
import m2m_postaviz.data_utils as du
import pandas as pd
import os.path
from shiny import App, render, run_app, ui, reactive
import m2m_postaviz.data_utils as du
from shinywidgets import output_widget
from shinywidgets import render_widget
import json

def sum_abundance_table(abundance_table: pd.DataFrame):
    new_dataframe =  {}
    for index,row in abundance_table.iterrows():
        new_dataframe[index] = row.values.sum()
    return pd.DataFrame(new_dataframe, index=["E1"])


def multiply_production_abundance(row: pd.Series, abundance_matrix: pd.DataFrame,sample_id):
    row = row.astype(float)
    abundance_value = abundance_matrix.at[row.name,sample_id]
    row *= abundance_value
    return row


def generate_stoichiometric_matrix(binary_matrix: pd.DataFrame, abundance_matrix: pd.DataFrame, sample_id: str):
    binary_matrix.set_index("Name",inplace=True)
    binary_matrix = binary_matrix.apply(lambda row: multiply_production_abundance(row, abundance_matrix,sample_id),axis=1) 
    return binary_matrix


### Import data

palleja_abundance_matrix = du.open_tsv("~/Downloads/matrix_palleja.tsv")
# palleja_abundance_matrix.insert(0,"bin_id",palleja_abundance_matrix.index)
# palleja_abundance_matrix.reset_index

########################################

### App object

# sample_choice_type = ui.input_select("dataframe_type_choice","Cscope or Iscope", choices=["cscope","iscope"], selected="cscope")
# sample_choice = ui.input_select("dataframe_sample_choice","Choose sample's dataframe to show", choices=[sample for sample in sample_data.keys()])

# sample_main_dataframe = ui.card(
#     ui.output_data_frame("main_dataframe"),
#     ui.output_data_frame("cross_dataframe")
# )

matrix_palleja = ui.output_data_frame("matrix_palleja")
# plot_palleja = output_widget("plot_palleja")

### App tree

app_ui = ui.page_fluid(
    ui.accordion(
        ui.accordion_panel("Matrix palleja dataframe",
            matrix_palleja,
            # plot_palleja
    )
))

def server(input, output, session):

    @render.data_frame
    def matrix_palleja():
        df = palleja_abundance_matrix
        return render.DataGrid(df)

    @output
    @output_widget
    def plot_palleja():
        df= palleja_abundance_matrix
        plot = px.box(df,y='ERAS1d0')
        return plot


app = App(app_ui, server)
run_app(app=app, launch_browser=True)





