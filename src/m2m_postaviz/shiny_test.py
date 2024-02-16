import plotly.express as px
import m2m_postaviz.data_utils as du
import pandas as pd
import os.path
from shiny import App, render, run_app, ui, reactive
import m2m_postaviz.data_utils as du
from shinywidgets import output_widget
from shinywidgets import render_widget
import json


def multiply_production_abundance(df_row,abundance_matrix,sample_id):
    print(df_row)
    return


def generate_abundance_matrix(sample_matrix: pd.DataFrame, abundance_matrix: pd.DataFrame, sample_id):
    abundance_matrix_col = abundance_matrix[sample_id]
    relative_abundance_matrix = sample_matrix.apply(lambda row: multiply_production_abundance(row,abundance_matrix,abundance_matrix_col),axis=1) 
    return


### Import data

dir_path = "/home/lbrindel/cscope_metadata/"
metadata = "/home/lbrindel/m2m-postaviz/tests/metadata_test.tsv"

data, sample_data = du.build_df(dir_path,metadata)

mgs_data = pd.read_csv("~/Downloads/specI.mat", sep="\t")
mgs = mgs_data.columns.values[0]

palleja_abundance_matrix = du.open_tsv("~/Downloads/matrix_palleja.tsv")
palleja_abundance_matrix.insert(0,"bin_id",palleja_abundance_matrix.index)
palleja_abundance_matrix.reset_index

########################################

### App object

sample_choice_type = ui.input_select("dataframe_type_choice","Cscope or Iscope", choices=["cscope","iscope"], selected="cscope")
sample_choice = ui.input_select("dataframe_sample_choice","Choose sample's dataframe to show", choices=[sample for sample in sample_data.keys()])

sample_main_dataframe = ui.card(
    ui.output_data_frame("main_dataframe"),
    ui.output_data_frame("cross_dataframe")
)

matrix_palleja = ui.output_data_frame("matrix_palleja")

### App tree

app_ui = ui.page_fluid(
    ui.row(sample_choice,sample_choice_type),
    ui.accordion(
        ui.accordion_panel("sample and SpecI cross dataframe",
            sample_main_dataframe
        ),
        ui.accordion_panel("Matrix palleja dataframe",
            matrix_palleja
        ),
        ui.accordion_panel("print",
            ui.output_text("printo")
        ),
        open=False
    )
)

def server(input, output, session):

    @render.data_frame
    def main_dataframe():
        return render.DataGrid(sample_data[input.dataframe_sample_choice()][input.dataframe_type_choice()],row_selection_mode="single")
    

    @render.data_frame
    def cross_dataframe():
        cross_df = mgs_data.loc[mgs_data[mgs].isin(sample_data[input.dataframe_sample_choice()][input.dataframe_type_choice()]["Name"])]
        return render.DataGrid(mgs_data.loc[mgs_data[mgs].isin(sample_data[input.dataframe_sample_choice()][input.dataframe_type_choice()]["Name"])])


    @render.data_frame
    def matrix_palleja():
        df = palleja_abundance_matrix
        return render.DataGrid(df)


    @render.text
    def printo():
        return generate_abundance_matrix(sample_data["E1_0"]["cscope"],palleja_abundance_matrix,"ERAS1d0")

app = App(app_ui, server)
run_app(app=app, launch_browser=True)





