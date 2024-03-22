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

dir_path = "/home/lbrindel/cscope_metadata/"
metadata = "/home/lbrindel/m2m-postaviz/tests/metadata_test.tsv"

data, sample_data = du.build_df(dir_path,metadata)

mgs_data = pd.read_csv("~/Downloads/specI.mat", sep="\t")
mgs = mgs_data.columns.values[0]

palleja_abundance_matrix = du.open_tsv("~/Downloads/matrix_palleja.tsv")
# palleja_abundance_matrix.insert(0,"bin_id",palleja_abundance_matrix.index)
# palleja_abundance_matrix.reset_index

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
        ui.accordion_panel("Abundance",
            ui.output_data_frame("abundance_table")
        ),
        ui.accordion_panel("Abundance Boxplot",
            output_widget("abundance_plot")
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
        return render.DataGrid(mgs_data.loc[mgs_data[mgs].isin(sample_data[input.dataframe_sample_choice()][input.dataframe_type_choice()]["Name"])])


    @render.data_frame
    def matrix_palleja():
        df = palleja_abundance_matrix
        return render.DataGrid(df)


    @render.data_frame
    def abundance_table():
        return render.DataGrid(generate_stoichiometric_matrix(sample_data["E1_0"]["cscope"],palleja_abundance_matrix,"ERAS1d0"))


    @output
    @render_widget
    def abundance_plot():
        df = generate_stoichiometric_matrix(sample_data["E1_0"]["cscope"],palleja_abundance_matrix,"ERAS1d0")
        df = sum_abundance_table(df)
        fig = px.box(df, x=df.index,y=df.columns.values)
        return fig

app = App(app_ui, server)
run_app(app=app, launch_browser=True)





