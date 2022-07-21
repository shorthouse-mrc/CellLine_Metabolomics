from dash import Dash, html, dcc, Output, Input
import dash, dash_table
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns

app = Dash(__name__)

## Read in data
metabolomics_TF_correlations_shorthouse = pd.read_csv("./Data/Progeny_correlations_shorthouse.csv", index_col = "Pathway")
metabolomics_TF_correlations_cherkaoui = pd.read_csv("./Data/Progeny_correlations_cherkaoui.csv", index_col = "Pathway")


# Get the top 50 pathways by correlation
metabolomics_TF_correlations_shorthouse = -metabolomics_TF_correlations_shorthouse
metabolomics_TF_correlations_cherkaoui = -metabolomics_TF_correlations_cherkaoui

# Rename columns (remove start of name)
newcolumns = []
for item in metabolomics_TF_correlations_shorthouse.columns:
    splititem = item.split("_")[1]
    newcolumns.append(splititem)
metabolomics_TF_correlations_shorthouse.columns = newcolumns

newcolumns = []
for item in metabolomics_TF_correlations_cherkaoui.columns:
    splititem = item.split("_")[1]
    newcolumns.append(splititem)
metabolomics_TF_correlations_cherkaoui.columns = newcolumns

## Generate colourmap
palette_cmap = sns.light_palette(sns.color_palette("Set3")[3], as_cmap = True)
palette_cmap._init()
rgbas = palette_cmap._lut
hexes = [matplotlib.colors.rgb2hex(x) for x in rgbas]

# Load in specific TF/pathway correlations
TF_correlations_shorthouse = pd.read_csv("./Data/TF_pathway_correlations_shorthouse.csv", index_col = "Pathway")
TF_columnnames = []
for item in TF_correlations_shorthouse.columns:
    TF_name = item.split("_")[1]
    TF_columnnames.append(TF_name)
TF_correlations_shorthouse.columns = TF_columnnames

TF_correlations_cerkaoui = pd.read_csv("./Data/TF_pathway_correlations_cherkaoui.csv", index_col = "Pathway")
TF_columnnames = []
for item in TF_correlations_cerkaoui.columns:
    TF_name = item.split("_")[1]
    TF_columnnames.append(TF_name)
TF_correlations_cerkaoui.columns = TF_columnnames



### ----------------------
# Layout
layout = html.Div(children=[
    html.Br() ,
    html.H1(children='Correlations between transcription factor (TF) activity and metabolic pathways',style={'textAlign': 'center'}),

    html.Div(children='''
        This page contains plots to explore the relationships between transcription factors (TFS) and SMPDB metabolic pathways.
        The top of the page is a heatmap of the top pathway/PROGENY associations, scroll down to explore the correlations between
        specific SMPDB pathways and transcription factors using the dropdown menus.
    ''',style={'textAlign': 'center'}),

    html.Br(),
    html.Div(children='''
        This data has been normalised in two differing ways - please see the relevant publications for details, but toggle between them below - default is Shorthouse et al.
    ''',style={'textAlign': 'center'}),

    dcc.RadioItems(
    options=[
        {'label': 'Shorthouse et al  ', 'value': 'shorthouse'},
        {'label': 'Cherkaoui et al', 'value': 'cherkaoui'}], value = 'shorthouse',
        inline=True, style={'textAlign': 'center'}, inputStyle={"margin-right": "5px", "margin-left": "5px"},
        id = 'dataset_type'),
    ## Heatmap figure
    dcc.Graph(
        id='heatmap_top',
        style = {'padding': 10}
    ),

    html.Div([

        # Graph container
        html.Div([
    ## Dropdown for pathways
            dcc.Dropdown(id="pathway"
                 , value="Citric Acid Cycle"
                 , searchable=True
                 , clearable=True),
                 ## Pathway ranking graph
            dcc.Graph(id="TF_ranking_per_pathway"),
            ],style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

        html.Div([
            ## Dropdown for TFs
            dcc.Dropdown(id="TF"
                 , value="AR"
                 , searchable=True
                 , placeholder="Choose Transcription Factor..."
                 , clearable=True),
    ## TF ranking graph
    dcc.Graph(id="pathway_ranking_per_TF"),
    ],style={'width': '52%', 'display': 'inline-block', 'padding': 10})
    ], style={'display': 'flex'}),
])

#-------------------------------------------
## Callbacks
## Callback for heatmap
@app.callback(
    Output(component_id = 'heatmap_top', component_property = 'figure'),
    Input(component_id = 'dataset_type', component_property = 'value')
)

def heatmap_TFS_plot(dataset):
    if dataset == "shorthouse":
        metabolomics_TF_correlations = metabolomics_TF_correlations_shorthouse
    elif dataset == "cherkaoui":
        metabolomics_TF_correlations = metabolomics_TF_correlations_cherkaoui
    heatmap_TFS = px.imshow(metabolomics_TF_correlations, color_continuous_scale=list(hexes[:-3]),
                     aspect="auto", labels={
                         "y": ""})
    heatmap_TFS.update_coloraxes(colorbar_exponentformat="power")

    heatmap_TFS.update_traces(hovertemplate='Progeny Signature: %{x} <br>SMPDB Pathway: %{y}<br>-log10(P value): %{z}', zhoverformat = "power")
    return heatmap_TFS

## Callback for Dropdown for pathways
@app.callback(
    Output(component_id = "pathway", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)
def set_dropdown_options_1(dataset):
    if dataset == "shorthouse":
        pathway_dropdown = [{'label': i, 'value': i} for i in TF_correlations_shorthouse.index.tolist()[1:]]
    elif dataset == "cherkaoui":
        pathway_dropdown = [{'label': i, 'value': i} for i in TF_correlations_cerkaoui.index.tolist()[1:]]
    return pathway_dropdown

## Callback for Dropdown for TFS
@app.callback(
    Output(component_id = "TF", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)
def set_dropdown_options_2(dataset):
    if dataset == "shorthouse":
        pathway_dropdown = [{'label': i, 'value': i} for i in TF_correlations_shorthouse.columns]
    elif dataset == "cherkaoui":
        pathway_dropdown = [{'label': i, 'value': i} for i in TF_correlations_cerkaoui.columns]
    return pathway_dropdown

## Callback for mutation rankings per metabolite
@app.callback(
    Output(component_id = "TF_ranking_per_pathway", component_property = "figure"),
    [Input(component_id = "pathway", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def TF_ranking_by_pathway_id_plot(pathway_name, dataset):
    if dataset == "shorthouse":
        pathway_dataframe = -TF_correlations_shorthouse
    elif dataset == "cherkaoui":
        pathway_dataframe = -TF_correlations_cerkaoui
    pathway = pd.DataFrame(pathway_dataframe.loc[pathway_name].sort_values())
    pathway.columns = ["-log10(Pvalue)"]
    pathway["Transcription Factor"] = pathway.index
    pathway = pathway.reset_index()
    pathway["TF Rank"] = pathway.index+1

    scatterplot = px.scatter(pathway, x = "TF Rank", y = "-log10(Pvalue)", hover_name = "Transcription Factor"
                             ,title = "Ranks of TFs against " + pathway_name + " activity", template = "simple_white")
    return scatterplot

@app.callback(
    Output(component_id = "pathway_ranking_per_TF", component_property = "figure"),
    [Input(component_id = "TF", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def pathway_ranking_by_TF_id_plot(TF_name, dataset):
    if dataset == "shorthouse":
        pathway_dataframe = -TF_correlations_shorthouse.T
    elif dataset == "cherkaoui":
        pathway_dataframe = -TF_correlations_cerkaoui.T
    pathway = pd.DataFrame(pathway_dataframe.loc[TF_name].sort_values())
    pathway.columns = ["-log10(Pvalue)"]
    pathway["Pathway Name"] = pathway.index
    pathway = pathway.reset_index()
    pathway["Pathway Rank"] = pathway.index+1

    scatterplot = px.scatter(pathway, x = "Pathway Rank", y = "-log10(Pvalue)", hover_name = "Pathway"
                             ,title = "Ranks of pathways against " + TF_name + " activity", template = "simple_white")
    return scatterplot

if __name__ == '__main__':
    app.run_server(debug=True, port = 8051)
