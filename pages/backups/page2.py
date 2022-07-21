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
metabolomics_TF_correlations = pd.read_csv("./Data/Progeny_correlations_shorthouse.csv", index_col = "Pathway")

# Get the top 50 pathways by correlation
metabolomics_TF_correlations = -metabolomics_TF_correlations
binaried = pd.DataFrame(np.where(metabolomics_TF_correlations > -np.log10(0.05), 1, 0))
binaried.columns = metabolomics_TF_correlations.columns
binaried.index = metabolomics_TF_correlations.index
binaried_nozeroes = binaried.loc[~(binaried==0).all(axis=1)]
binaried_nozeroes["sum"] = binaried_nozeroes.sum(axis =1)
top_50 = binaried_nozeroes[binaried_nozeroes["sum"] >1].index
# Use the top correlations to subset the original dataframe
top_correlations = metabolomics_TF_correlations.loc[top_50]
#print(top_correlations)

# Rename columns (remove start of name)
newcolumns = []
for item in metabolomics_TF_correlations.columns:
    splititem = item.split("_")[1]
    newcolumns.append(splititem)
metabolomics_TF_correlations.columns = newcolumns

# Transform data
top_correlations_unlog = -top_correlations**10
## Generate colourmap
palette_cmap = sns.light_palette(sns.color_palette("Set3")[3], as_cmap = True)
palette_cmap._init()
rgbas = palette_cmap._lut
hexes = [matplotlib.colors.rgb2hex(x) for x in rgbas]

heatmap_TFS = px.imshow(metabolomics_TF_correlations, color_continuous_scale=list(hexes[:-3]),
                 aspect="auto", labels={
                     "y": ""})
heatmap_TFS.update_coloraxes(colorbar_exponentformat="power")

heatmap_TFS.update_traces(hovertemplate='Progeny Signature: %{x} <br>SMPDB Pathway: %{y}<br>-log10(P value): %{z}', zhoverformat = "power")

# Load in specific TF/pathway correlations
TF_correlations = pd.read_csv("./Data/TF_pathway_correlations.csv", index_col = "Pathway")
TF_columns = TF_correlations.columns
TF_columnnames = []
for item in TF_columns:
    TF_name = item.split("_")[1]
    TF_columnnames.append(TF_name)

TF_correlations.columns = TF_columnnames

pathway_dropdown = [{'label': i, 'value': i} for i in TF_correlations.index.tolist()[1:]]
TF_dropdown = [{'label': i, 'value': i} for i in TF_columnnames]

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
    ## Heatmap figure
    dcc.Graph(
        id='heatmap-top',
        figure=heatmap_TFS,
        style = {'padding': 10}
    ),

    html.Div([

        # Graph container
        html.Div([
    ## Dropdown for pathways
            dcc.Dropdown(id="pathway"
                 , options=pathway_dropdown
                 , value="Citric Acid Cycle"
                 , searchable=True
                 , placeholder="Choose pathway..."
                 , clearable=True),
                 ## Pathway ranking graph
            dcc.Graph(id="TF_ranking_per_pathway"),
            ],style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

        html.Div([
            ## Dropdown for TFs
            dcc.Dropdown(id="TF"
                 , options=TF_dropdown
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
## Callback for mutation rankings per metabolite
@app.callback(
    Output(component_id = "TF_ranking_per_pathway", component_property = "figure"),
    [Input(component_id = "pathway", component_property = "value")]
)

def TF_ranking_by_pathway_id_plot(pathway_name):
    pathway_dataframe = -TF_correlations
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
    [Input(component_id = "TF", component_property = "value")]
)

def pathway_ranking_by_TF_id_plot(TF_name):
    pathway_dataframe = -TF_correlations.T
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
