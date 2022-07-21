from dash import Dash, html, dcc, Output, Input
import dash, dash_table
import dash_bootstrap_components as dbc
import dash_daq as daq
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns

from matplotlib.patches import Rectangle

app = Dash(__name__)
## Read in data
shorthouse_data = pd.read_csv("./Data/Mutation_metabolite_associations_ordered_shorthouse.csv", index_col = "Unnamed: 0")
## Manipulate dataframe
original_index = shorthouse_data.index.tolist()
original_index = [str(i) for i in original_index]
shorthouse_data = shorthouse_data.reset_index()
new_index = shorthouse_data.index.tolist()
shorthouse_data = shorthouse_data.drop("index", axis =1)
## Generate colour scheme
palette_cmap = sns.diverging_palette(237, 8.7, s= 99, l = 50, as_cmap = True)
palette_cmap._init()
rgbas = palette_cmap._lut
hexes = [matplotlib.colors.rgb2hex(x) for x in rgbas]

## Plot heatmap
heatmap_big = px.imshow(shorthouse_data,y=original_index, x = shorthouse_data.columns, zmin = -15, zmax = 15, color_continuous_scale=list(hexes[:-3]),
                labels=dict(x="Gene", y="Metabolite", color="T-Statistic"), aspect="auto", title = "T-statistics for Mutation/Metabolite Pairings")

## Read in metabolite table for display
metabolite_lookup = pd.read_csv("./Data/Metabolite_reference_table_shorthouse.csv")
metabolite_lookup = metabolite_lookup[["ionIdx","id", "score", "name"]]

## Read in original data
metabolite_levels = pd.read_csv("./Data/Metabolite_levels_shorthouse.csv", index_col = "ionIdx")
metabolite_levels = metabolite_levels.drop("ionMz", axis =1)

## function for generating metabolite name string from id
def metaboname(metaboid, headnumber =2):
    metabolite_dataframe = metabolite_lookup
    metabolite_names = metabolite_dataframe[metabolite_dataframe["ionIdx"] == metaboid]
    metabolite_names = metabolite_names["name"].head(headnumber).str.cat(sep = "/")
    metabolite_string = str(metaboid) + ": " + metabolite_names

    return metabolite_string[:50]

## Generating dictionary for dropdown:
metabolite_dropdown = [{'label': metaboname(i), 'value': i} for i in shorthouse_data.index.tolist()[1:]]

## Generating name series for labels

metabolite_series = pd.Series([metaboname(i, 1) for i in shorthouse_data.index.tolist()], name = "metabolite")

## Reading in differential expression metabolite data
metabolite_diff_expr = pd.read_csv("./Data/Mutation_differential_expression_shorthouse.csv", index_col = "ionIdx")

## Reading in mutational data for celllines
celline_mutation_database = pd.read_csv("./Data/Mutations_in_celllines.csv")

## Read in mapping for celllines
celline_mapping = pd.read_csv("./Data/Cellline_mappings.csv")

#-------------------------------------
# Layout


layout = html.Div(children=[
    html.Br() ,
    html.H1(children='Influence of Mutations on Metabolite Abundance',style={'textAlign': 'center'}),

    html.Div(children='''
        This page contains plots to explore the relationships between nonsynonymous mutations and metabolites. Included are a heatmap of T-statistics (proportional to a p-value) for a logistic regression run on every nonsynonymous mutation/metabolite pairing. Click and drag to zoom.
        Scroll down for the metabolite table, and to explore the relationships between specific metabolites and mutations using the dropdown menus.
    ''',style={'textAlign': 'center'}),

    html.Br(),
    html.Div(children='''
        This data has been normalised in two differing ways - please see the relevant publications for details, but toggle between them below - default is Shorthouse et al.
        ''',style={'textAlign': 'center'}),

    dcc.RadioItems(
    options=[
       {'label': 'Shorthouse et al  ', 'value': 'shorthouse'},
       {'label': 'Cherkaoui et al', 'value': 'Cherkaoui'}], value = 'shorthouse',
    inline=True, style={'textAlign': 'center'}, inputStyle={"margin-right": "5px", "margin-left": "5px"},
    id = 'dataset_type'),

    ## Heatmap figure
    dcc.Graph(
        id='heatmap-top',
        figure=heatmap_big,
        ## add box-shadow below
        style = {'padding': 10}
    ),

    html.Div([

        # Graph container
        html.Div([

            ## Metabolite table
            dash_table.DataTable(
                id='table1',
                columns = [{"name": i, "id": i} for i in metabolite_lookup.columns],
                style_cell={'textAlign':'center','minWidth': 95, 'maxWidth': 95, 'width': 95,'font_size': '12px','whiteSpace':'normal','height':'auto'},
                data = metabolite_lookup.to_dict("rows"),
                #style_table={'overflow':'scroll','height':550},
                fixed_rows={'headers': True, 'data': 0},
                fixed_columns={'headers': True, 'data': 0},
                fill_width=False,
                style_cell_conditional=[
                    {
                        'if': {'column_id': "name"},
                        'minWidth': 400,
                        'width': 400,
                        "maxWidth": 400
                    },
                ],
            ),

        ], style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

        # Table container
        html.Div([

    ## Dropdown for metabolites
            dcc.Dropdown(id = "metabolite-id"
                 #options = [{'label': i, 'value': i} for i in shorthouse_data.index.tolist()],
                 ,options = metabolite_dropdown
                 ,value = 1
                 ,searchable = True
                 ,placeholder = "Peak id..."
                 ,clearable = True
                 ),
    ## Metabolite ranking graph
            dcc.Graph(id = "mutation_ranking_per_metabolite"),

        ], style={'width': '52%', 'display': 'inline-block', 'padding': 10}),

    ], style={'display': 'flex'}),




    ## Dropdown for mutations
    html.Div([
        html.Div([
            dcc.Dropdown(id = "mutation_id",
                 options = [{"label": i, "value": i} for i in metabolite_diff_expr.columns],
                 value = "A1CF",
                 placeholder = "Gene",
                 clearable = True,
                 searchable = True),
    ## Mutation volcano plot graph
            dcc.Graph(id = "mutation_volcano_plot")], style={'width': '49%', 'display': 'inline-block', 'padding': 10}),

    ## Metabolite/Mutation swarmplot
        html.Br() ,
        html.Div([
            dcc.Graph(id = "swarmplot_metabolite"),
            ], style={'width': '52%', 'display': 'inline-block', 'padding': 10}),
        ], style = {'display': 'flex'})

])

#-------------------------------------------
## Callbacks

## Callback for dynamic dropdown ## Doesnt seem to work
#@app.callback(
#    Output("metabolite-id", "metabolite_dropdown"),
#    Input("metabolite-id", "search_value")
#)
#def update_options(search_value):
#    if not search_value:
##        raise PreventUpdate
#    return [o for o in metabolite_dropdown if search_value in o["label"]]

## Callback for mutation rankings per metabolite
@app.callback(
    Output(component_id = "mutation_ranking_per_metabolite", component_property = "figure"),
    [Input(component_id = "metabolite-id", component_property = "value")]
)

## Function to extract metabolite column from dataframe and plot scatterplot
def mutation_ranking_per_metabolite_plot(metabolite_id_value):
    tstat_dataframe = shorthouse_data
    metabolite_names = metaboname(metabolite_id_value)
    metabolite = pd.DataFrame(tstat_dataframe.loc[metabolite_id_value].sort_values())
    metabolite.columns = ["T-Statistic"]
    metabolite["Gene"] = metabolite.index
    metabolite = metabolite.reset_index()
    metabolite["Mutation Rank"] = metabolite.index+1

    scatterplot = px.scatter(metabolite, x = "Mutation Rank", y = "T-Statistic", hover_name = "Gene"
                             ,title = "Mutation rankings for " + metabolite_names, template="simple_white")

    return scatterplot

## Callback for volcano plot based on mutations
@app.callback(
    Output(component_id = "mutation_volcano_plot", component_property = "figure"),
    [Input(component_id = "mutation_id", component_property = "value")]
)

## Function to generate data for plotting volcano plot
def volcano_plot_per_mutation(mutation_name):
    # Get data
    tstat_dataframe = shorthouse_data
    tstat_sorted = tstat_dataframe.sort_index()
    diff_expr = metabolite_diff_expr
    metabolite_series_names = metabolite_series
# get mutation of interest and merge frames
    mutation_tstat = tstat_sorted[mutation_name]
    diff_expr_mutation = diff_expr[mutation_name]
    mutation_tstat = abs(mutation_tstat)
    plotting_frame = pd.concat([mutation_tstat, diff_expr_mutation, metabolite_series_names], axis=1)
    plotting_frame.columns = ["tstat", "diffexpr", "metabolite"]

# Generate color column
    plotting_frame["colour"] = np.where((plotting_frame["tstat"] >= 5) & (plotting_frame["diffexpr"] >0), "Highly Increased",
                                        np.where((plotting_frame["tstat"] >= 5) & (plotting_frame["diffexpr"] <0), "Highly Decreased",
                                                 "Neutral"))
    volcanoplot = px.scatter(plotting_frame, y = "tstat", x = "diffexpr", color = "colour",
                             color_discrete_sequence = ["grey", "blue", "red"], hover_name = "metabolite",
                             labels={
                                 "tstat": "T-Statistic",
                                 "diffexpr": "Metabolite log(10) Difference",
                                 "colour": "Effect"
                             }, title = "Volcano plot for metabolite changes associated with " + mutation_name,
                             template="simple_white"
                             )
    return volcanoplot

## Callback for mutation swarmplot
@app.callback(
    Output(component_id = "swarmplot_metabolite", component_property = "figure"),
    [Input(component_id = "metabolite-id", component_property = "value")
        , Input(component_id = "mutation_id", component_property = "value")]
)

def swarmplot_per_metabolite_permutation(metabolite_id_value, mutation_gene):
# Copy dataframes
    metabolite_levels_base = metabolite_levels
    cellline_mutations = celline_mutation_database
    cellline_mappings_2 = celline_mapping

    cellline_mappings_short = cellline_mappings_2[["dsIdx", "ID"]]
    cellline_mappings_short["dsIdx"] = cellline_mappings_short["dsIdx"].apply(pd.to_numeric)
# Get only the metabolite we want
    metabolite_only = pd.DataFrame(metabolite_levels_base.T[metabolite_id_value])
    metabolite_only["dsIdx"] = metabolite_only.index
    metabolite_only["dsIdx"] = metabolite_only["dsIdx"].apply(pd.to_numeric)
    metabolite_cellines = metabolite_only.merge(cellline_mappings_short, left_on = "dsIdx", right_on = "dsIdx")
#Get list of celllines
    celllines_list = metabolite_cellines["ID"].unique().tolist()

    ## now have metabolite level for each cell line TAG - need to map to mutational status.

#Get only the celllines which have mutations in gene of interest
    gene_of_interest = cellline_mutations[cellline_mutations["HGNC"] == mutation_gene]
    gene_of_interest = gene_of_interest[gene_of_interest["MutationType"] != "Silent"]
    gene_of_interest = gene_of_interest[gene_of_interest["MutationType"].notna()]
    #print(gene_of_interest)
## Get only celllines in our dataset
    celllines_with_gene_of_interest = gene_of_interest[gene_of_interest["CellLineName_Cellosaurus"].isin(celllines_list)]
## Merge mutations
    celllines_with_gene_of_interest["Mutation"] = celllines_with_gene_of_interest.groupby(['CellLineName_Cellosaurus'])['AA_Mutation'].transform(lambda x : ', '.join(x.unique()))
    celllines_with_gene_of_interest = celllines_with_gene_of_interest.drop_duplicates(subset = "Mutation")[["CellLineName_Cellosaurus", "Mutation"]]
    celllines_with_gene_of_interest.columns = ["ID", "Mutation"]
    #print(celllines_with_gene_of_interest)
## Add a column that counts up the mutations for colouring
    celllines_with_gene_of_interest["Mutant"] = range(1,len(celllines_with_gene_of_interest)+1)
## Merge our values with our labels
    plotting_dataframe = metabolite_cellines.merge(celllines_with_gene_of_interest, how = "outer", left_on = 'ID',
                                                   right_on = "ID")
## Fill the dataframe NaNs with 0 or "-"
    plotting_dataframe["Mutant"] = plotting_dataframe["Mutant"].fillna(0)
    plotting_dataframe["Mutation"] = plotting_dataframe["Mutation"].fillna(" -")
## Log 10 transform the data
    plotting_dataframe[metabolite_id_value] = np.log10(plotting_dataframe[metabolite_id_value])
## Get the name of the metabolite
    metabolite_name = metaboname(metabolite_id_value)
    swarmplot = px.strip(plotting_dataframe[metabolite_id_value], color = plotting_dataframe["Mutation"],
                         stripmode = "overlay", hover_name = plotting_dataframe["ID"],
                         labels=dict(color= mutation_gene + " mutant", value="log(10) metabolite expression",
                                     variable = metabolite_name),
                                     template="simple_white",
                                     title = "Expression of " + metabolite_name + " in comparison to mutations in " + mutation_gene)
    return(swarmplot)


if __name__ == '__main__':
    app.run_server(debug=True)
