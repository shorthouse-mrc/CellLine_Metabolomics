# Import necessary libraries
from dash import html, dcc
from dash.dependencies import Input, Output

# Connect to main app.py file
from app import app

# Connect to your app pages
from pages import page1, page2, page3

# Connect the navbar to the index
from components import navbar

# Define the navbar
nav = navbar.Navbar()

app.config.suppress_callback_exceptions = True


# Define the index page layout
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    nav,
    html.Div(id='page-content', children=[]),
])

## Callbacks for page 1

## Callback for tstat heatmap
@app.callback(
    Output(component_id = 'heatmap_top1', component_property = 'figure'),
    Input(component_id = 'dataset_type', component_property = 'value')
)
def plot_heatmap_tstats(dataset):
    graph = page1.plot_heatmap_tstats(dataset)
    return graph

## Callback for table1
@app.callback(
    Output(component_id = "metabolite_table", component_property = "data"),
    [Input(component_id = "dataset_type", component_property = "value")]
)
def generate_table(dataset):
    tabledata = page1.generate_tabledata(dataset)
    return tabledata

## Callback for updating table 1
@app.callback(Output('table1', 'data'),
              Input('metabolite_table', 'data'))
def update_table(dataset):
    tableupdate = page1.on_data_set_table(dataset)
    return tableupdate

## Callback for updating dropdowns
@app.callback(
    Output(component_id = "metabolite_id", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)
def update_dropdown_page1_1(dataset):
    dropdown_values = page1.dropdown_update(dataset)
    return dropdown_values

## Callback for mutation rankings per metabolite
@app.callback(
    Output(component_id = "mutation_ranking_per_metabolite", component_property = "figure"),
    [Input(component_id = "metabolite_id", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)
def update_graph1(metabolite_id_value, dataset):
    graph = page1.mutation_ranking_per_metabolite_plot(metabolite_id_value, dataset)
    return graph

## Callback for volcano plot of metabolites against mutations
@app.callback(
    Output(component_id = "mutation_volcano_plot", component_property = "figure"),
    [Input(component_id = "mutation_id", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)
def update_graph2(mutation_name, dataset):
    graph = page1.volcano_plot_per_mutation(mutation_name, dataset)
    return graph

## Callback for mutation swarmplot
@app.callback(
    Output(component_id = "swarmplot_metabolite", component_property = "figure"),
    [Input(component_id = "metabolite_id", component_property = "value"),
    Input(component_id = "mutation_id", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def update_graph3(metabolite_id_value, mutation_gene, dataset):
    graph = page1.swarmplot_per_metabolite_permutation(metabolite_id_value, mutation_gene, dataset)
    return graph

## Callbacks for page 2

## Callback for TF heatmap_big
@app.callback(
    Output(component_id = 'heatmap_top', component_property = 'figure'),
    Input(component_id = 'dataset_type', component_property = 'value')
)
def update_graph_TFheatmap(dataset):
    graph = page2.heatmap_TFS_plot(dataset)
    return graph

## Callback for dropdown for pathways
@app.callback(
    Output(component_id = "pathway", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)

def update_dropdownpathways(dataset):
    dropdownlist = page2.set_dropdown_options_1(dataset)
    return dropdownlist

## Callback for dropdown for TFs
@app.callback(
    Output(component_id = "TF", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)

def update_dropdownTFS(dataset):
    dropdownlist = page2.set_dropdown_options_2(dataset)
    return dropdownlist

## Callback for mutation rankings per metabolite
@app.callback(
    Output(component_id = "TF_ranking_per_pathway", component_property = "figure"),
    [Input(component_id = "pathway", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)
def update_graph4(pathway_name, dataset):
    graph = page2.TF_ranking_by_pathway_id_plot(pathway_name, dataset)
    return graph

## Callback for TF ranking per pathway
@app.callback(
    Output(component_id = "pathway_ranking_per_TF", component_property = "figure"),
    [Input(component_id = "TF", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def update_graph5(TF_name, dataset):
    graph = page2.pathway_ranking_by_TF_id_plot(TF_name, dataset)
    return graph

## Callbacks for page 3

## Callbacks for dropdowns
@app.callback(
    Output(component_id = "pathway2", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)

def update_dropdown_page3_1(dataset):
    dropdownlist = page3.set_dropdown_options_page3_1(dataset)
    return dropdownlist

@app.callback(
    Output(component_id = "drug", component_property = "options"),
    [Input(component_id = "dataset_type", component_property = "value")]
)

def update_dropdown_page3_2(dataset):
    dropdownlist = page3.set_dropdown_options_page3_2(dataset)
    return dropdownlist

## Callback for drug sensitivity per pathway
@app.callback(
    Output(component_id = "drug_sensitivity_by_pathway", component_property = "figure"),
    [Input(component_id = "pathway2", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def update_graph6(pathway_name, dataset):
    graph = page3.drug_sensitivity_by_pathway_plot(pathway_name, dataset)
    return graph

## Callback for pathway rankings per drug_name
@app.callback(
    Output(component_id = "pathway_ranking_by_drug", component_property = "figure"),
    [Input(component_id = "drug", component_property = "value"),
    Input(component_id = "dataset_type", component_property = "value")]
)

def update_graph7(drug_name, dataset):
    graph = page3.pathway_ranking_by_drug_plot(drug_name, dataset)
    return graph




@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/page1':
        return page1.layout
    if pathname == '/page2':
        return page2.layout
    if pathname == '/page3':
        return page3.layout
    else: # if redirected to unknown link
        return "404 Page Error! Please choose a link"








# Run the app on localhost:8050
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port='80')
