# Import necessary libraries
from dash import html
import dash_bootstrap_components as dbc

# Define the navbar structure
def Navbar():

    layout = html.Div([
        dbc.NavbarSimple(
            children=[
                dbc.NavItem(dbc.NavLink("Mutations", href="/page1")),
                dbc.NavItem(dbc.NavLink("Transcription Factors", href="/page2")),
            ] ,
            brand="Heterogeneity of the Cancer Cell Line Metabolic Landscape",
            brand_href="/page1",
            color="dark",
            dark=True,
        ),
    ])

    return layout
