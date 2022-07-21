# third party
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc

# project
from app import app, prefix_url
from .sidebar import sidebar
from .page1 import page1
from .page2 import page2

# css styling for page
CONTENT_STYLE = {
    # "margin-top": 0,
    #"margin-left": sidebar_width,
    # "margin-right": 0,
    # "margin-bottom": 0,
    "padding": "1rem 1rem",
}

dashboard_layout = html.Div(
    [
        dcc.Location(id="promotion-url", refresh=False),
        dcc.Store(id="query-string-store"),
        dbc.Container(
            dbc.Row(
                [
                    dbc.Col(html.Div(id='sidebar-content'), width=2),
                    dbc.Col(html.Div(id='page-content', style=CONTENT_STYLE), width=10),
                ]
            ),
            fluid=True

        ),
    ],
    className="w-100 pl-2",
)

def fetch_query_string(href: str) -> str:
    """
    Parses the href for any query string, returning it if found
    """
    split_url = href.split("?")
    if len(split_url) == 2:
        return f"?{split_url[1]}"
    else:
        return ""

# callback for sidebar
@app.callback(
    Output("sidebar-content", "children"),
    [Input("query-string-store", "data")],
)
def render_sidebar(query_string):
    children = sidebar(query_string)
    return children

# callback for sidebar navigation
@app.callback(
    [Output("page-content", "children"), Output("query-string-store", "data")],
    [Input("promotion-url", "href"), Input("promotion-url", "pathname")],
)
def render_page_content(href, pathname):
    query_string = fetch_query_string(href)
    if f"{prefix_url}page1" in pathname:
        layout = page1()
    elif f"{prefix_url}page2" in pathname:
        layout = page2()
    else:
        layout = page1()

    return layout, query_string

