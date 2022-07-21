# third party
import dash_bootstrap_components as dbc
from dash import html
from app import prefix_url


SIDEBAR_STYLE = {
    # "position": "fixed",
    # "width": "15%",
    "padding": "0.25rem 0.25rem",
    "background-color": "#f9f9f9",
    "border-right": "1px solid #cacaca",
    "overflow": "scroll",
}


def sidebar(query_string=""):
    return html.Div(
        [
            dbc.Nav(
                [
                    html.H5(
                        "Some Analysis",
                        style={'color': "#000000"}
                    ),
                    dbc.NavLink(
                        "📉 Page 1",
                        href=f"{prefix_url}page1/{query_string}",
                        id="page1-nav",
                    ),
                    dbc.NavLink(
                        "🧮 Page 2",
                        href=f"{prefix_url}page2/{query_string}",
                        id="page2-nav",
                    ),
                ],
                vertical=True,
                pills=True,
            )
        ],
        style=SIDEBAR_STYLE,
    )

