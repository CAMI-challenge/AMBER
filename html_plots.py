#!/usr/bin/env python

import pandas as pd
import bokeh.models.widgets.tables

from bokeh.layouts import widgetbox, layout
from bokeh.models import DataTable
from bokeh.palettes import d3
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models.widgets import Div
from bokeh.models import HoverTool
import numpy as np
import argparse
from utils.labels import abbreviations
import datetime

SCATTER_ELEMENT_WIDTH = 1500
SCATTER_ELEMENT_HEIGHT = 500

TABLE_ELEMENT_WIDTH = 1500
TABLE_ELEMENT_HEIGHT = 300

DIV_WIDTH = 1200
DIV_HEIGHT = None

DESCRIPTION_WIDTH = 600
DESCRIPTION_HEIGHT = 60

A_WIDTH = 1200
A_HEIGHT = 15

COLORS_20 = d3['Category20'][20]
COLORS_10 = d3['Category10'][10]

ID_SUMMARY = "summary"
ID_PRECISION_VS_RECALL_TOOLS = "precision_vs_recall_tools"
ID_RAND_INDEX_ASSIGNED_BPS = "rand_index_assigned_bps"
ID_ALL_GENOMES = "all_genomes"
ID_SINGLE_TOOL = "single_tool"

DESCRIPTION_PRECISION_RECALL_TOOLS = """
Scatter plot showing average precision vs average recall including standard error
of the mean.
"""

DESCRIPTION_PRECISION_RECALL_GENOME = """
Scatter plot showing precision vs recall per genome per tool.
"""

DESCRIPTION_PRECISION_RECALL_TOOL = """
Precision vs. recall plots for every tool. The first plot is sorted by precision and the second one by recall.
"""

DESCRIPTION_ADJUSTED_RAND_INDEX_TOOLS = """
Scatter plot showing adjusted rand index vs percentage of assigned base pairs.
"""


def get_color(number_of_tools):
    if number_of_tools > 10:
        return COLORS_20
    else:
        return COLORS_10


def create_title_div(id, name, info):
    div = Div(text="""<h1 id="{0}">{1}</h1><p>{2}</p>""".format(id, name, info),
              width=DIV_WIDTH, height=DIV_HEIGHT)
    return [div]


def _set_default_figure_properties(figure, x_label, y_label):
    figure.legend.click_policy = "hide"
    figure.legend.location = "top_right"
    figure.xaxis.axis_label = x_label
    figure.yaxis.axis_label = y_label
    figure.xaxis.major_label_orientation = "vertical"
    return figure


def create_description(content, width=DESCRIPTION_WIDTH, height=DESCRIPTION_HEIGHT):
    div = Div(text="""<h3>Description:</h3><p>{0}</p>""".format(content),
              width=width, height=height)
    return [div]


def create_subtitle_div(id, name):
    div = Div(text="""<h2 id="{0}">{1}</h2>""".format(id, name),
              width=DIV_WIDTH, height=DIV_HEIGHT)
    return [div]


def create_subtitle_a(href, name):
    div = Div(text="""<a href=#{0} style='padding-left:30px'>{1}</h2>""".format(href, name),
              width=A_WIDTH, height=A_HEIGHT)
    return [div]


def create_title_a(href, number, name):
    div = Div(text="""<a href=#{0}>{2}. {1}</h2>""".format(href, name, number),
              width=A_WIDTH, height=A_HEIGHT)
    return [div]


def create_summary_table(df):
    dt = DataTable(source=bokeh.models.ColumnDataSource(df),
                   columns=list(map(lambda x: bokeh.models.widgets.TableColumn(title=x, field=x), df.columns.values)),
                   width=TABLE_ELEMENT_WIDTH,
                   height=TABLE_ELEMENT_HEIGHT,
                   reorderable=True,
                   selectable=True)
    return [widgetbox(dt)]


def errorbar(fig, x, y, xerr=None, yerr=None, color='red',
             point_kwargs={}, error_kwargs={}):
    """
    Errorbar based on
    https://stackoverflow.com/questions/29166353/how-do-you-add-error-bars-to-bokeh-plots-in-python
    :param fig: Bokeh figure
    :param x: x coordinate
    :param y: y coordinate
    :param xerr: x axis error bar
    :param yerr: y axis error bar
    :param color: color of the error bar
    :param point_kwargs: further arguments provided to bokeh circle function
    :param error_kwargs: further arguments provided to bokeh multiline function
    """
    fig.circle(x, y, color=color, **point_kwargs)

    if xerr:
        x_err_x = []
        x_err_y = []
        for px, py, err in zip(x, y, xerr):
            x_err_x.append((px - err, px + err))
            x_err_y.append((py, py))
        fig.multi_line(x_err_x, x_err_y, color=color, **error_kwargs)

    if yerr:
        y_err_x = []
        y_err_y = []
        for px, py, err in zip(x, y, yerr):
            y_err_x.append((px, px))
            y_err_y.append((py - err, py + err))
        fig.multi_line(y_err_x, y_err_y, color=color, **error_kwargs)


def create_recall_precision_scatter(df):
    """
    Creates average precision vs average recall plot with error bars for all tools.
    :param df: Dataframe
    :return: bokeh figure
    """
    p = figure(x_range=[0, 1], y_range=[0 ,1], plot_width=SCATTER_ELEMENT_WIDTH, plot_height=SCATTER_ELEMENT_HEIGHT)

    for i, (idx, row) in zip(np.arange(len(df.index)), df.iterrows()):
        errorbar(p, [row["avg_precision"]], [row["avg_recall"]], xerr=[row["sem_precision"]],
                 yerr=[row["sem_recall"]], color=get_color(len(df.index))[i], point_kwargs={"legend": idx, "size": 10},
                 error_kwargs={"legend": idx})

    p = _set_default_figure_properties(p, "Precision", "Recall")
    return [p]


def create_rand_index_assigned_bps_scatter(df):
    """
    Creates average rand index vs percent of assigned base pairs plot with error bars for all tools.
    :param df: Dataframe
    :return: bokeh figure
    """
    p = figure(plot_width=SCATTER_ELEMENT_WIDTH, plot_height=SCATTER_ELEMENT_HEIGHT)

    for i, (idx, row) in zip(np.arange(len(df.index)), df.iterrows()):
        p.circle(row["a_rand_index_by_bp"], row["percent_assigned_bps"], alpha=0.8, color=get_color(len(df.index))[i], size=10,
                 legend=idx)

    p = _set_default_figure_properties(p, "Average Rand Index by base pair", "Percent of Assigned base pairs")
    return [p]


def create_precision_recall_all_genomes_scatter(paths, names):
    """
    Creates precision vs recall scatter plot for all tools and genomes.
    :param paths: Path to all precision and recall tables
    :return: bokeh figure
    """

    def assign_tool_name(path, name):
        df = pd.DataFrame.from_csv(path, sep='\t', header=0)
        df["name"] = name
        df["genome"] = df.index
        return df

    df = pd.concat(map(lambda x: assign_tool_name(*x), zip(paths, names)))

    p = figure(plot_width=SCATTER_ELEMENT_WIDTH, plot_height=SCATTER_ELEMENT_HEIGHT)
    p.add_tools(HoverTool(tooltips=[
        ("genome", "@genome")
    ]))
    names = list(df.name.unique())

    def create_circles(df, plot, name, color):
        df = df[df["name"] == name]
        source = ColumnDataSource(data=df)
        plot.circle("precision", "recall", color=color, alpha=0.8, legend=name, source=source)
        return plot

    for idx, name in enumerate(names):
        p = create_circles(df, p, name, get_color(len(names))[idx])

    p = _set_default_figure_properties(p, "Precision", "Recall")
    return [p]


def save_html_file(path, elements):
    html = file_html(layout(elements, sizing_mode='scale_width'), CDN, "AMBER: Assessment of Metagenome BinnERs")
    file = open(path, "w+")
    file.write(html)
    file.close()


def build_html(precision_recall_paths, names, summary, html_output):
    element_column = list()

    df = summary
    # if this script is run directly (not from evaluate.py), summary is a file path
    if isinstance(summary, str):
        df = pd.DataFrame.from_csv(summary, sep='\t', header=0)
    df.insert(0, "Tool", df.index)

    element_column.append(create_title_div("main", "AMBER: Assessment of Metagenome BinnERs", " produced on {0} ".format(
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))))
    element_column.append(create_subtitle_div("contents", "Contents"))

    element_column.append(create_title_a(ID_SUMMARY, "1", "Summary"))
    element_column.append(create_title_a(ID_PRECISION_VS_RECALL_TOOLS, "2", "Precision vs. Recall all Tools"))
    element_column.append(
        create_title_a(ID_RAND_INDEX_ASSIGNED_BPS, "3", "Adjusted Rand index vs. Percentage of Assigned base pairs"))
    element_column.append(create_title_a(ID_ALL_GENOMES, "4", "Precision vs. Recall all Genomes"))

    element_column.append(create_subtitle_div(ID_SUMMARY, "Summary"))

    def create_list(k, v):
        return "<li><strong>{0}: </strong>{1}</li>".format(k, v)

    html_listing = list(map(lambda k: create_list(k, abbreviations[k]), abbreviations.keys()))

    element_column.append(
        create_description("Table columns: <ul>{0}</ul>".format(" ".join(html_listing)), height=None))

    element_column.append(create_summary_table(df))

    element_column.append(create_subtitle_div(ID_PRECISION_VS_RECALL_TOOLS, "Average Precision vs. Average Recall"))
    element_column.append(create_description(DESCRIPTION_PRECISION_RECALL_TOOLS))
    element_column.append(create_recall_precision_scatter(df))

    element_column.append(
        create_subtitle_div(ID_RAND_INDEX_ASSIGNED_BPS, "Adjusted Rand Index vs. Percentage of Assigned base pairs"))
    element_column.append(create_description(DESCRIPTION_ADJUSTED_RAND_INDEX_TOOLS))
    element_column.append(create_rand_index_assigned_bps_scatter(df))

    element_column.append(create_subtitle_div(ID_ALL_GENOMES, "Precision vs. Recall per Genome"))
    element_column.append(create_description(DESCRIPTION_PRECISION_RECALL_GENOME))
    element_column.append(create_precision_recall_all_genomes_scatter(precision_recall_paths, names))

    save_html_file(html_output, element_column)


def main():
    parser = argparse.ArgumentParser(description="Create html-based plots.")
    parser.add_argument('-o', '--output_file', help="Directory to write the results to", required=True)
    parser.add_argument('-p', '--precision_recall_files', nargs='+', help='<Required> Set flag', required=True)
    parser.add_argument('-n', '--names', nargs='+', help='<Required> Set flag', required=True)
    parser.add_argument('-s', '--summary', help='Summary of all metrics', required=True)
    args = parser.parse_args()
    build_html(args.precision_recall_files, args.names, args.summary, args.output_file)

if __name__ == "__main__":
    main()
