#!/usr/bin/env python

import argparse
import datetime
import re
from math import pi
from collections import OrderedDict

import bokeh.models.widgets.tables
import numpy as np
import pandas as pd
from bokeh import palettes
from bokeh.embed import file_html
from bokeh.layouts import widgetbox, layout
from bokeh.models import (ColorBar,
                          Text,
                          BasicTicker,
                          HoverTool,
                          FuncTickFormatter,
                          DataTable,
                          widgets)
from bokeh.models.mappers import LinearColorMapper
from bokeh.models.widgets import Div
from bokeh.plotting import figure, ColumnDataSource
from bokeh.resources import CDN

import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src.utils import labels
from version import __version__

SCATTER_ELEMENT_WIDTH = 1500
SCATTER_ELEMENT_HEIGHT = 500

TABLE_ELEMENT_WIDTH = 750
TABLE_ELEMENT_HEIGHT = 300

DIV_WIDTH = 1200
DIV_HEIGHT = None

DESCRIPTION_WIDTH = 600
DESCRIPTION_HEIGHT = 60

A_WIDTH = 1200
A_HEIGHT = 15

COLORS_20 = palettes.d3['Category20'][20]
COLORS_10 = palettes.d3['Category10'][10]
HEATMAP_COLORS = list(reversed(palettes.RdYlBu[10])) + ['white']

NO_COLOR_COLUMNS = [labels.SEM_PRECISION, labels.SEM_RECALL, labels.STD_DEV_PRECISION, labels.STD_DEV_RECALL]

COL_DESCRIPTION = "description"
COL_TITLE = "title"
COL_FILE_ID = "file_id"

ID_SUMMARY = "summary"
ID_CONTAMINATION_COMPLETENESS = "contamination_completeness"
ID_PRECISION_VS_RECALL_TOOLS = "purity_vs_completeness_tools"
ID_PRECISION_VS_RECALL_TOOLS_BASE_PAIR = "purity_vs_completeness_tools_base_pair"
ID_RAND_INDEX_ASSIGNED_BPS = "rand_index_assigned_bps"
ID_ALL_GENOMES = "all_genomes"
ID_SINGLE_TOOL = "single_tool"

DESCRIPTION_PRECISION_RECALL_TOOLS = """
Scatter plot showing the average purity per bin vs. the average completeness per genome including the standard error
of the mean.
"""

DESCRIPTION_PRECISION_RECALL_TOOLS_BASE_PAIR = """
Scatter plot showing the average purity vs. the average completeness per base pair.
"""

DESCRIPTION_PRECISION_RECALL_GENOME = """
Scatter plot showing purity vs. completeness per genome per tool.
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
    figure.axis.axis_label_text_font_size = "15pt"
    figure.axis.major_label_text_font_size = "13pt"
    figure.xaxis.axis_label = x_label
    figure.yaxis.axis_label = y_label
    figure.xaxis.major_label_orientation = "vertical"
    return figure


def create_description(content):
    div = Div(text="""<h3>Description:</h3><p>{0}</p>""".format(content))
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


def get_contamination_completeness_thresholds(df):
    contamination_completeness_cols = [c for c in df.columns]
    completeness_thr = []
    contamination_thr = []
    for col in contamination_completeness_cols:
        values = col.split('<')
        completeness_item = values[0]
        contamination_item = '<' + values[1]
        if completeness_item not in completeness_thr:
            completeness_thr.append(completeness_item)
        if contamination_item not in contamination_thr:
            contamination_thr.append(contamination_item)
    return completeness_thr, contamination_thr


def create_contamination_completeness_table(df, completeness_thr, contamination_thr):
    contamination_completenes_row_arr = []
    for index, row in df.iterrows():
        for contamination_item in contamination_thr:
            cols = OrderedDict()
            contamination_val = format(float(re.match("\d+\.\d+", contamination_item[1:]).group()) * 100, '.0f')
            cols['Tool'] = '{} <{}% contamination'.format(row.name, contamination_val)
            for completeness_item in completeness_thr:
                cols[completeness_item] = row['{}{}'.format(completeness_item, contamination_item)]
            contamination_completenes_row_arr.append(cols)

    df = pd.DataFrame(contamination_completenes_row_arr)

    def create_table_column(field):
        if field == "Tool":
            return widgets.TableColumn(title=field, field=field, width=600)
        else:
            completeness_val = format(float(re.match("\d+\.\d+", field[1:]).group()) * 100, '.0f')
            return widgets.TableColumn(title='>{}% completeness'.format(completeness_val), field=field)

    dt = DataTable(source=bokeh.models.ColumnDataSource(df),
                   columns=list(map(lambda x: create_table_column(x), df.columns.values)),
                   width=TABLE_ELEMENT_WIDTH,
                   reorderable=True,
                   selectable=True)
    return [widgetbox(dt)]


def create_summary_heatmap(df, std_dev_sem_columns):
    df = df.set_index('Tool').iloc[::-1]
    tools = list(df.index)
    metrics = list(reversed(list(df.columns)))

    # unweighted columns should be at the right side
    for column in std_dev_sem_columns:
        metrics.append(metrics.pop(metrics.index(column)))

    df = df[metrics]

    UNWEIGHTED_NUMBER = 1.10001
    WEIGHTING_COLUMN = 'rate_extended'
    DEFAULT_TOOL_HEIGHT = 10
    COLORBAR_HEIGHT = 150
    ALPHA_COLOR=0.85

    df.columns.name = 'Metrics'
    df = pd.DataFrame(df.stack(), columns=['rate']).reset_index()
    df['rate'] = df['rate'].map('{:,.5f}'.format)
    df[WEIGHTING_COLUMN] = df['rate']

    for column in std_dev_sem_columns:
        df.loc[df.Metrics == column, WEIGHTING_COLUMN] = UNWEIGHTED_NUMBER

    mapper = LinearColorMapper(palette=HEATMAP_COLORS, low=0, high=UNWEIGHTED_NUMBER)
    source = ColumnDataSource(df)

    p = figure(x_range=metrics, y_range=tools,
               x_axis_location="above", plot_height=len(tools) * DEFAULT_TOOL_HEIGHT + COLORBAR_HEIGHT,
               tools="hover,save,box_zoom,reset,wheel_zoom", toolbar_location='below')

    p = _set_default_figure_properties(p, "Metrics", "Tools")
    p.xaxis.major_label_orientation = pi / 2.5

    p.rect(x="Metrics", y="Tool",
           width=1,
           height=1,
           source=source,
           alpha=ALPHA_COLOR,
           fill_color={'field': WEIGHTING_COLUMN, 'transform': mapper},
           line_color="black")

    glyph = Text(x="Metrics", y="Tool", text_align="center",
                 text_font_size="10pt",
                 text_baseline="middle", text="rate", text_color="black")
    p.add_glyph(source, glyph)

    tickFormatter = FuncTickFormatter(code="""
    if(tick==1){
        return tick + " (good)"
    } else if(tick==0){
        return tick + " (bad)"
    } else if(tick==1.1){
        return "Unweighted"
    } else {
        return tick.toLocaleString(
  undefined, // use a string like 'en-US' to override browser locale
  { minimumFractionDigits: 2 }
);
    }
    """)

    color_bar = ColorBar(color_mapper=mapper,
                         major_label_text_font_size="12pt",
                         ticker=BasicTicker(desired_num_ticks=len(HEATMAP_COLORS)),
                         scale_alpha=ALPHA_COLOR,
                         major_label_text_align="right",
                         major_label_text_baseline="middle",
                         bar_line_color="black",
                         formatter=tickFormatter,
                         label_standoff=13,
                         orientation="horizontal",
                         location=(-250, 0))

    p.add_layout(color_bar, 'above')
    p.select_one(HoverTool).tooltips = [
        ('Metric', '@Metrics'),
        ('Tool', '@Tool'),
        ('Value', '@rate'),
    ]
    return [p]


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


def create_scatter(df, x, y, x_title, y_title):
    """
    Creates scatter plot.
    :param df: Dataframe
    :param x: column name for x variable
    :param y: column name for y variable
    :param x_title: Title for x axis
    :param y_title: Title for y axis
    :return: bokeh figure
    """
    p = figure(plot_width=SCATTER_ELEMENT_WIDTH, plot_height=SCATTER_ELEMENT_HEIGHT)

    for i, (idx, row) in zip(np.arange(len(df.index)), df.iterrows()):
        p.circle(row[x], row[y], alpha=0.8, color=get_color(len(df.index))[i],
                 size=10,
                 legend=idx)

    p = _set_default_figure_properties(p, x_title, y_title)

    return [p]


def create_precision_recall_all_genomes_scatter(paths, names):
    """
    Creates purity vs completeness scatter plot for all tools and genomes.
    :param paths: Path to all purity and completeness tables
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
        plot.circle("purity", "completeness", color=color, alpha=0.8, legend=name, source=source)
        return plot

    for idx, name in enumerate(names):
        p = create_circles(df, p, name, get_color(len(names))[idx])

    p = _set_default_figure_properties(p, "Purity", "Completeness")
    return [p]


def save_html_file(path, elements):
    html = file_html(layout(elements, sizing_mode='scale_width'), CDN, "AMBER: Assessment of Metagenome BinnERs")
    file = open(path, "w+")
    file.write(html)
    file.close()


def build_html_summary_path(precision_recall_paths, names, summary, html_output):
    summary_df = pd.DataFrame.from_csv(summary, sep='\t', header=0)
    build_html(precision_recall_paths, names, summary_df, html_output, NO_COLOR_COLUMNS)


def build_html(precision_recall_paths, names, summary, html_output, std_dev_sem_columns=NO_COLOR_COLUMNS):
    contamination_completeness_cols = [c for c in summary.columns if c.startswith('>')]
    element_column = list()

    summary.insert(0, "Tool", summary.index)
    element_column.append(
        create_title_div("main", "AMBER: Assessment of Metagenome BinnERs", " produced on {0} with AMBER version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__)))
    element_column.append(create_subtitle_div("contents", "Contents"))

    element_column.append(create_title_a(ID_SUMMARY, "1", "Summary"))
    element_column.append(create_title_a(ID_CONTAMINATION_COMPLETENESS, "2", "Contamination and completeness"))
    element_column.append(create_title_a(ID_PRECISION_VS_RECALL_TOOLS, "3", "Average purity per bin vs. average completeness per genome"))
    element_column.append(create_title_a(ID_PRECISION_VS_RECALL_TOOLS_BASE_PAIR, "4", "Average purity vs. average completeness per base pair"))
    element_column.append(
        create_title_a(ID_RAND_INDEX_ASSIGNED_BPS, "5", "Adjusted Rand Index vs. percentage of assigned base pairs"))
    element_column.append(create_title_a(ID_ALL_GENOMES, "6", "Purity vs. completeness of all genome bins"))

    element_column.append(create_subtitle_div(ID_SUMMARY, "Summary"))

    def create_entry(k, v):
        return "<li><strong>{0}: </strong>{1}</li>".format(k, v)

    html_listing = list(map(lambda k: create_entry(k, labels.abbreviations[k]),
                            list(filter(lambda k: k not in contamination_completeness_cols,
                                        labels.abbreviations.keys()))))

    element_column.append(
        create_description("<ul>{0}</ul>".format(" ".join(html_listing))))
    df_without_contam_complete = summary.drop(contamination_completeness_cols, axis=1)

    element_column.append(create_summary_heatmap(df_without_contam_complete, std_dev_sem_columns))

    completeness_thr, contamination_thr = get_contamination_completeness_thresholds(summary[contamination_completeness_cols])

    completeness_thr_formatted = []
    for completeness_item in completeness_thr:
        value = format(float(re.match("\d+\.\d+", completeness_item[1:]).group()) * 100, '.0f')
        completeness_thr_formatted.append('{}% completeness'.format(value))

    contamination_thr_formatted = []
    for contamination_item in contamination_thr:
        value = format(float(re.match("\d+\.\d+", contamination_item[1:]).group()) * 100, '.0f')
        contamination_thr_formatted.append('{}% contamination'.format(value))

    def create_entry_completeness(v):
        return "<li><strong>>{0}: </strong>number of bins with more than {1}</li>".format(v, v)

    def create_entry_contamination(v):
        return "<li><strong><{0}: </strong>number of bins with less than {1}</li>".format(v, v)

    completion_contamination_col = [create_entry_completeness(v) for v in completeness_thr_formatted]
    completion_contamination_row = [create_entry_contamination(v) for v in contamination_thr_formatted]

    element_column.append(create_subtitle_div(ID_CONTAMINATION_COMPLETENESS, "Contamination and completeness"))

    element_column.append(
        create_description("Table columns: <ul>{0}</ul> "
                           "Table rows: <ul>{1}</ul>".format(" ".join(completion_contamination_col),
                                                             " ".join(completion_contamination_row))))
    element_column.append(create_contamination_completeness_table(summary[contamination_completeness_cols], completeness_thr, contamination_thr))

    element_column.append(create_subtitle_div(ID_PRECISION_VS_RECALL_TOOLS, "Average purity per bin vs. average completeness per genome"))
    element_column.append(create_description(DESCRIPTION_PRECISION_RECALL_TOOLS))
    element_column.append(create_scatter(summary, "avg_purity", "avg_completeness", "Average purity per bin", "Average completeness per genome"))

    element_column.append(create_subtitle_div(ID_PRECISION_VS_RECALL_TOOLS_BASE_PAIR, "Average purity per base pair vs. average completeness per base pair"))
    element_column.append(create_description(DESCRIPTION_PRECISION_RECALL_TOOLS_BASE_PAIR))
    element_column.append(create_scatter(summary, "avg_purity_per_bp", "avg_completeness_per_bp", "Average purity per base pair", "Average completeness per base pair"))

    element_column.append(
        create_subtitle_div(ID_RAND_INDEX_ASSIGNED_BPS, "Adjusted Rand Index vs. percentage of assigned base pairs"))
    element_column.append(create_description(DESCRIPTION_ADJUSTED_RAND_INDEX_TOOLS))
    element_column.append(create_scatter(summary, "a_rand_index_by_bp", "percent_assigned_bps",
                                         "Adjusted Rand Index by base pair", "Percentage of assigned base pairs"))

    element_column.append(create_subtitle_div(ID_ALL_GENOMES, "Purity vs. completeness of all genome bins"))
    element_column.append(create_description(DESCRIPTION_PRECISION_RECALL_GENOME))
    element_column.append(create_precision_recall_all_genomes_scatter(precision_recall_paths, names))

    save_html_file(html_output, element_column)


def main():
    parser = argparse.ArgumentParser(description="Create HTML-based plots.")
    parser.add_argument('-o', '--output_file', help="Directory to write the results to", required=True)
    parser.add_argument('-p', '--purity_completeness_files', nargs='+', help='<Required> Set flag', required=True)
    parser.add_argument('-n', '--names', nargs='+', help='<Required> Set flag', required=True)
    parser.add_argument('-s', '--summary', help='Summary of all metrics', required=True)
    args = parser.parse_args()
    build_html_summary_path(args.purity_completeness_files, args.names, args.summary, args.output_file)


if __name__ == "__main__":
    main()
