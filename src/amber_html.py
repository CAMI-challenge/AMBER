import os
from collections import defaultdict
import pandas as pd
import datetime
import numpy as np
import math
import re
from jinja2 import Template
from statistics import median
from src import plots
from src.utils import labels as utils_labels
from src.utils import load_ncbi_taxinfo

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as pltx
from matplotlib.colors import rgb2hex
from matplotlib.colors import Normalize

from version import __version__

from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.models.widgets import TableColumn, Slider, Div, Select, Panel, Tabs, CheckboxGroup
from bokeh.models import (DataTable,
                          CustomJS)
from bokeh.embed import components
from bokeh.resources import INLINE
from bokeh.plotting import ColumnDataSource


TEMPLATE = Template('''<!DOCTYPE html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <title>AMBER: Assessment of Metagenome BinnERs</title>
            {{ js_resources }}
            {{ css_resources }}
            <style>.bk-fit-content {width: fit-content; width: -moz-fit-content;}
            .bk-width-auto {width: auto !important;}
            .bk-width-auto-main>div {width: -webkit-fill-available !important;}
            div.bk-width-auto-main {width: -webkit-fill-available !important;}
            .bk-tabs-margin{margin-top: 20px !important;}
            .bk-tabs-margin-lr{margin-left: 10px; margin-right: 10px}
            .bk-root {display: flex; justify-content: center;}
            .bk-padding-top {padding-top: 10px;}
            .bk-shannon-plots > div:first-child {
                padding-bottom: 0px;
                padding-left: 20px;
                padding-top: 14px;
            }
            .bk-shannon-plots > div:last-child {
                padding-top: 0px;
            }
            .bk-checkbox-shannon input[type="checkbox"] {
                margin: 2px 0 0;
            } 
            html {overflow: -moz-scrollbars-vertical; overflow-y: scroll;}
            .tooltip {
                position: relative;
                display: inline-block;
                border-bottom: 1px dashed lightgray;
                cursor: help;
            }
            .tooltip-sem {
                border-bottom: 1px dotted black;
            }
            .tooltip .tooltiptext {
                visibility: hidden;
                width: 120px;
                background-color: #555;
                color: #fff;
                text-align: center;
                border-radius: 6px;
                padding: 5px 0;
                position: absolute;
                z-index: 1;
                bottom: 125%;
                left: 50%;
                margin-left: -60px;
                opacity: 0;
            }
            .tooltip .tooltiptext::after {
                content: "";
                position: absolute;
                top: 100%;
                left: 50%;
                margin-left: -5px;
                border-width: 5px;
                border-style: solid;
                border-color: #555 transparent transparent transparent;
            }
            .tooltip:hover .tooltiptext {
                visibility: visible;
                opacity: 1;
            }
            .proportions {
                cursor: pointer;
            }
            .legend {
                position:absolute;
                cursor: move;
                z-index: 1;
            }
            </style>
        </head>
        <body>
            {{ div }}
            {{ script }}
        </body>
    </html>
    ''')


def create_title_div(id, name, info):
    div = Div(text="""<div style="text-align:left;font-size: 20pt;font-weight: bold;">{1}<span style="float: right;font-size: 10pt;font-weight:normal;">{2}</span>""".format(id, name, info), css_classes=['bk-width-auto']) # width=DIV_WIDTH, height=DIV_HEIGHT)
    return div


def create_table_html(df_summary):
    metrics1 = [utils_labels.AVG_PRECISION,
                utils_labels.AVG_PRECISION_STD,
                utils_labels.AVG_PRECISION_SEM,
                utils_labels.AVG_RECALL,
                utils_labels.AVG_RECALL_STD,
                utils_labels.AVG_RECALL_SEM]
    metrics2 = [utils_labels.AVG_PRECISION_PER_BP,
                utils_labels.AVG_RECALL_PER_BP,
                utils_labels.RI_BY_SEQ,
                utils_labels.ARI_BY_SEQ,
                utils_labels.RI_BY_BP,
                utils_labels.ARI_BY_BP,
                utils_labels.PERCENTAGE_ASSIGNED_BPS,
                utils_labels.ACCURACY,
                utils_labels.MISCLASSIFICATION]
    all_metrics = [metrics1, metrics2]
    metrics1_label = 'Quality of bins: all bins have the same weight'
    metrics2_label = 'Quality for sample: quality weighted by bin sizes'
    all_metrics_labels = [metrics1_label, metrics2_label]

    styles = [{'selector': 'td', 'props': [('width', '90pt')]},
              {'selector': 'th', 'props': [('width', '90pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '205pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    df_summary.index.name = None

    html = ''
    first_metrics = True
    for metrics, metrics_label in zip(all_metrics, all_metrics_labels):
        html += '<p style="margin-bottom: auto"><b>{}</b></p>'.format(metrics_label)
        df_metrics = df_summary.loc[metrics]
        df_metrics.index = [x[:1].upper() + x[1:] for x in df_metrics.index.values]
        sorted_columns = df_metrics.columns.tolist()
        df_metrics = df_metrics.loc[:, sorted_columns]

        if first_metrics:
            this_style = styles
            first_metrics = False
        else:
            this_style = styles_hidden_thead
        html += df_metrics.style.set_table_styles(this_style).render()

    return html


def create_genome_binning_html(df_summary):
    df_summary_g = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'genome']
    if df_summary_g.size == 0:
        return None

    colors_list = plots.create_colors_list()
    bokeh_colors = [matplotlib.colors.to_hex(c) for c in colors_list]

    p = figure(title=None, plot_width=430, plot_height=400, x_range=(0, 1), y_range=(0, 1))
    p.xaxis.axis_label = utils_labels.AVG_PRECISION[:1].upper() + utils_labels.AVG_PRECISION[1:]
    p.yaxis.axis_label = utils_labels.AVG_RECALL[:1].upper() + utils_labels.AVG_RECALL[1:]
    for color, (tool, pd_bins_tool) in zip(bokeh_colors, df_summary_g.groupby(utils_labels.TOOL)):
        source = ColumnDataSource(data=pd_bins_tool)
        p.circle(utils_labels.AVG_PRECISION, utils_labels.AVG_RECALL, source=source, legend=tool, color=color, fill_alpha=0.2, size=10)

    genome_html = create_table_html(df_summary_g.rename(columns={'tool': 'Tool'}).set_index('Tool').T)
    genome_div = Div(text="""<div>{}</div>""".format(genome_html), css_classes=['bk-width-auto'])
    metrics_row_g = row(column(genome_div, sizing_mode='scale_width', css_classes=['bk-width-auto']), column(column(p), sizing_mode='scale_width', css_classes=['bk-width-auto']), css_classes=['bk-width-auto'], sizing_mode='scale_width')

    return metrics_row_g


def create_taxonomic_binning_html(df_summary):
    df_summary_t = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'taxonomic'].rename(columns={'tool': 'Tool'})
    rank_to_html = {}
    for rank, pd_group in df_summary_t.groupby('rank'):
        pd_rank = pd_group.set_index('Tool').T
        rank_to_html[rank] = [create_table_html(pd_rank)]

    if not rank_to_html:
        return None

    taxonomic_div = Div(text="""<div>{}</div>""".format(rank_to_html[load_ncbi_taxinfo.RANKS[0]][0]), css_classes=['bk-width-auto'])
    source = ColumnDataSource(data=rank_to_html)

    select_rank = Select(title="Taxonomic rank:", value=load_ncbi_taxinfo.RANKS[0], options=load_ncbi_taxinfo.RANKS, css_classes=['bk-fit-content'])
    select_rank_callback = CustomJS(args=dict(source=source), code="""
        mytable.text = source.data[select_rank.value][0];
    """)
    select_rank.js_on_change('value', select_rank_callback)
    select_rank_callback.args["mytable"] = taxonomic_div
    select_rank_callback.args["select_rank"] = select_rank
    metrics_row_t = row(column(select_rank, taxonomic_div, sizing_mode='scale_width', css_classes=['bk-width-auto']), css_classes=['bk-width-auto'], sizing_mode='scale_width')
    return metrics_row_t


def create_html(df_summary, pd_bins, output_dir, desc_text):
    tabs_list = []

    metrics_row_g = create_genome_binning_html(df_summary)
    if metrics_row_g:
        tabs_list.append(Panel(child=metrics_row_g, title="Genome binning"))

    metrics_row_t = create_taxonomic_binning_html(df_summary)
    if metrics_row_t:
        tabs_list.append(Panel(child=metrics_row_t, title="Taxonomic binning"))

    tabs = Tabs(tabs=tabs_list, css_classes=['bk-tabs-margin'])

    title = create_title_div("main", "AMBER: Assessment of Metagenome BinnERs", " produced on {0} with AMBER version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))

    if desc_text:
        data_desc_div = Div(text="""<div style="text-align:left;font-size: 11pt;font-weight: bold;">{}""".format(desc_text), css_classes=['bk-width-auto'])
        html_columns = column(title, data_desc_div, tabs, sizing_mode='scale_width', css_classes=['bk-width-auto-main'])
    else:
        html_columns = column(title, tabs, sizing_mode='scale_width', css_classes=['bk-width-auto-main'])
    script, div = components(html_columns)
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    html = TEMPLATE.render(js_resources=js_resources,
            css_resources=css_resources,
            script=script,
            div=div)

    with open(os.path.join(output_dir, "index.html"), 'w') as f:
        f.write(html)
