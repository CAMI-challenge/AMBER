import os
from collections import OrderedDict
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
from bokeh.layouts import column, row, widgetbox
from bokeh.models.widgets import TableColumn, Div, Select, Panel, Tabs, CheckboxGroup
from bokeh.models import (DataTable,
                          CustomJS,
                          Legend,
                          Band,
                          FuncTickFormatter,
                          HoverTool)
from bokeh.models.tickers import FixedTicker
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
            .bk-width-auto {width: auto !important; height: auto !important;}
            .bk-display-block {display: block !important;}
            .bk-float-left {float: left;}
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


def create_heatmap_bar(output_dir):
    fig = pltx.figure(figsize=(2, 2))
    x = np.arange(25).reshape(5, 5)
    cm = sns.diverging_palette(12, 250, sep=30, l=50, n=15, s=90, as_cmap=True)
    ax = sns.heatmap(x, cmap=cm, cbar=False)
    cbar = pltx.colorbar(ax.get_children()[0], orientation='horizontal')

    # hide border
    cbar.outline.set_visible(False)

    # hide ticks
    cbar.set_ticks([])

    # hide heatmap
    pltx.gca().set_visible(False)

    fig.savefig(os.path.join(output_dir, 'heatmap_bar.png'), dpi=100, format='png', bbox_inches='tight', pad_inches=-.035, transparent=True)
    pltx.close(fig)


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def upper1(x):
        return x[:1].upper() + x[1:]


def get_colors_and_ranges(name):
    color1 = 'dodgerblue'
    color2 = 'red'
    hue1 = 12
    hue2 = 240

    if name == upper1(utils_labels.MISCLASSIFICATION_PER_BP) or name == upper1(utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN) or name == upper1(utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN):
        return color2, color1, hue2, hue1, 0, 1
    return color1, color2, hue1, hue2, 0, 1


def get_heatmap_colors(pd_series, **args):
    values = pd_series.tolist()
    notnan_values = [x for x in values if isinstance(x, (float, int)) and not np.isnan(x)]

    if pd_series.name == upper1(utils_labels.AVG_PRECISION_BP_SEM) or pd_series.name == upper1(utils_labels.AVG_RECALL_BP_SEM) or\
        pd_series.name == upper1(utils_labels.AVG_PRECISION_SEQ_SEM) or pd_series.name == upper1(utils_labels.AVG_RECALL_SEQ_SEM):
        return ['background-color: white' for x in values]

    color1, color2, hue1, hue2, min_value, max_value = get_colors_and_ranges(pd_series.name)

    cm = sns.diverging_palette(hue1, hue2, sep=50, l=80, n=15, s=90, as_cmap=True)
    norm = MidpointNormalize(vmin=min_value, vmax=max_value, midpoint=median(notnan_values))

    normed = norm([round(x, 3) if not math.isnan(x) else max_value for x in values])
    heatmap_colors = [rgb2hex(x) for x in pltx.cm.get_cmap(cm)(normed)]
    return_colors = []
    for val, x in zip(values, heatmap_colors):
        if val <= min_value:
            return_colors.append('background-color: {}'. format(color2))
        elif val >= max_value:
            return_colors.append('background-color: {}'. format(color1))
        elif math.isnan(val):
            return_colors.append('background-color: red')
        else:
            return_colors.append('background-color: {}'. format(x))
    return return_colors


def create_title_div(id, name, info):
    div = Div(text="""<div style="text-align:left;font-size: 20pt;font-weight: bold;">{1}<span style="float: right;font-size: 10pt;font-weight:normal;">{2}</span>""".format(id, name, info), css_classes=['bk-width-auto']) # width=DIV_WIDTH, height=DIV_HEIGHT)
    return div


def create_metrics_per_bin_panel(pd_bins, bins_columns):
    styles = [{'selector': 'td', 'props': [('width', '80pt')]},
              {'selector': 'th', 'props': [('width', '80pt'), ('text-align', 'left')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]

    tools = pd_bins.tool.unique().tolist()
    tool_to_html = OrderedDict([(x, '') for x in tools])

    for tool, pd_group in pd_bins.groupby(utils_labels.TOOL):
        table = pd_group[list(bins_columns.keys())].rename(columns=dict(bins_columns))
        tool_to_html[tool] = [table.style.set_table_styles(styles).hide_index().render()]

    table_div = Div(text="""<div>{}</div>""".format(tool_to_html[tools[0]][0]), css_classes=['bk-width-auto'])

    source = ColumnDataSource(data=tool_to_html)

    select_tool = Select(title="Binner:", value=tools[0], options=tools, css_classes=['bk-fit-content'])
    select_tool_callback = CustomJS(args=dict(source=source), code="""
        mytable.text = source.data[select_tool.value][0];
    """)
    select_tool.js_on_change('value', select_tool_callback)
    select_tool_callback.args["mytable"] = table_div
    select_tool_callback.args["select_tool"] = select_tool

    table_column = row(column(select_tool, table_div, sizing_mode='scale_width', css_classes=['bk-width-auto']), css_classes=['bk-width-auto'], sizing_mode='scale_width') # column(table_div, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-width-auto-main'])
    metrics_bins_panel = Panel(child=table_column, title="Metrics per bin")
    return metrics_bins_panel


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


def create_contamination_completeness_table(gold_standard, df):
    contamination_completeness_cols = [c for c in df.columns if c.startswith('>')]
    completeness_thr, contamination_thr = get_contamination_completeness_thresholds(df[contamination_completeness_cols])

    contamination_completenes_row_arr = []

    cols = OrderedDict()
    cols['Tool'] = 'Gold standard'
    # num_genomes = len(gold_standard.genome_query.get_bin_ids()) # Number of unfiltered genomes
    gs_genome_bins_metrics = gold_standard.genome_query.get_bins_metrics(gold_standard)
    num_genomes = len(gs_genome_bins_metrics)
    for completeness_item in completeness_thr:
        cols[completeness_item] = num_genomes
    contamination_completenes_row_arr.append(cols)

    for index, row in df.iterrows():
        for contamination_item in contamination_thr:
            cols = OrderedDict()
            contamination_val = format(float(re.match("\d+\.\d+", contamination_item[1:]).group()) * 100, '.0f')
            cols['Tool'] = '{} <{}% contamination'.format(row[utils_labels.TOOL], contamination_val)
            for completeness_item in completeness_thr:
                cols[completeness_item] = row['{}{}'.format(completeness_item, contamination_item)]
            contamination_completenes_row_arr.append(cols)

    df = pd.DataFrame(contamination_completenes_row_arr)

    def create_table_column(field):
        if field == "Tool":
            return TableColumn(title=field, field=field, width=600)
        else:
            completeness_val = format(float(re.match("\d+\.\d+", field[1:]).group()) * 100, '.0f')
            return TableColumn(title='>{}% completeness'.format(completeness_val), field=field)

    dt = DataTable(source=ColumnDataSource(df),
                   columns=list(map(lambda x: create_table_column(x), df.columns.values)),
                   width=750,
                   height=1000,
                   reorderable=True,
                   selectable=True)
    return [widgetbox(dt)]


def create_heatmap_div():
    heatmap_legend = '<img src="heatmap_bar.png" /><div style="text-align:left;font-size: 11px;">Worst<span style="float:right;">Best</span><span style="margin-right: 36px;float:right;">Median</span></div>'
    heatmap_legend_div = Div(text=heatmap_legend, style={"width": "155px", "margin-bottom": "-10px"})
    return  heatmap_legend_div


def create_table_html(df_summary, is_genome):
    metrics1 = [utils_labels.AVG_PRECISION_BP,
                utils_labels.AVG_PRECISION_SEQ,
                utils_labels.AVG_RECALL_BP,
                utils_labels.AVG_RECALL_SEQ,
                utils_labels.AVG_PRECISION_BP_SEM,
                utils_labels.AVG_PRECISION_SEQ_SEM,
                utils_labels.AVG_RECALL_BP_SEM,
                utils_labels.AVG_RECALL_SEQ_SEM]
    metrics2 = [utils_labels.ACCURACY_PER_BP,
                utils_labels.ACCURACY_PER_SEQ,
                utils_labels.MISCLASSIFICATION_PER_BP,
                utils_labels.MISCLASSIFICATION_PER_SEQ,
                utils_labels.AVG_PRECISION_PER_BP,
                utils_labels.AVG_PRECISION_PER_SEQ,
                utils_labels.AVG_RECALL_PER_BP,
                utils_labels.AVG_RECALL_PER_SEQ,
                utils_labels.RI_BY_SEQ,
                utils_labels.ARI_BY_SEQ,
                utils_labels.RI_BY_BP,
                utils_labels.ARI_BY_BP,
                utils_labels.PERCENTAGE_ASSIGNED_BPS,
                utils_labels.PERCENTAGE_ASSIGNED_SEQS,
                utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN,
                utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN]
    if is_genome:
        metrics2.remove(utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN)
        metrics2.remove(utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN)
    all_metrics = [metrics1, metrics2]
    metrics1_label = 'Quality of bins: all bins have the same weight'
    metrics2_label = 'Quality for sample: quality weighted by bin sizes'
    all_metrics_labels = [metrics1_label, metrics2_label]

    styles = [{'selector': 'td', 'props': [('width', '95pt')]},
              {'selector': 'th', 'props': [('width', '95pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '205pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    d = {utils_labels.ACCURACY_PER_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.ACCURACY_PER_BP, utils_labels.ACCURACY_PER_BP, utils_labels.ACCURACY_PER_BP),
         utils_labels.MISCLASSIFICATION_PER_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.MISCLASSIFICATION_PER_BP, utils_labels.MISCLASSIFICATION_PER_BP, utils_labels.MISCLASSIFICATION_PER_BP),
         utils_labels.AVG_PRECISION_PER_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.AVG_PRECISION_PER_BP, utils_labels.AVG_PRECISION_PER_BP, utils_labels.AVG_PRECISION_PER_BP),
         utils_labels.AVG_RECALL_PER_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.AVG_RECALL_PER_BP, utils_labels.AVG_RECALL_PER_BP, utils_labels.AVG_RECALL_PER_BP),
         utils_labels.RI_BY_SEQ: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.RI_BY_SEQ, utils_labels.RI_BY_SEQ, utils_labels.RI_BY_SEQ),
         utils_labels.ARI_BY_SEQ: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.ARI_BY_SEQ, utils_labels.ARI_BY_SEQ, utils_labels.ARI_BY_SEQ),
         utils_labels.RI_BY_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.RI_BY_BP, utils_labels.RI_BY_BP, utils_labels.RI_BY_BP),
         utils_labels.ARI_BY_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.ARI_BY_BP, utils_labels.ARI_BY_BP, utils_labels.ARI_BY_BP),
         utils_labels.PERCENTAGE_ASSIGNED_BPS: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.PERCENTAGE_ASSIGNED_BPS, utils_labels.PERCENTAGE_ASSIGNED_BPS, utils_labels.PERCENTAGE_ASSIGNED_BPS),
         utils_labels.PERCENTAGE_ASSIGNED_SEQS: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.PERCENTAGE_ASSIGNED_SEQS, utils_labels.PERCENTAGE_ASSIGNED_SEQS, utils_labels.PERCENTAGE_ASSIGNED_SEQS),
         utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN, utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN, utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN),
         utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN, utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN, utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN),
         utils_labels.AVG_PRECISION_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.AVG_PRECISION_BP, utils_labels.AVG_PRECISION_BP, utils_labels.AVG_PRECISION_BP),
         utils_labels.AVG_PRECISION_BP_SEM: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.AVG_PRECISION_BP_SEM, utils_labels.AVG_PRECISION_BP_SEM, utils_labels.AVG_PRECISION_BP_SEM),
         utils_labels.AVG_RECALL_BP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.AVG_RECALL_BP, utils_labels.AVG_RECALL_BP, utils_labels.AVG_RECALL_BP),
         utils_labels.AVG_RECALL_BP_SEM: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(utils_labels.AVG_RECALL_BP_SEM, utils_labels.AVG_RECALL_BP_SEM, utils_labels.AVG_RECALL_BP_SEM)}
    pattern = re.compile('|'.join(map(re.escape, d)))

    def translate(match):
        return d[match.group(0)]

    df_summary.index.name = None

    html = ''
    first_metrics = True
    for metrics, metrics_label in zip(all_metrics, all_metrics_labels):
        html += '<p style="margin-bottom: auto"><b>{}</b></p>'.format(metrics_label)
        df_metrics = df_summary.loc[metrics]
        sorted_columns = df_metrics.columns.tolist()
        df_metrics = df_metrics.loc[:, sorted_columns].fillna(0)

        if first_metrics:
            this_style = styles
            first_metrics = False
        else:
            this_style = styles_hidden_thead
        html += df_metrics.style.apply(get_heatmap_colors, df_metrics=df_metrics, axis=1).set_precision(3).set_table_styles(this_style).render()
    html = pattern.sub(translate, html)

    return '<div style="margin-bottom:10pt;">{}</div>'.format(html)


def create_precision_recall_figure(df_summary, xname, yname, title):
    colors_list = plots.create_colors_list()
    bokeh_colors = [matplotlib.colors.to_hex(c) for c in colors_list]

    legend_it = []
    p = figure(title=title, plot_width=580, plot_height=400, x_range=(0, 1), y_range=(0, 1), toolbar_location="below")
    p.xaxis.axis_label = upper1(xname)
    p.yaxis.axis_label = upper1(yname)
    for color, (index, row) in zip(bokeh_colors, df_summary.iterrows()):
        tool = row[utils_labels.TOOL]
        source = ColumnDataSource(data=row.to_frame().T)
        pcircle = p.circle(xname, yname, source=source, color=color, fill_alpha=0.2, size=10)
        legend_it.append((tool, [pcircle]))
    p.add_layout(Legend(items=legend_it), 'right')
    p.legend.click_policy = 'hide'
    return p


def create_precision_recall_all_genomes_scatter(pd_bins, tools):
    colors_list = plots.create_colors_list()
    bokeh_colors = [matplotlib.colors.to_hex(c) for c in colors_list]

    p = figure(plot_width=580, plot_height=400, x_range=(0, 1), y_range=(0, 1), toolbar_location="below")
    p.add_tools(HoverTool(tooltips=[
        ("genome", "@mapping_id")
    ]))

    pd_genome_bins = pd_bins[pd_bins['rank'] == 'NA']
    legend_it = []
    for color, tool in zip(bokeh_colors, tools):
        source = ColumnDataSource(data=pd_genome_bins[pd_genome_bins[utils_labels.TOOL] == tool])
        pcircle = p.circle('purity_bp', 'completeness_bp', color=color, alpha=0.8, source=source)
        legend_it.append((tool, [pcircle]))
    p.add_layout(Legend(items=legend_it), 'right')
    p.xaxis.axis_label = 'Purity per bin'
    p.yaxis.axis_label = 'Completeness per bin'
    p.legend.click_policy = 'hide'
    return p


def create_tax_figure(tool, df_summary, metrics_list, errors_list):
    rank_indices = list(range(len(load_ncbi_taxinfo.RANKS)))

    df = df_summary.set_index("rank").loc[load_ncbi_taxinfo.RANKS]
    df["x"] = rank_indices
    for metric, error in zip(metrics_list, errors_list):
        if error:
            df[metric + "upper"] = df[metric] + df[error]
            df[metric + "lower"] = df[metric] - df[error]
    source = ColumnDataSource(df.reset_index())

    p = figure(plot_width=500, plot_height=550, x_range=(0, 7), y_range=(0, 1))
    plines = []
    line_colors = ["#006cba", "#008000", "#ba9e00"]
    for metric, color in zip(metrics_list, line_colors):
        plines.append([p.line(x='x', y=metric, line_color=color, source=source, line_width=2)])

    p.xaxis.ticker = FixedTicker(ticks=rank_indices)
    p.xaxis.formatter = FuncTickFormatter(code="""
        var mapping = {0: "superkingdom", 1: "phylum", 2: "class", 3: "order", 4: "family", 5: "genus", 6: "species", 7: "strain"};
        return mapping[tick];
    """)
    p.xaxis.major_label_orientation = 1.2

    for metric, error, color in zip(metrics_list, errors_list, line_colors):
        if error:
            p.add_layout(Band(base='x', lower=metric + "lower", upper=metric + "upper", source=source, level='underlay',
                              fill_alpha=.3, line_width=1, line_color='black', fill_color=color))
    legend = Legend(items=list(zip(metrics_list, plines)))
    p.add_layout(legend, 'above')

    p.xgrid[0].grid_line_color=None
    p.ygrid[0].grid_line_alpha=0.5
    p.title.text = tool
    return p


def create_rankings_table(df_summary, show_rank=False):
    columns= [utils_labels.TOOL,
              utils_labels.AVG_PRECISION_BP,
              utils_labels.AVG_RECALL_BP,
              utils_labels.AVG_PRECISION_PER_BP,
              utils_labels.AVG_RECALL_PER_BP,
              utils_labels.ARI_BY_SEQ,
              utils_labels.ARI_BY_BP,
              utils_labels.PERCENTAGE_ASSIGNED_BPS,
              utils_labels.ACCURACY_PER_BP]
    if show_rank:
        columns.insert(1, utils_labels.RANK)
    pd_rankings = df_summary[columns].rename(columns={utils_labels.RANK: 'Taxonomic rank'}).round(decimals=5)
    pd_rankings.columns = list(map(upper1, pd_rankings.columns))

    def create_table_column(field):
        return TableColumn(title=field, field=field, width=100)

    dt = DataTable(source=ColumnDataSource(pd_rankings),
                   columns=list(map(lambda x: create_table_column(x), pd_rankings.columns.values)),
                   width=1500,
                   height=1000,
                   reorderable=True,
                   selectable=True)
    return [widgetbox(dt)]


def create_genome_binning_html(gold_standard, df_summary, pd_bins):
    df_summary_g = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'genome']
    if df_summary_g.size == 0:
        return None

    purity_completeness_plot = create_precision_recall_figure(df_summary_g, utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP, None)
    purity_completeness_bp_plot = create_precision_recall_figure(df_summary_g, utils_labels.AVG_PRECISION_PER_BP, utils_labels.AVG_RECALL_PER_BP, None)
    all_genomes_plot = create_precision_recall_all_genomes_scatter(pd_bins, df_summary_g[utils_labels.TOOL].tolist())

    genome_html = create_table_html(df_summary_g.rename(columns={'tool': 'Tool'}).set_index('Tool').T, True)
    genome_div = Div(text="""<div style="margin-bottom:10pt;">{}</div>""".format(genome_html), css_classes=['bk-width-auto'])
    metrics_row_g = column(column(create_heatmap_div(), genome_div, sizing_mode='scale_width', css_classes=['bk-width-auto']), column(row(purity_completeness_plot, purity_completeness_bp_plot, all_genomes_plot), sizing_mode='scale_width', css_classes=['bk-width-auto']), css_classes=['bk-width-auto'], sizing_mode='scale_width')
    metrics_panel = Panel(child=metrics_row_g, title="Metrics")

    bins_columns = OrderedDict([('id', 'Bin ID'), ('mapping_id', 'Mapped genome'), ('purity_bp', 'Purity'), ('completeness_bp', 'Completeness'), ('predicted_size', 'Predicted size'), ('true_positive_bps', 'True positives'), ('true_size', 'True size')])
    metrics_bins_panel = create_metrics_per_bin_panel(pd_bins[pd_bins['rank'] == 'NA'], bins_columns)

    cc_table = create_contamination_completeness_table(gold_standard, df_summary_g)
    cc_panel = Panel(child=row(cc_table), title="#Recovered genomes")

    rankings_panel = Panel(child=column([Div(text="Click on the columns header for sorting.", style={"width": "500px", "margin-top": "20px"}),
                                        row(create_rankings_table(df_summary_g))]), title="Rankings")

    tabs = Tabs(tabs=[metrics_panel, metrics_bins_panel, rankings_panel, cc_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])

    return tabs


def create_plots_per_binner(df_summary_t):
    tools = df_summary_t.tool.unique().tolist()

    metrics_list = [utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP]
    errors_list = [utils_labels.AVG_PRECISION_BP_SEM, utils_labels.AVG_RECALL_BP_SEM]
    tools_figures = [create_tax_figure(tool, df_summary_t[df_summary_t[utils_labels.TOOL] == tool], metrics_list, errors_list) for tool in tools]
    tools_figures_columns = [column(x, css_classes=['bk-width-auto', 'bk-float-left']) for x in tools_figures]
    tools_column_unweighted = column(tools_figures_columns, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-display-block'])

    metrics_list = [utils_labels.AVG_PRECISION_PER_BP, utils_labels.AVG_RECALL_PER_BP, utils_labels.ACCURACY_PER_BP]
    errors_list = ["", "", ""]
    tools_figures_weighted = [create_tax_figure(tool, df_summary_t[df_summary_t[utils_labels.TOOL] == tool], metrics_list, errors_list) for tool in tools]
    tools_figures_columns = [column(x, css_classes=['bk-width-auto', 'bk-float-left']) for x in tools_figures_weighted]
    tools_column_weighted = column(tools_figures_columns, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-display-block'])

    tools_unweighted_panel = Panel(child=tools_column_unweighted, title="Unweighted by bin sizes")
    tools_weighted_panel = Panel(child=tools_column_weighted, title="Weighted by bin sizes")
    tools_tabs = Tabs(tabs=[tools_unweighted_panel, tools_weighted_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])
    return tools_tabs


def create_taxonomic_binning_html(df_summary, pd_bins):
    df_summary_t = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'taxonomic']
    rank_to_html = {}
    plots_list = []
    pd_groups_rank = df_summary_t.groupby('rank')
    available_ranks = list(pd_groups_rank.groups.keys())
    available_ranks_sorted = [rank for rank in load_ncbi_taxinfo.RANKS if rank in available_ranks]

    for rank in load_ncbi_taxinfo.RANKS:
        if rank not in available_ranks:
            continue
        pd_group = pd_groups_rank.get_group(rank)
        pd_rank = pd_group.rename(columns={'tool': 'Tool'}).set_index('Tool').T
        rank_to_html[rank] = [create_table_html(pd_rank, False)]
        purity_completeness_plot = create_precision_recall_figure(pd_group, utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP, rank)
        purity_completeness_bp_plot = create_precision_recall_figure(pd_group, utils_labels.AVG_PRECISION_PER_BP, utils_labels.AVG_RECALL_PER_BP, rank)
        plots_list.append(row(purity_completeness_plot, purity_completeness_bp_plot))

    if not rank_to_html:
        return None

    taxonomic_div = Div(text="""<div style="margin-bottom:10pt;">{}</div>""".format(rank_to_html[load_ncbi_taxinfo.RANKS[0]][0]), css_classes=['bk-width-auto'])
    source = ColumnDataSource(data=rank_to_html)

    select_rank = Select(title="Taxonomic rank:", value=load_ncbi_taxinfo.RANKS[0], options=available_ranks_sorted, css_classes=['bk-fit-content'])
    select_rank_callback = CustomJS(args=dict(source=source), code="""
        mytable.text = source.data[select_rank.value][0];
    """)
    select_rank.js_on_change('value', select_rank_callback)
    select_rank_callback.args["mytable"] = taxonomic_div
    select_rank_callback.args["select_rank"] = select_rank

    metrics_column = column(column(select_rank, create_heatmap_div(), taxonomic_div, sizing_mode='scale_width', css_classes=['bk-width-auto']), column(plots_list, sizing_mode='scale_width', css_classes=['bk-width-auto']), css_classes=['bk-width-auto'], sizing_mode='scale_width')
    metrics_panel = Panel(child=metrics_column, title="Metrics")

    bins_columns = OrderedDict([('mapping_id', 'Taxon ID'), ('name', 'Scientific name'), ('rank', 'Taxonomic rank'), ('purity_bp', 'Purity'), ('completeness_bp', 'Completeness'), ('predicted_size', 'Predicted size'), ('true_positive_bps', 'True positives'), ('true_size', 'True size')])
    tax_bins = pd_bins[pd_bins['rank'] != 'NA']
    if tax_bins['name'].isnull().any():
        del bins_columns['name']
    metrics_bins_panel = create_metrics_per_bin_panel(tax_bins, bins_columns)

    rankings_panel = Panel(child=column([Div(text="Click on the columns header for sorting.", style={"width": "500px", "margin-top": "20px"}),
                                        row(create_rankings_table(df_summary_t, True))]), title="Rankings")

    tools_panel = Panel(child=create_plots_per_binner(df_summary_t), title="Plots per binner")

    tabs = Tabs(tabs=[metrics_panel, metrics_bins_panel, rankings_panel, tools_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])

    return tabs


def create_html(gold_standard, df_summary, pd_bins, output_dir, desc_text):
    create_heatmap_bar(output_dir)
    tabs_list = []

    metrics_row_g = create_genome_binning_html(gold_standard, df_summary, pd_bins)
    if metrics_row_g:
        tabs_list.append(Panel(child=metrics_row_g, title="Genome binning"))

    metrics_row_t = create_taxonomic_binning_html(df_summary, pd_bins)
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
