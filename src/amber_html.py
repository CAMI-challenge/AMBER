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
                          FuncTickFormatter)
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

    if name == upper1(utils_labels.MISCLASSIFICATION):
        return color2, color1, hue2, hue1, 0, 1
    return color1, color2, hue1, hue2, 0, 1


def get_heatmap_colors(pd_series, **args):
    values = pd_series.tolist()
    notnan_values = [x for x in values if isinstance(x, (float, int)) and not np.isnan(x)]

    if pd_series.name == upper1(utils_labels.AVG_PRECISION_SEM) or pd_series.name == upper1(utils_labels.AVG_PRECISION_STD) or\
        pd_series.name == upper1(utils_labels.AVG_RECALL_SEM) or pd_series.name == upper1(utils_labels.AVG_RECALL_STD):
        return ['background-color: white' for x in values]

    color1, color2, hue1, hue2, min_value, max_value = get_colors_and_ranges(pd_series.name)

    cm = sns.diverging_palette(hue1, hue2, sep=50, l=80, n=15, s=90, as_cmap=True)
    norm = MidpointNormalize(vmin=min_value, vmax=max_value, midpoint=median(notnan_values))

    normed = norm([round(x, 3) if not math.isnan(x) else max_value for x in values])
    heatmap_colors = [rgb2hex(x) for x in pltx.cm.get_cmap(cm)(normed)]
    return_colors = []
    for val, x in zip(values, heatmap_colors):
        if val == min_value:
            return_colors.append('background-color: {}'. format(color2))
        elif val == max_value:
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

    pd_groups = pd_bins.groupby(utils_labels.TOOL)
    tools = list(pd_groups.groups.keys())
    tool_to_html = OrderedDict([(x, '') for x in tools])

    for tool, pd_group in pd_groups:
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


def create_contamination_completeness_table(df):
    contamination_completeness_cols = [c for c in df.columns if c.startswith('>')]
    completeness_thr, contamination_thr = get_contamination_completeness_thresholds(df[contamination_completeness_cols])

    contamination_completenes_row_arr = []
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
                   reorderable=True,
                   selectable=True)
    return [widgetbox(dt)]


def create_heatmap_div():
    heatmap_legend = '<img src="heatmap_bar.png" /><div style="text-align:left;font-size: 11px;">Worst<span style="float:right;">Best</span><span style="margin-right: 36px;float:right;">Median</span></div>'
    heatmap_legend_div = Div(text=heatmap_legend, style={"width": "155px", "margin-bottom": "-10px"})
    return  heatmap_legend_div


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

    styles = [{'selector': 'td', 'props': [('width', '100pt')]},
              {'selector': 'th', 'props': [('width', '100pt'), ('text-align', 'left')]},
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
        html += df_metrics.style.apply(get_heatmap_colors, df_metrics=df_metrics, axis=1).set_precision(3).set_table_styles(this_style).render()

    return html


def create_figure(df_summary, xname, yname):
    colors_list = plots.create_colors_list()
    bokeh_colors = [matplotlib.colors.to_hex(c) for c in colors_list]

    legend_it = []
    fig = figure(title=None, plot_width=450, plot_height=500, x_range=(0, 1), y_range=(0, 1), toolbar_location="right")
    fig.xaxis.axis_label = upper1(xname)
    fig.yaxis.axis_label = upper1(yname)
    for color, (tool, pd_bins_tool) in zip(bokeh_colors, df_summary.groupby(utils_labels.TOOL)):
        source = ColumnDataSource(data=pd_bins_tool)
        p = fig.circle(xname, yname, source=source, color=color, fill_alpha=0.2, size=10)
        legend_it.append((tool, [p]))
    legend = Legend(items=legend_it)
    fig.add_layout(legend, 'above')
    return fig


def create_tax_figure(tool, df_summary):
    rank_indices = list(range(len(load_ncbi_taxinfo.RANKS)))

    rank_to_precision = OrderedDict([(k, .0) for k in load_ncbi_taxinfo.RANKS])
    rank_to_precision_error = OrderedDict([(k, .0) for k in load_ncbi_taxinfo.RANKS])
    rank_to_recall = OrderedDict([(k, .0) for k in load_ncbi_taxinfo.RANKS])
    rank_to_recall_error = OrderedDict([(k, .0) for k in load_ncbi_taxinfo.RANKS])
    for index, row in df_summary.iterrows():
        rank_to_precision[row[utils_labels.RANK]] = .0 if np.isnan(row[utils_labels.AVG_PRECISION]) else row[utils_labels.AVG_PRECISION]
        rank_to_precision_error[row[utils_labels.RANK]] = .0 if np.isnan(row[utils_labels.AVG_PRECISION_SEM]) else row[utils_labels.AVG_PRECISION_SEM]
        rank_to_recall[row[utils_labels.RANK]] = .0 if np.isnan(row[utils_labels.AVG_RECALL]) else row[utils_labels.AVG_RECALL]
        rank_to_recall_error[row[utils_labels.RANK]] = .0 if np.isnan(row[utils_labels.AVG_RECALL_SEM]) else row[utils_labels.AVG_RECALL_SEM]

    y_values1 = list(rank_to_precision.values())
    y_values2 = list(rank_to_recall.values())
    sem1 = list(rank_to_precision_error.values())
    sem2 = list(rank_to_recall_error.values())
    lower1 = np.subtract(y_values1, sem1).tolist()
    lower2 = np.subtract(y_values2, sem2).tolist()
    upper1 = np.add(y_values1, sem1).tolist()
    upper2 = np.add(y_values2, sem2).tolist()

    df = pd.DataFrame({'x': rank_indices, 'y1': y_values1, 'y2': y_values2, 'lower1': lower1, 'upper1': upper1, 'lower2': lower2, 'upper2': upper2})
    source = ColumnDataSource(df.reset_index())

    p = figure(plot_width=450, plot_height=500)
    pline1 = p.line(x='x', y='y1', line_color='#006cba', source=source, line_width=2)
    pline2 = p.line(x='x', y='y2', line_color='#008000', source=source, line_width=2)
    p.xaxis.ticker = FixedTicker(ticks=rank_indices)
    p.xaxis.formatter = FuncTickFormatter(code="""
        var mapping = {0: "superkingdom", 1: "phylum", 2: "class", 3: "order", 4: "family", 5: "genus", 6: "species", 7: "strain"};
        return mapping[tick];
    """)
    band1 = Band(base='x', lower='lower1', upper='upper1', source=source, level='underlay',
                fill_alpha=.3, line_width=1, line_color='black', fill_color='#006cba')
    band2 = Band(base='x', lower='lower2', upper='upper2', source=source, level='underlay',
                fill_alpha=.3, line_width=1, line_color='black', fill_color='#008000')
    p.add_layout(band1)
    p.add_layout(band2)
    legend_it = [('Average purity', [pline1]), ('Average completeness', [pline2])]
    legend = Legend(items=legend_it)
    p.add_layout(legend, 'above')

    p.xgrid[0].grid_line_color=None
    p.ygrid[0].grid_line_alpha=0.5
    p.title.text = tool
    return p


def create_genome_binning_html(df_summary, pd_bins):
    df_summary_g = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'genome']
    if df_summary_g.size == 0:
        return None

    purity_completeness_plot = create_figure(df_summary_g, utils_labels.AVG_PRECISION, utils_labels.AVG_RECALL)
    purity_completeness_bp_plot = create_figure(df_summary_g, utils_labels.AVG_PRECISION_PER_BP, utils_labels.AVG_RECALL_PER_BP)

    genome_html = create_table_html(df_summary_g.rename(columns={'tool': 'Tool'}).set_index('Tool').T)
    genome_div = Div(text="""<div>{}</div>""".format(genome_html), css_classes=['bk-width-auto'])
    metrics_row_g = row(column(create_heatmap_div(), genome_div, sizing_mode='scale_width', css_classes=['bk-width-auto']), column(row(purity_completeness_plot, purity_completeness_bp_plot), sizing_mode='scale_width', css_classes=['bk-width-auto']), css_classes=['bk-width-auto'], sizing_mode='scale_width')
    metrics_panel = Panel(child=metrics_row_g, title="Metrics")

    bins_columns = OrderedDict([('id', 'Bin ID'), ('mapping_id', 'Mapped genome'), ('purity', 'Purity'), ('completeness', 'Completeness'), ('predicted_size', 'Predicted size'), ('true_positives', 'True positives'), ('real_size', 'Real size')])
    metrics_bins_panel = create_metrics_per_bin_panel(pd_bins[pd_bins['rank'] == 'NA'], bins_columns)

    cc_table = create_contamination_completeness_table(df_summary_g)
    cc_panel = Panel(child=row(cc_table), title="#Recovered bins")

    tabs = Tabs(tabs=[metrics_panel, metrics_bins_panel, cc_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])

    return tabs


def create_taxonomic_binning_html(df_summary, pd_bins):
    df_summary_t = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'taxonomic'].rename(columns={'tool': 'Tool'})
    rank_to_html = {}
    for rank, pd_group in df_summary_t.groupby('rank'):
        pd_rank = pd_group.set_index('Tool').T
        rank_to_html[rank] = [create_table_html(pd_rank)]

    if not rank_to_html:
        return None

    figures = [create_tax_figure(tool, df_tool) for tool, df_tool in df_summary_t.groupby('Tool')]

    taxonomic_div = Div(text="""<div>{}</div>""".format(rank_to_html[load_ncbi_taxinfo.RANKS[0]][0]), css_classes=['bk-width-auto'])
    source = ColumnDataSource(data=rank_to_html)

    select_rank = Select(title="Taxonomic rank:", value=load_ncbi_taxinfo.RANKS[0], options=load_ncbi_taxinfo.RANKS, css_classes=['bk-fit-content'])
    select_rank_callback = CustomJS(args=dict(source=source), code="""
        mytable.text = source.data[select_rank.value][0];
    """)
    select_rank.js_on_change('value', select_rank_callback)
    select_rank_callback.args["mytable"] = taxonomic_div
    select_rank_callback.args["select_rank"] = select_rank
    metrics_row_t = row(column(select_rank, create_heatmap_div(), taxonomic_div, sizing_mode='scale_width', css_classes=['bk-width-auto']), column(row(figures), sizing_mode='scale_width', css_classes=['bk-width-auto']),css_classes=['bk-width-auto'], sizing_mode='scale_width')
    metrics_panel = Panel(child=metrics_row_t, title="Metrics")

    bins_columns = OrderedDict([('id', 'Taxon ID'), ('rank', 'Rank'), ('purity', 'Purity'), ('completeness', 'Completeness'), ('predicted_size', 'Predicted size'), ('true_positives', 'True positives'), ('real_size', 'Real size')])
    metrics_bins_panel = create_metrics_per_bin_panel(pd_bins[pd_bins['rank'] != 'NA'], bins_columns)

    tabs = Tabs(tabs=[metrics_panel, metrics_bins_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])

    return tabs


def create_html(df_summary, pd_bins, output_dir, desc_text):
    create_heatmap_bar(output_dir)
    tabs_list = []

    metrics_row_g = create_genome_binning_html(df_summary, pd_bins)
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
