# Copyright 2023 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import io
import os
from collections import OrderedDict
from collections import defaultdict
import datetime
import numpy as np
import pandas as pd
import math
import re
import logging
from jinja2 import Template
from statistics import median
from src import plots
from src import binning_classes
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
from bokeh.models.widgets import TableColumn, Div, Select
from bokeh.models import TabPanel, Tabs
from bokeh.models import (DataTable,
                          CustomJS,
                          Legend, LegendItem,
                          Band,
                          CustomJSTickFormatter,
                          HoverTool)
from bokeh.models.tickers import FixedTicker
from bokeh.resources import INLINE
from bokeh.embed import file_html
from bokeh.plotting import ColumnDataSource


AVG_OVER_SAMPLES = '[average over samples]'
CLICK_ON_LEGENDS_DIV = '<div>' + AVG_OVER_SAMPLES + '</div><div style="padding-top: 10px;">Click on the legends to enable/disable the results for a tool.</div>'

TEMPLATE = Template('''
{% from macros import embed %}
<!DOCTYPE html>
<html lang="en">
  {% block head %}
  <head>
  {% block inner_head %}
    <meta charset="utf-8">
    <title>{% block title %}{{ title | e if title else "Bokeh Plot" }}{% endblock %}</title>
  {%  block preamble -%}{%- endblock %}
  {%  block resources %}
    <style>
      html, body {
        box-sizing: border-box;
        height: 100%;
        margin: 8px;
        padding: 0;
      }
      .bk-fit-content {width: fit-content; width: -moz-fit-content;}
      .bk-width-auto {width: auto !important; height: auto !important;}
      .bk-display-block {display: block !important;}
      .bk-float-left {float: left;}
      bk-width-auto-main>div {width: -webkit-fill-available !important;}
      div.bk-width-auto-main {width: -webkit-fill-available !important;}
      .bk-tabs-margin{margin-top: 20px !important;}
      .bk-tabs-margin-lr{margin-left: 10px; margin-right: 10px}
      .bk-root {display: flex; justify-content: center;}
      .bk-padding-top {padding-top: 10px;}
      .bk-padding-top2 {padding-top: 20px;}
      .bk-combo-box > div:first-child {
        width: auto !important;
        padding-right: 40px;
      }
    </style>
  {%   block css_resources -%}
    {{- bokeh_css if bokeh_css }}
  {%-  endblock css_resources %}
  {%   block js_resources -%}
    {{  bokeh_js if bokeh_js }}
  {%-  endblock js_resources %}
  {%  endblock resources %}
  {%  block postamble %}{% endblock %}
  {% endblock inner_head %}
  </head>
  {% endblock head%}
  {% block body %}
  <body>
  {%  block inner_body %}
  {%    block contents %}
  {%      for doc in docs %}
  {{        embed(doc) if doc.elementid }}
  {%-       for root in doc.roots %}
  {%          block root scoped %}
  {{            embed(root) }}
  {%          endblock %}
  {%        endfor %}
  {%      endfor %}
  {%    endblock contents %}
  {{ plot_script | indent(4) }}
  {%  endblock inner_body %}
  </body>
  {% endblock body%}
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

    fig.savefig(os.path.join(output_dir, 'heatmap_bar.png'), dpi=100, format='png', bbox_inches='tight', pad_inches=-.001, transparent=True)
    pltx.close(fig)


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def get_colors_and_ranges(name, all_values):
    color1 = 'dodgerblue'
    color2 = 'red'
    hue1 = 12
    hue2 = 240

    metrics = [utils_labels.MISCLASSIFICATION_PER_BP, utils_labels.MISCLASSIFICATION_PER_SEQ]
    if name in metrics + [x + utils_labels.UNFILTERED for x in metrics]:
        return color2, color1, hue2, hue1, 0, 1
    metrics = [utils_labels.UNIFRAC_BP, utils_labels.UNIFRAC_SEQ]
    if name in metrics + [x + utils_labels.UNFILTERED for x in metrics]:
        return color2, color1, hue2, hue1, 0, max(all_values)
    return color1, color2, hue1, hue2, 0, 1


def get_heatmap_colors(pd_series, **args):
    values = pd_series.tolist()

    metrics = [utils_labels.AVG_PRECISION_BP_SEM, utils_labels.AVG_RECALL_BP_SEM, utils_labels.AVG_PRECISION_SEQ_SEM,
               utils_labels.AVG_RECALL_SEQ_SEM, utils_labels.AVG_RECALL_BP_SEM_CAMI1, utils_labels.AVG_RECALL_SEQ_SEM_CAMI1]
    if pd_series.name in metrics + [x + utils_labels.UNFILTERED for x in metrics]:
        return ['background-color: white' for x in values]

    dropped_gs = False
    if pd_series.index[0] == utils_labels.GS:
        values = values[1:]
        dropped_gs = True
    if len(values) == 0:
        return ['']

    notnan_values = [x for x in values if isinstance(x, (float, int)) and not np.isnan(x)]
    if not notnan_values:
        red = 'background-color: red'
        return [red for x in values] if not dropped_gs else [''] + [red for x in values]

    color1, color2, hue1, hue2, min_value, max_value = get_colors_and_ranges(pd_series.name, values)

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

    if dropped_gs:
        return [''] + return_colors
    else:
        return return_colors


def create_title_div(id, name, info):
    div = Div(text="""<div style="text-align:left;font-size: 20pt;font-weight: bold;">{1}<span style="float: right;font-size: 10pt;font-weight:normal;">{2}</span>""".format(id, name, info), styles={'width': 'auto'})
    div.stylesheets = ["""
    .bk-clearfix {
      width: -webkit-fill-available;
      width: -moz-available;
    }
    """]
    return div


def create_metrics_per_bin_panel(pd_bins, bins_columns, sample_ids_list, output_dir, binning_type):
    styles = [{'selector': 'td', 'props': [('width', '92pt')]},
              {'selector': 'th', 'props': [('width', '92pt'), ('text-align', 'left')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]

    tools = pd_bins[utils_labels.TOOL].unique().tolist()
    tool_to_sample_to_html = defaultdict(list)

    for tool, pd_tool_groupby in pd_bins.groupby(utils_labels.TOOL):
        pd_tool_sample_groupby = pd_tool_groupby.groupby('sample_id')
        available_samples = pd_tool_sample_groupby.groups.keys()
        for sample_id in sample_ids_list:
            if sample_id not in available_samples:
                tool_to_sample_to_html[tool].append('')
                continue
            pd_tool_sample = pd_tool_sample_groupby.get_group(sample_id)
            pd_tool_sample = pd_tool_sample[list(bins_columns.keys())].rename(columns=dict(bins_columns))
            if binning_type == 'taxonomic':
                pd_tool_sample['Taxon ID'] = pd_tool_sample['Taxon ID'].astype('Int64')
            tool_sample_html = pd_tool_sample.head(500).style.set_table_styles(styles).format(precision=3).hide(axis='index').to_html()
            tool_sample_html += '<div style="padding-top: 20px; padding-bottom: 20px;">{}</div>'.format('Complete table available in: ' + os.path.join(output_dir, binning_type, tool, 'metrics_per_bin.tsv'))
            tool_to_sample_to_html[tool].append(tool_sample_html)
    tool_to_sample_to_html['all_samples'] = sample_ids_list

    table_div = Div(text="""<div>{}</div>""".format(tool_to_sample_to_html[tools[0]][0]))

    source = ColumnDataSource(data=tool_to_sample_to_html)

    select_tool = Select(title="Binner:", value=tools[0], options=tools)
    select_sample = Select(title="Sample:", value='0', options=list(zip(map(str, range(len(sample_ids_list))), sample_ids_list)))
    select_tool_sample_callback = CustomJS(args=dict(source=source), code="""
        select_sample.options = [];
        const options_array = [];
        for (let index in source.data[select_tool.value]) {
            if (source.data[select_tool.value][index] != "") {
                options_array.push([index, source.data["all_samples"][index]])
            }
        }
        select_sample.options = options_array;
        if (source.data[select_tool.value][select_sample.value] != "") {
            mytable.text = source.data[select_tool.value][select_sample.value];
        } else {
            mytable.text = source.data[select_tool.value][options_array[0][0]];
        }
    """)
    select_tool.js_on_change('value', select_tool_sample_callback)
    select_sample.js_on_change('value', select_tool_sample_callback)
    select_tool_sample_callback.args["mytable"] = table_div
    select_tool_sample_callback.args["select_tool"] = select_tool
    select_tool_sample_callback.args["select_sample"] = select_sample

    table_column = row(column(row(select_tool, select_sample), table_div, sizing_mode='scale_width'), sizing_mode='scale_width')
    metrics_bins_panel = TabPanel(child=table_column, title="Metrics per bin")
    return metrics_bins_panel


def create_contamination_completeness_table(pd_bins, min_completeness, max_contamination):
    df = binning_classes.GenomeQuery.calc_num_recovered_genomes(pd_bins, min_completeness, max_contamination)

    def create_table_column(field):
        return TableColumn(title=field, field=field, width=600)

    dt = DataTable(source=ColumnDataSource(df),
                   columns=list(map(lambda x: create_table_column(x), df.columns.values)),
                   width=800,
                   height=1000,
                   reorderable=True,
                   selectable=True)
    return [column(dt)]


def create_heatmap_div():
    heatmap_legend = '<img src="heatmap_bar.png" /><div style="text-align:left;font-size: 11px;">Worst<span style="float:right;">Best</span><span style="margin-right: 36px;float:right;">Median</span></div>'
    heatmap_legend_div = Div(text=heatmap_legend, styles={"width": "155px", "margin-bottom": "-10px"})
    return heatmap_legend_div


def get_labels_genome():
    return ([(utils_labels.ACCURACY_PER_BP, utils_labels.TOOLTIP_ACCURACY_PER_BP),
             (utils_labels.ACCURACY_PER_SEQ, utils_labels.TOOLTIP_ACCURACY_PER_SEQ),
             (utils_labels.MISCLASSIFICATION_PER_BP, utils_labels.TOOLTIP_MISCLASSIFICATION_PER_BP),
             (utils_labels.MISCLASSIFICATION_PER_SEQ, utils_labels.TOOLTIP_MISCLASSIFICATION_PER_SEQ),
             (utils_labels.PRECISION_PER_BP, utils_labels.TOOLTIP_PRECISION_PER_BP),
             (utils_labels.PRECISION_PER_SEQ, utils_labels.TOOLTIP_PRECISION_PER_SEQ),
             (utils_labels.RECALL_PER_BP, utils_labels.TOOLTIP_RECALL_PER_BP),
             (utils_labels.RECALL_PER_SEQ, utils_labels.TOOLTIP_RECALL_PER_SEQ),
             (utils_labels.F1_SCORE_PER_BP, utils_labels.TOOLTIP_F1_SCORE_PER_BP),
             (utils_labels.F1_SCORE_PER_SEQ, utils_labels.TOOLTIP_F1_SCORE_PER_SEQ),
             (utils_labels.RI_BY_BP, utils_labels.TOOLTIP_RI_BY_BP),
             (utils_labels.RI_BY_SEQ, utils_labels.TOOLTIP_RI_BY_SEQ),
             (utils_labels.ARI_BY_BP, utils_labels.TOOLTIP_ARI_BY_BP),
             (utils_labels.ARI_BY_SEQ, utils_labels.TOOLTIP_ARI_BY_SEQ),
             (utils_labels.AVG_PRECISION_BP, utils_labels.TOOLTIP_AVG_PRECISION_BP),
             (utils_labels.AVG_PRECISION_BP_SEM, utils_labels.TOOLTIP_AVG_PRECISION_BP_SEM),
             (utils_labels.AVG_PRECISION_SEQ, utils_labels.TOOLTIP_AVG_PRECISION_SEQ),
             (utils_labels.AVG_PRECISION_SEQ_SEM, utils_labels.TOOLTIP_AVG_PRECISION_SEQ_SEM),
             (utils_labels.AVG_RECALL_BP, utils_labels.TOOLTIP_AVG_RECALL_BP),
             (utils_labels.AVG_RECALL_BP_CAMI1, utils_labels.TOOLTIP_AVG_RECALL_BP_CAMI1),
             (utils_labels.AVG_RECALL_BP_SEM, utils_labels.TOOLTIP_AVG_RECALL_BP_SEM),
             (utils_labels.AVG_RECALL_BP_SEM_CAMI1, utils_labels.TOOLTIP_AVG_RECALL_BP_SEM),
             (utils_labels.AVG_RECALL_SEQ, utils_labels.TOOLTIP_AVG_RECALL_SEQ),
             (utils_labels.AVG_RECALL_SEQ_CAMI1, utils_labels.TOOLTIP_AVG_RECALL_SEQ_CAMI1),
             (utils_labels.F1_SCORE_BP, utils_labels.TOOLTIP_F1_SCORE_BP),
             (utils_labels.F1_SCORE_SEQ, utils_labels.TOOLTIP_F1_SCORE_SEQ),
             (utils_labels.F1_SCORE_BP_CAMI1, utils_labels.TOOLTIP_F1_SCORE_BP),
             (utils_labels.F1_SCORE_SEQ_CAMI1, utils_labels.TOOLTIP_F1_SCORE_SEQ),
             (utils_labels.AVG_RECALL_SEQ_SEM, utils_labels.TOOLTIP_AVG_RECALL_SEQ_SEM),
             (utils_labels.AVG_RECALL_SEQ_SEM_CAMI1, utils_labels.TOOLTIP_AVG_RECALL_SEQ_SEM),
             (utils_labels.PERCENTAGE_ASSIGNED_BPS, utils_labels.TOOLTIP_PERCENTAGE_ASSIGNED_BPS),
             (utils_labels.PERCENTAGE_ASSIGNED_SEQS, utils_labels.TOOLTIP_PERCENTAGE_ASSIGNED_SEQS),
             (utils_labels.UNIFRAC_BP, utils_labels.TOOLTIP_UNIFRAC_BP),
             (utils_labels.UNIFRAC_SEQ, utils_labels.TOOLTIP_UNIFRAC_SEQ)])


def get_labels_taxonomic():
    return ([(utils_labels.ACCURACY_PER_BP, utils_labels.TOOLTIP_ACCURACY_PER_BP_TAX),
             (utils_labels.ACCURACY_PER_SEQ, utils_labels.TOOLTIP_ACCURACY_PER_SEQ_TAX),
             (utils_labels.MISCLASSIFICATION_PER_BP, utils_labels.TOOLTIP_MISCLASSIFICATION_PER_BP),
             (utils_labels.MISCLASSIFICATION_PER_SEQ, utils_labels.TOOLTIP_MISCLASSIFICATION_PER_SEQ),
             (utils_labels.PRECISION_PER_BP, utils_labels.TOOLTIP_PRECISION_PER_BP_TAX),
             (utils_labels.PRECISION_PER_SEQ, utils_labels.TOOLTIP_PRECISION_PER_SEQ_TAX),
             (utils_labels.RECALL_PER_BP, utils_labels.TOOLTIP_RECALL_PER_BP_TAX),
             (utils_labels.RECALL_PER_SEQ, utils_labels.TOOLTIP_RECALL_PER_SEQ_TAX),
             (utils_labels.F1_SCORE_PER_BP, utils_labels.TOOLTIP_F1_SCORE_PER_BP),
             (utils_labels.F1_SCORE_PER_SEQ, utils_labels.TOOLTIP_F1_SCORE_PER_SEQ),
             (utils_labels.RI_BY_BP, utils_labels.TOOLTIP_RI_BY_BP_TAX),
             (utils_labels.RI_BY_SEQ, utils_labels.TOOLTIP_RI_BY_SEQ_TAX),
             (utils_labels.ARI_BY_BP, utils_labels.TOOLTIP_ARI_BY_BP),
             (utils_labels.ARI_BY_SEQ, utils_labels.TOOLTIP_ARI_BY_SEQ),
             (utils_labels.AVG_PRECISION_BP, utils_labels.TOOLTIP_AVG_PRECISION_BP_TAX),
             (utils_labels.AVG_PRECISION_BP_SEM, utils_labels.TOOLTIP_AVG_PRECISION_BP_SEM),
             (utils_labels.AVG_PRECISION_SEQ, utils_labels.TOOLTIP_AVG_PRECISION_SEQ_TAX),
             (utils_labels.AVG_PRECISION_SEQ_SEM, utils_labels.TOOLTIP_AVG_PRECISION_SEQ_SEM),
             (utils_labels.AVG_RECALL_BP, utils_labels.TOOLTIP_AVG_RECALL_BP_TAX),
             (utils_labels.AVG_RECALL_BP_CAMI1, utils_labels.TOOLTIP_AVG_RECALL_BP_CAMI1),
             (utils_labels.AVG_RECALL_BP_SEM, utils_labels.TOOLTIP_AVG_RECALL_BP_SEM),
             (utils_labels.AVG_RECALL_SEQ, utils_labels.TOOLTIP_AVG_RECALL_SEQ_TAX),
             (utils_labels.AVG_RECALL_SEQ_CAMI1, utils_labels.TOOLTIP_AVG_RECALL_SEQ_CAMI1),
             (utils_labels.F1_SCORE_BP, utils_labels.TOOLTIP_F1_SCORE_BP),
             (utils_labels.F1_SCORE_SEQ, utils_labels.TOOLTIP_F1_SCORE_SEQ),
             (utils_labels.AVG_RECALL_SEQ_SEM, utils_labels.TOOLTIP_AVG_RECALL_SEQ_SEM),
             (utils_labels.PERCENTAGE_ASSIGNED_BPS, utils_labels.TOOLTIP_PERCENTAGE_ASSIGNED_BPS),
             (utils_labels.PERCENTAGE_ASSIGNED_SEQS, utils_labels.TOOLTIP_PERCENTAGE_ASSIGNED_SEQS),
             (utils_labels.UNIFRAC_BP, utils_labels.TOOLTIP_UNIFRAC_BP),
             (utils_labels.UNIFRAC_SEQ, utils_labels.TOOLTIP_UNIFRAC_SEQ)])


def create_table_html(df_summary, is_taxonomic=False, include_cami1=False):
    metrics1 = [utils_labels.AVG_PRECISION_BP,
                utils_labels.AVG_PRECISION_SEQ,
                utils_labels.AVG_RECALL_BP,
                utils_labels.AVG_RECALL_SEQ,
                utils_labels.F1_SCORE_BP,
                utils_labels.F1_SCORE_SEQ,
                utils_labels.AVG_PRECISION_BP_SEM,
                utils_labels.AVG_PRECISION_SEQ_SEM,
                utils_labels.AVG_RECALL_BP_SEM,
                utils_labels.AVG_RECALL_SEQ_SEM]

    if is_taxonomic and utils_labels.AVG_PRECISION_BP + utils_labels.UNFILTERED in df_summary.index:
        metrics1 += [x + utils_labels.UNFILTERED for x in metrics1]

    if include_cami1:
        metrics1 += [utils_labels.AVG_RECALL_BP_CAMI1,
                     utils_labels.AVG_RECALL_SEQ_CAMI1,
                     utils_labels.F1_SCORE_BP_CAMI1,
                     utils_labels.F1_SCORE_SEQ_CAMI1,
                     utils_labels.AVG_RECALL_BP_SEM_CAMI1,
                     utils_labels.AVG_RECALL_SEQ_SEM_CAMI1]
    metrics2 = [utils_labels.ACCURACY_PER_BP,
                utils_labels.ACCURACY_PER_SEQ,
                utils_labels.MISCLASSIFICATION_PER_BP,
                utils_labels.MISCLASSIFICATION_PER_SEQ,
                utils_labels.PRECISION_PER_BP,
                utils_labels.PRECISION_PER_SEQ,
                utils_labels.RECALL_PER_BP,
                utils_labels.RECALL_PER_SEQ,
                utils_labels.F1_SCORE_PER_BP,
                utils_labels.F1_SCORE_PER_SEQ,
                utils_labels.RI_BY_BP,
                utils_labels.RI_BY_SEQ,
                utils_labels.ARI_BY_BP,
                utils_labels.ARI_BY_SEQ,
                utils_labels.PERCENTAGE_ASSIGNED_BPS,
                utils_labels.PERCENTAGE_ASSIGNED_SEQS]
    if is_taxonomic:
        metrics2.append(utils_labels.UNIFRAC_BP)
        metrics2.append(utils_labels.UNIFRAC_SEQ)

    if is_taxonomic and utils_labels.AVG_PRECISION_BP + utils_labels.UNFILTERED in df_summary.index:
        metrics2 += [x + utils_labels.UNFILTERED for x in metrics2]

    all_metrics = [metrics1, metrics2]
    metrics1_label = utils_labels.QUALITY_OF_BINS
    metrics2_label = utils_labels.QUALITY_OF_SAMPLE
    all_metrics_labels = [metrics1_label, metrics2_label]

    styles = [{'selector': 'td', 'props': [('width', '115pt')]},
              {'selector': 'th', 'props': [('width', '115pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '210pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    def get_html_dict(metrics):
        d_dict = {}
        for tuple in metrics:
            d_dict[tuple[0] + '<'] = '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div><'.format(utils_labels.LABELS[tuple[0]], utils_labels.LABELS[tuple[0]], tuple[1])
            d_dict[tuple[0] + utils_labels.UNFILTERED + '<'] = '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div><'.format(utils_labels.LABELS[tuple[0]] + utils_labels.UNFILTERED, utils_labels.LABELS[tuple[0]], tuple[1])
        return d_dict

    d = get_html_dict(get_labels_taxonomic() if is_taxonomic else get_labels_genome())

    pattern = re.compile('|'.join(map(re.escape, d)))

    def translate(match):
        return d[match.group(0)]

    df_summary.index.name = None

    html = io.StringIO()
    first_metrics = True
    for metrics, metrics_label in zip(all_metrics, all_metrics_labels):
        html.write('<p style="margin-bottom: auto"><b>{}</b></p>'.format(metrics_label))
        df_metrics = df_summary.loc[metrics]
        sorted_columns = df_metrics.columns.tolist()
        df_metrics = df_metrics.loc[:, sorted_columns]

        if first_metrics:
            this_style = styles
            first_metrics = False
        else:
            this_style = styles_hidden_thead
        html.write(df_metrics.style.apply(get_heatmap_colors, df_metrics=df_metrics, axis=1).format(
            precision=3).set_table_styles(this_style).to_html())
    html = html.getvalue()
    html = pattern.sub(translate, html)

    tooltips = """
    <style>
        .tooltip {
            position: relative;
            display: inline-block;
            border-bottom: 1px dashed lightgray;
            cursor: help;
        }
        div.tooltip span.tooltiptext {
            visibility: hidden;
            width: 280px;
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
        div.tooltip span.tooltiptext::after {
            content: "";
            position: absolute;
            top: 100%;
            left: 50%;
            margin-left: -5px;
            border-width: 5px;
            border-style: solid;
            border-color: #555 transparent transparent transparent;
        }
        div.tooltip:hover span.tooltiptext {
            visibility: visible;
            opacity: 1;
        }
        tr:hover {
            outline: 1px solid black;
        }
    </style>"""

    return '{}<div style="margin-bottom:10pt;">{}</div>'.format(tooltips, html)


def create_precision_recall_figure(df_summary, xname1, yname1, xname2, yname2, title):
    colors_list = plots.create_colors_list()
    bokeh_colors = [matplotlib.colors.to_hex(c) for c in colors_list]

    legend_it = []
    tooltips1 = [(utils_labels.TOOL, '@index'),
                 (xname1, '@{' + xname1 + '}'),
                 (yname1, '@{' + yname1 + '}')]
    tooltips2 = [(utils_labels.TOOL, '@index'),
                 (xname2, '@{' + xname2 + '}'),
                 (yname2, '@{' + yname2 + '}')]
    p = figure(title=title, width=580, height=400, x_range=(0, 1), y_range=(0, 1), toolbar_location="below")
    p.xaxis.axis_label = utils_labels.LABELS[xname1]
    p.yaxis.axis_label = utils_labels.LABELS[yname1]
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'
    for color, (index, row) in zip(bokeh_colors, df_summary.iterrows()):
        tool = index
        source = ColumnDataSource(data=row.to_frame().T)
        pcircle = p.circle(xname1, yname1, source=source, color=color, fill_alpha=0.2, size=10)
        px = p.x(xname2, yname2, source=source, color=color, size=10)
        p.add_tools(HoverTool(tooltips=tooltips1, renderers=[pcircle], toggleable=False))
        p.add_tools(HoverTool(tooltips=tooltips2, renderers=[px], toggleable=False))
        legend_it.append(LegendItem(label=tool, renderers=[pcircle, px]))
    pcircle = p.circle([0], [0], color='black', fill_alpha=0, size=0)
    legend_it.append(LegendItem(label='by bp', renderers=[pcircle]))
    px = p.x([0], [0], color='black', fill_alpha=0, size=0)
    legend_it.append(LegendItem(label='by seq', renderers=[px]))

    p.add_layout(Legend(items=legend_it), 'right')
    p.legend.click_policy = 'hide'
    return p


def create_precision_recall_all_genomes_scatter(pd_genome_bins, tools):
    colors_list = plots.create_colors_list()
    bokeh_colors = [matplotlib.colors.to_hex(c) for c in colors_list]

    p = figure(title='Quality per bin', width=580, height=400, x_range=(0, 1), y_range=(0, 1), toolbar_location="below")
    p.add_tools(HoverTool(tooltips=[('Sample', '@sample_id'),
                                    ('Genome', '@mapping_id'),
                                    ('Purity of bin (bp)', '@precision_bp'),
                                    ('Completeness of bin (bp)', '@recall_bp')], toggleable=False))

    legend_it = []
    for color, tool in zip(bokeh_colors, tools):
        source = ColumnDataSource(data=pd_genome_bins[pd_genome_bins[utils_labels.TOOL] == tool])
        pcircle = p.circle('precision_bp', 'recall_bp', color=color, alpha=0.8, source=source)
        legend_it.append((tool, [pcircle]))
    p.add_layout(Legend(items=legend_it), 'right')
    p.xaxis.axis_label = 'Purity per bin (bp)'
    p.yaxis.axis_label = 'Completeness per bin (bp)'
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'
    p.legend.click_policy = 'hide'
    return p


def create_contamination_plot(pd_bins, tools, title, xlabel, ylabel, create_column_function):
    pd_bins_copy = pd_bins[[utils_labels.TOOL, 'precision_bp', 'recall_bp']].copy().dropna(subset=['precision_bp'])
    create_column_function(pd_bins_copy)

    colors_list = plots.create_colors_list()
    bokeh_colors = [matplotlib.colors.to_hex(c) for c in colors_list]

    p = figure(title=title, width=580, height=400, toolbar_location="below")
    p.x_range.start = 0
    # p.y_range.start = pd_bins_copy['newcolumn'].min()
    p.y_range.end = 1
    legend_it = []
    for color, tool in zip(bokeh_colors, tools):
        pd_tool_bins = pd_bins_copy[pd_bins_copy[utils_labels.TOOL] == tool]
        pd_tool_bins = pd_tool_bins.sort_values(by='newcolumn', ascending=False).reset_index()
        pd_tool_bins = pd_tool_bins.drop(['index'], axis=1)

        source = ColumnDataSource(data=pd_tool_bins)

        pline = p.line('index', 'newcolumn', color=color, line_width=2, source=source)
        legend_it.append((tool, [pline]))
    p.add_layout(Legend(items=legend_it), 'right')
    p.xaxis.axis_label = xlabel
    p.yaxis.axis_label = ylabel
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'
    p.legend.click_policy = 'hide'
    return p


def create_tax_figure(tool, df_summary, metrics_list, errors_list):
    rank_indices = list(range(len(load_ncbi_taxinfo.RANKS)))

    df = df_summary.set_index("rank").reindex(load_ncbi_taxinfo.RANKS)
    df["x"] = rank_indices
    for metric, error in zip(metrics_list, errors_list):
        if error:
            df[metric + "upper"] = df[metric] + df[error]
            df[metric + "lower"] = df[metric] - df[error]
    source = ColumnDataSource(df.reset_index())

    p = figure(width=500, height=550, x_range=(0, 7), y_range=(0, 1))
    line_colors = ["#006cba", "#008000", "#ba9e00", "red"]
    legend_it = []
    for i, (metric, error, color) in enumerate(zip(metrics_list, errors_list, line_colors)):
        pline = p.line(x='x', y=metric, line_color=color, source=source, line_width=2)
        legend_it.append(LegendItem(label=utils_labels.LABELS[metric], renderers=[pline]))
        if error:
            band = Band(base='x', lower=metric + "lower", upper=metric + "upper", source=source, level='underlay',
                              fill_alpha=.3, line_width=1, line_color='black', fill_color=color)
            p.add_layout(band)
            callback = CustomJS(args = dict(band=band), code="""
                if (band.visible == false)
                    band.visible = true;
                else
                    band.visible = false; """)
            pline.js_on_change('visible', callback)
            if i > 1:
                pline.visible = False
                band.visible = False

    p.xaxis.ticker = FixedTicker(ticks=rank_indices)
    p.xaxis.formatter = CustomJSTickFormatter(code="""
        var mapping = {0: "superkingdom", 1: "phylum", 2: "class", 3: "order", 4: "family", 5: "genus", 6: "species", 7: "strain"};
        return mapping[tick];
    """)
    p.xaxis.major_label_orientation = 1.2
    p.add_layout(Legend(items=legend_it, location=(10, 10)), 'above')
    p.xgrid[0].grid_line_color=None
    p.ygrid[0].grid_line_alpha=0.5
    p.title.text = tool
    p.legend.click_policy = 'hide'
    return p


def create_rankings_table(df_summary, show_rank=False):
    columns = [utils_labels.AVG_PRECISION_BP,
               utils_labels.AVG_RECALL_BP,
               utils_labels.PRECISION_PER_BP,
               utils_labels.RECALL_PER_BP,
               utils_labels.ARI_BY_SEQ,
               utils_labels.ARI_BY_BP,
               utils_labels.PERCENTAGE_ASSIGNED_BPS,
               utils_labels.ACCURACY_PER_BP]
    if show_rank:
        columns.insert(0, utils_labels.RANK)
    labels_dict = utils_labels.LABELS.copy()
    labels_dict[utils_labels.RANK] = 'Taxonomic rank'
    pd_rankings = df_summary[columns].rename(columns=labels_dict).round(decimals=5).reset_index()

    def create_table_column(field):
        return TableColumn(title=field, field=field, width=100)

    dt = DataTable(source=ColumnDataSource(pd_rankings),
                   columns=list(map(lambda x: create_table_column(x), pd_rankings.columns.values)),
                   width=1500,
                   height=1000,
                   reorderable=True,
                   selectable=True)
    return [column(dt)]


def create_genome_binning_plots_panel(pd_bins, pd_mean):
    click_div = Div(text=CLICK_ON_LEGENDS_DIV, styles={"width": "500px", "margin-top": "15px; margin-bottom:5px;"})
    purity_completeness_plot = column(create_precision_recall_figure(pd_mean, utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP, utils_labels.AVG_PRECISION_SEQ, utils_labels.AVG_RECALL_SEQ, utils_labels.QUALITY_OF_BINS))
    purity_recall_bp_plot = column(create_precision_recall_figure(pd_mean, utils_labels.PRECISION_PER_BP, utils_labels.RECALL_PER_BP, utils_labels.PRECISION_PER_SEQ, utils_labels.RECALL_PER_SEQ, utils_labels.QUALITY_OF_SAMPLE))

    all_samples_div = Div(text='<div style="padding-top: 20px;">All samples</div>', styles={"width": "500px", "margin-top": "15px; margin-bottom:5px;"})
    all_bins_plot = column(create_precision_recall_all_genomes_scatter(pd_bins, pd_mean.index.tolist()))
    completeness_contamination_plot = column(create_contamination_plot(pd_bins, pd_mean.index.tolist(), 'Completeness - contamination', 'Index of bin (sorted by completeness - contamination (bp))', 'Completeness - contamination (bp)', plots.create_completeness_minus_contamination_column))
    contamination_plot = column(create_contamination_plot(pd_bins, pd_mean.index.tolist(), 'Contamination', 'Index of bin (sorted by contamination (bp))', 'Contamination (bp)', plots.create_contamination_column))

    return TabPanel(child=column([click_div, purity_completeness_plot, purity_recall_bp_plot, all_samples_div, all_bins_plot, completeness_contamination_plot, contamination_plot], sizing_mode='scale_width'), title='Plots')


def create_genome_binning_html(df_summary, pd_bins, labels, sample_ids_list, options):
    if pd_bins.empty:
        return None
    df_summary_g = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'genome']
    if df_summary_g.empty:
        return None

    sample_to_html = {}
    available_samples = list(df_summary_g[utils_labels.SAMPLE].unique())
    available_samples = [sample_id for sample_id in sample_ids_list if sample_id in available_samples]
    available_tools = list(df_summary_g[utils_labels.TOOL].unique())
    available_tools = [tool for tool in labels if tool in available_tools]

    pd_mean = df_summary_g.groupby(utils_labels.TOOL).mean(numeric_only=True).reindex(available_tools)
    pd_mean[utils_labels.SAMPLE] = AVG_OVER_SAMPLES
    sample_to_html[AVG_OVER_SAMPLES] = [create_table_html(pd_mean.T, include_cami1=True)]

    for sample_id, pd_sample in df_summary_g.groupby(utils_labels.SAMPLE):
        sample_to_html[sample_id] = [create_table_html(pd_sample.set_index(utils_labels.TOOL).T, include_cami1=True)]

    sample_ids_list_combo = [AVG_OVER_SAMPLES] + available_samples

    genome_div = Div(text="""<div style="margin-bottom:10pt;">{}</div>""".format(sample_to_html[sample_ids_list_combo[0]][0]))
    source = ColumnDataSource(data=sample_to_html)

    select_sample = Select(title="Sample:", value=sample_ids_list_combo[0], options=sample_ids_list_combo)
    select_sample_callback = CustomJS(args=dict(source=source), code="""
        mytable.text = source.data[select_sample.value][0];
    """)
    select_sample.js_on_change('value', select_sample_callback)
    select_sample_callback.args["mytable"] = genome_div
    select_sample_callback.args["select_sample"] = select_sample

    metrics_column = column(column(select_sample, create_heatmap_div(), genome_div, sizing_mode='scale_width'), sizing_mode='scale_width')
    metrics_panel = TabPanel(child=metrics_column, title="Metrics")

    plots_panel = create_genome_binning_plots_panel(pd_bins, pd_mean)

    bins_columns = utils_labels.get_genome_bins_columns()
    metrics_bins_panel = create_metrics_per_bin_panel(pd_bins, bins_columns, sample_ids_list, options.output_dir, 'genome')

    cc_table = create_contamination_completeness_table(pd_bins, options.min_completeness, options.max_contamination)
    cc_panel = TabPanel(child=row(cc_table), title="#Recovered genomes")

    rankings_panel = TabPanel(child=column([Div(text="Click on the columns header for sorting.", styles={"width": "500px", "margin-top": "20px"}),
                                        row(create_rankings_table(pd_mean.reset_index().set_index([utils_labels.SAMPLE, utils_labels.TOOL])))]), title="Rankings")

    tabs = Tabs(tabs=[metrics_panel, plots_panel, metrics_bins_panel, rankings_panel, cc_panel])

    return tabs


def create_plots_per_binner(df_summary_t):
    tools = df_summary_t[utils_labels.TOOL].unique().tolist()
    sample_div = Div(text=AVG_OVER_SAMPLES, styles={"margin-top": "15px; margin-bottom:5px;"})
    # using the same div breaks the layout
    sample_div2 = Div(text=AVG_OVER_SAMPLES, styles={"margin-top": "15px; margin-bottom:5px;"})

    metrics_list = [utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP, utils_labels.AVG_PRECISION_SEQ, utils_labels.AVG_RECALL_SEQ]
    errors_list = [utils_labels.AVG_PRECISION_BP_SEM, utils_labels.AVG_RECALL_BP_SEM, utils_labels.AVG_PRECISION_SEQ_SEM, utils_labels.AVG_RECALL_SEQ_SEM]
    tools_figures = [create_tax_figure(tool, df_summary_t[df_summary_t[utils_labels.TOOL] == tool], metrics_list, errors_list) for tool in tools]
    tools_figures_columns = [column(x) for x in tools_figures]
    tools_column_unweighted = column([sample_div] + tools_figures_columns, sizing_mode='scale_width')

    metrics_list = [utils_labels.PRECISION_PER_BP, utils_labels.RECALL_PER_BP, utils_labels.PRECISION_PER_SEQ, utils_labels.RECALL_PER_SEQ]
    errors_list = ['', '', '', '']
    tools_figures_weighted = [create_tax_figure(tool, df_summary_t[df_summary_t[utils_labels.TOOL] == tool], metrics_list, errors_list) for tool in tools]
    tools_figures_columns = [column(x) for x in tools_figures_weighted]
    tools_column_weighted = column([sample_div2] + tools_figures_columns, sizing_mode='scale_width')

    tools_unweighted_panel = TabPanel(child=tools_column_unweighted, title=utils_labels.QUALITY_OF_BINS)
    tools_weighted_panel = TabPanel(child=tools_column_weighted, title=utils_labels.QUALITY_OF_SAMPLE)
    tools_tabs = Tabs(tabs=[tools_unweighted_panel, tools_weighted_panel])
    return tools_tabs


def create_tax_ranks_panel(qbins_plots_list, qsamples_plots_list, cc_plots_dict_list, contamination_plots_list):
    sample_div = Div(text=CLICK_ON_LEGENDS_DIV, styles={"width": "500px", "margin-top": "20px"})
    all_samples_div = Div(text='<div style="padding-top: 20px;">All samples</div>', styles={"width": "500px", "margin-top": "15px; margin-bottom:5px;"})
    all_samples_div2 = Div(text='<div style="padding-top: 20px;">All samples</div>', styles={"width": "500px", "margin-top": "15px; margin-bottom:5px;"})
    unweighted_panel = TabPanel(child=column([sample_div] + qbins_plots_list + [all_samples_div] + cc_plots_dict_list + [all_samples_div2] + contamination_plots_list, sizing_mode='scale_width'), title=utils_labels.QUALITY_OF_BINS)

    sample_div2 = Div(text=CLICK_ON_LEGENDS_DIV, styles={"width": "500px", "margin-top": "20px"})
    weighted_panel = TabPanel(child=column([sample_div2] + qsamples_plots_list, sizing_mode='scale_width'), title=utils_labels.QUALITY_OF_SAMPLE)
    tax_ranks_tabs = Tabs(tabs=[unweighted_panel, weighted_panel])
    return TabPanel(child=tax_ranks_tabs, title="Plots per taxonomic rank")


def create_taxonomic_binning_html(df_summary, pd_bins, labels, sample_ids_list, options):
    if pd_bins.empty:
        return None
    df_summary_t = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'taxonomic'].copy()
    df_summary_t[[utils_labels.UNIFRAC_BP, utils_labels.UNIFRAC_SEQ]] = df_summary_t[[utils_labels.UNIFRAC_BP, utils_labels.UNIFRAC_SEQ]].apply(pd.to_numeric)

    rank_to_sample_to_html = defaultdict(list)
    qbins_plots_dict = OrderedDict([(rank, '') for rank in load_ncbi_taxinfo.RANKS])
    qsamples_plots_dict = OrderedDict([(rank, '') for rank in load_ncbi_taxinfo.RANKS])
    cc_plots_dict = OrderedDict([(rank, '') for rank in load_ncbi_taxinfo.RANKS])
    contamination_plots_dict = OrderedDict([(rank, '') for rank in load_ncbi_taxinfo.RANKS])

    pd_mean = df_summary_t.groupby([utils_labels.RANK, utils_labels.TOOL]).mean(numeric_only=True).reset_index()
    pd_mean[utils_labels.SAMPLE] = AVG_OVER_SAMPLES
    for rank, pd_mean_rank in pd_mean.groupby(utils_labels.RANK):
        available_tools = list(pd_mean_rank[utils_labels.TOOL].unique())
        available_tools = [tool for tool in labels if tool in available_tools]
        pd_mean_rank = pd_mean_rank.set_index(utils_labels.TOOL).reindex(available_tools)

        purity_completeness_plot = create_precision_recall_figure(pd_mean_rank, utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP, utils_labels.AVG_PRECISION_SEQ, utils_labels.AVG_RECALL_SEQ, rank)
        purity_recall_bp_plot = create_precision_recall_figure(pd_mean_rank, utils_labels.PRECISION_PER_BP, utils_labels.RECALL_PER_BP, utils_labels.PRECISION_PER_SEQ, utils_labels.RECALL_PER_SEQ, rank)
        qbins_plots_dict[rank] = column(purity_completeness_plot)
        qsamples_plots_dict[rank] = column(purity_recall_bp_plot)

        pd_bins_rank = pd_bins[pd_bins['rank'] == rank]
        completeness_minus_contamination_plot = create_contamination_plot(pd_bins_rank, available_tools, rank + ' | Completeness - contamination', 'Index of bin (sorted by completeness - contamination (bp))', 'Completeness - contamination (bp)', plots.create_completeness_minus_contamination_column)
        contamination_plot = create_contamination_plot(pd_bins_rank, available_tools, rank + ' | Contamination', 'Index of bin (sorted by contamination (bp))', 'Contamination (bp)', plots.create_contamination_column)
        cc_plots_dict[rank] = column(completeness_minus_contamination_plot)
        contamination_plots_dict[rank] = column(contamination_plot)

        rank_to_sample_to_html[rank].append(create_table_html(pd_mean_rank.T, is_taxonomic=True))

    qbins_plots_list = [v for k, v in qbins_plots_dict.items() if v]
    qsamples_plots_list = [v for k, v in qsamples_plots_dict.items() if v]
    cc_plots_dict_list = [v for k, v in cc_plots_dict.items() if v]
    contamination_plots_list = [v for k, v in contamination_plots_dict.items() if v]

    pd_groupby_rank = df_summary_t.groupby(utils_labels.RANK)
    for rank, pd_rank in pd_groupby_rank:
        pd_rank_groupby_sample = pd_rank.groupby(utils_labels.SAMPLE)
        available_samples = pd_rank_groupby_sample.groups.keys()
        for sample_id in sample_ids_list:
            if sample_id not in available_samples:
                rank_to_sample_to_html[rank].append('')
                continue
            pd_rank_sample = pd_rank_groupby_sample.get_group(sample_id)
            pd_rank_sample = pd_rank_sample.set_index(utils_labels.TOOL).T
            rank_to_sample_to_html[rank].append(create_table_html(pd_rank_sample, is_taxonomic=True))

    if not rank_to_sample_to_html:
        return None

    sample_ids_list_combo = [AVG_OVER_SAMPLES] + sample_ids_list

    taxonomic_div = Div(text="""<div style="margin-bottom:10pt;">{}</div>""".format(rank_to_sample_to_html[load_ncbi_taxinfo.RANKS[0]][0]))
    source = ColumnDataSource(data=rank_to_sample_to_html)

    available_ranks = list(pd_groupby_rank.groups.keys())
    available_ranks_sorted = [rank for rank in load_ncbi_taxinfo.RANKS if rank in available_ranks]

    select_rank = Select(title="Taxonomic rank:", value=load_ncbi_taxinfo.RANKS[0], options=available_ranks_sorted)
    select_sample = Select(title="Sample:", value='0', options=list(zip(map(str, range(len(sample_ids_list_combo))), sample_ids_list_combo)))
    select_rank_sample_callback = CustomJS(args=dict(source=source), code="""
        mytable.text = source.data[select_rank.value][select_sample.value];
    """)
    select_rank.js_on_change('value', select_rank_sample_callback)
    select_sample.js_on_change('value', select_rank_sample_callback)
    select_rank_sample_callback.args["mytable"] = taxonomic_div
    select_rank_sample_callback.args["select_rank"] = select_rank
    select_rank_sample_callback.args["select_sample"] = select_sample

    metrics_column = column(column(row(select_sample, select_rank), create_heatmap_div(), taxonomic_div, sizing_mode='scale_width'), sizing_mode='scale_width')
    metrics_panel = TabPanel(child=metrics_column, title="Metrics")

    bins_columns = utils_labels.get_tax_bins_columns()
    if 'name' not in pd_bins.columns or pd_bins['name'].isnull().any():
        del bins_columns['name']

    metrics_bins_panel = create_metrics_per_bin_panel(pd_bins, bins_columns, sample_ids_list, options.output_dir, 'taxonomic')

    rankings_panel = TabPanel(child=column([Div(text="Click on the columns header for sorting.", styles={"width": "500px", "margin-top": "20px"}),
                                        row(create_rankings_table(pd_mean.reset_index().set_index([utils_labels.SAMPLE, utils_labels.TOOL]), True))]), title="Rankings")

    tax_ranks_panel = create_tax_ranks_panel(qbins_plots_list, qsamples_plots_list, cc_plots_dict_list, contamination_plots_list)

    tools_panel = TabPanel(child=create_plots_per_binner(pd_mean), title="Plots per binner")

    tabs = Tabs(tabs=[metrics_panel, tax_ranks_panel, tools_panel, metrics_bins_panel, rankings_panel])

    return tabs


def create_html(df_summary, pd_bins, labels, sample_ids_list, options, desc_text):
    logging.getLogger('amber').info('Creating HTML page')
    create_heatmap_bar(options.output_dir)
    tabs_list = []

    metrics_row_g = create_genome_binning_html(df_summary, pd_bins[pd_bins['rank'] == 'NA'], labels, sample_ids_list, options)
    if metrics_row_g:
        tabs_list.append(TabPanel(child=metrics_row_g, title="Genome binning"))

    metrics_row_t = create_taxonomic_binning_html(df_summary, pd_bins[pd_bins['rank'] != 'NA'], labels, sample_ids_list, options)
    if metrics_row_t:
        tabs_list.append(TabPanel(child=metrics_row_t, title="Taxonomic binning"))

    tabs = Tabs(tabs=tabs_list)

    title = create_title_div("main", "AMBER: Assessment of Metagenome BinnERs", " produced on {0} with AMBER version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))

    if desc_text:
        data_desc_div = Div(text="""<div style="text-align:left;font-size: 11pt;font-weight: bold;">{}""".format(desc_text))
        html_columns = column(title, data_desc_div, tabs, sizing_mode='scale_width')
    else:
        html_columns = column(title, tabs, sizing_mode='scale_width')

    html = file_html(models=html_columns, resources=INLINE, title='AMBER: Assessment of Metagenome BinnERs', template=TEMPLATE)

    with open(os.path.join(options.output_dir, "index.html"), 'w') as f:
        f.write(html)
