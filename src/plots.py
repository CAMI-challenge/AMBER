#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import numpy as np
import math
import re
import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src.utils import load_data
from src.utils import labels

LEGEND2 = False


def create_colors_list():
    colors_list = []
    for color in plt.cm.tab10(np.linspace(0, 1, 10))[:-1]:
        colors_list.append(tuple(color))
    colors_list.append("black")
    for color in plt.cm.Set2(np.linspace(0, 1, 8)):
        colors_list.append(tuple(color))
    for color in plt.cm.Set3(np.linspace(0, 1, 12)):
        colors_list.append(tuple(color))
    return colors_list


def load_results(files):
    table = []
    for file in files:
        f = open(file, 'r')
        field_names = f.readline().rstrip('\n').split('\t')
        field_names[field_names.index('tool')] = 'binning_label'
        for line in f:
            values = line.rstrip('\n').split('\t')
            results_dict = dict(zip(field_names, values))
            table.append(results_dict)
        f.close()
    return table


def scan_dir(output_dir):
    p_order = re.compile('^#\([0-9]+\)')
    order = []
    data_list = []
    binning_labels = []
    for path in [d for d in (os.path.join(output_dir, d1) for d1 in os.listdir(output_dir)) if os.path.isdir(d)]:
        f = open(os.path.join(path, 'purity_completeness.tsv'), 'r')
        data_list.append(load_data.load_tsv_table(f))

        # load label and order
        f = open(os.path.join(path, 'label.txt'), 'r')
        line = f.readline().rstrip('\n')
        match = p_order.match(line)
        match_string = match.group()
        order.append(int(match_string[2:match.end() - 1]))
        binning_labels.append(line[match.end():])
        f.close()

    return data_list, binning_labels, order


def plot_heatmap(df_confusion, output_dir, separate_bar=False):
    fig, axs = plt.subplots(figsize=(10, 8))

    # replace columns and rows labels by numbers
    # d = {value: key for (key, value) in enumerate(df_confusion.columns.tolist(), 1)}
    # df_confusion = df_confusion.rename(index=str, columns=d)
    # df_confusion.index = range(1, len(df_confusion.index) + 1)
    # sns_plot = sns.heatmap(df_confusion, ax=axs, annot=False, cmap="YlGnBu_r", xticklabels=3, yticklabels=3, cbar=False)

    sns_plot = sns.heatmap(df_confusion, ax=axs, annot=False, cmap="YlGnBu_r", xticklabels=False, yticklabels=False, cbar=True)
    sns_plot.set_xlabel("Genomes", fontsize=20)
    sns_plot.set_ylabel("Predicted bins", fontsize=20)
    plt.yticks(fontsize=8, rotation=0)
    plt.xticks(fontsize=8)
    # plt.yticks(fontsize=14, rotation=0)
    # plt.xticks(fontsize=14)

    fig.savefig(os.path.normpath(output_dir + '/heatmap.eps'), dpi=100, format='eps', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/heatmap.png'), dpi=100, format='png', bbox_inches='tight')
    plt.close(fig)

    if not separate_bar:
        return

    # create separate figure for bar
    fig = plt.figure(figsize=(6, 6))
    mappable = sns_plot.get_children()[0]
    fmt = lambda x, pos: '{:.0f}'.format(x / 1000000)

    cbar = plt.colorbar(mappable, orientation='vertical', label='[millions of base pairs]', format=ticker.FuncFormatter(fmt))

    text = cbar.ax.yaxis.label
    font = matplotlib.font_manager.FontProperties(size=16)
    text.set_font_properties(font)

    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=14)

    # store separate bar figure
    plt.gca().set_visible(False)
    fig.savefig(os.path.normpath(output_dir + '/heatmap_bar.eps'), dpi=100, format='eps', bbox_inches='tight')

    plt.close(fig)


def plot_boxplot(data_list, binning_labels, metric_name, output_dir, order=None):
    precision_all = []
    for metrics in data_list:
        precision = []
        for metric in metrics:
            if not math.isnan(metric[metric_name]):
                precision.append(metric[metric_name])
        precision_all.append(precision)

    if order:
        # sort binning_labels and precision_all by order
        enum_order = [(v, k) for k, v in enumerate(order)]
        enum_order = sorted(enum_order, key=lambda x: x[0])
        binning_labels = [binning_labels[i[1]] for i in enum_order]
        precision_all = [precision_all[i[1]] for i in enum_order]

    fig, axs = plt.subplots(figsize=(6, 5))

    medianprops = dict(linewidth=2.5, color='gold')
    bplot = axs.boxplot(precision_all, notch=0, vert=0, patch_artist=True, labels=binning_labels, medianprops=medianprops, sym='k.')
    colors_iter = iter(create_colors_list())

    # turn on grid
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    # force axes to be from 0 to 100%
    axs.set_xlim([-0.01, 1.01])

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}'.format(x * 100) for x in vals])

    # enable code to rotate labels
    tick_labels = axs.get_yticklabels()
    plt.setp(tick_labels, fontsize=14) ## rotation=55

    for box in bplot['boxes']:
        box.set(facecolor=next(colors_iter), linewidth=0.1)
    plt.ylim(plt.ylim()[::-1])

    if metric_name == 'purity':
        axs.set_xlabel('Purity per bin $p$ (%)' if LEGEND2 else 'Purity per bin (%)', fontsize=14)
    else:
        axs.set_xlabel('Completeness per genome $r$ (%)' if LEGEND2 else 'Completeness per genome (%)', fontsize=14)

    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '.png'), dpi=100, format='png', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')

    # remove labels but keep grid
    axs.get_yaxis().set_ticklabels([])
    for tic in axs.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False
    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '_wo_legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


def plot_summary(df_results, output_dir, plot_type, file_name, xlabel, ylabel):
    colors_list = create_colors_list()
    df_groups = df_results.groupby(labels.TOOL)

    if len(df_groups) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axes to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    plot_labels = []
    if plot_type == 'e':
        for i, (tool, pd_summary) in enumerate(df_groups):
            axs.errorbar(float(pd_summary[labels.AVG_PRECISION][0]), float(pd_summary[labels.AVG_RECALL][0]), xerr=float(pd_summary[labels.AVG_PRECISION_SEM][0]), yerr=float(pd_summary[labels.AVG_RECALL_SEM][0]),
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
            plot_labels.append(tool)
    if plot_type == 'w':
        for i, (tool, pd_summary) in enumerate(df_groups):
            axs.plot(float(pd_summary[labels.AVG_PRECISION_PER_BP][0]), float(pd_summary[labels.AVG_RECALL_PER_BP][0]), marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(tool)
    elif plot_type == 'p':
        for i, (tool, pd_summary) in enumerate(df_groups):
            axs.plot(float(pd_summary[labels.ARI_BY_BP][0]), float(pd_summary[labels.PERCENTAGE_ASSIGNED_BPS][0]), marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(tool)

    # turn on grid
    axs.minorticks_on()
    axs.grid(which='major', linestyle='-', linewidth='0.5')
    axs.grid(which='minor', linestyle=':', linewidth='0.5')

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}'.format(x * 100) for x in vals])
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}'.format(x * 100) for x in vals])

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.tight_layout()
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')

    colors_iter = iter(colors_list)
    circles = []
    for x in range(len(df_groups)):
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=11, markerfacecolor=next(colors_iter)))

    lgd = plt.legend(circles, plot_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False, fontsize=12)

    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def plot_avg_precision_recall(df_results_g, output_dir):
    plot_summary(df_results_g,
                 output_dir,
                 'e',
                 'avg_purity_completeness',
                 'Truncated average purity per bin $\overline{p}_{99}$ (%)' if LEGEND2 else 'Average purity per bin (%)',
                 'Average completeness per genome $\overline{r}$ (%)' if LEGEND2 else 'Average completeness per genome (%)')


def plot_weighed_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'w',
                 'avg_purity_completeness_per_bp',
                 'Average purity per base pair $\overline{p}_{bp}$ (%)' if LEGEND2 else 'Average purity per base pair (%)',
                 'Average completeness per base pair $\overline{r}_{bp}$ (%)' if LEGEND2 else 'Average completeness per base pair (%)')


def plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'p',
                 'ari_vs_assigned_bps',
                 'Adjusted Rand Index (%)' if LEGEND2 else 'Adjusted Rand Index',
                 'Percentage of assigned base pairs (%)' if LEGEND2 else 'Percentage of assigned base pairs')
