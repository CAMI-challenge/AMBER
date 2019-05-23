#!/usr/bin/env python

from src.utils import load_data
from src.utils import labels as utils_labels
from src.utils import load_ncbi_taxinfo
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import numpy as np
import re
import os, sys, inspect
from collections import OrderedDict
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)


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


def create_legend(df_results, output_dir):
    colors_iter = iter(create_colors_list())
    labels = list(df_results.groupby(utils_labels.TOOL).groups.keys())
    circles = [Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=10, markerfacecolor=next(colors_iter)) for label in labels]

    fig = plt.figure(figsize=(0.5, 0.5))
    fig.legend(circles, labels, loc='center', frameon=False, ncol=5, handletextpad=0.1)
    fig.savefig(os.path.join(output_dir, 'genome', 'legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


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


def plot_heatmap(df_confusion, sample_id, output_dir, separate_bar=False):
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

    fig.savefig(os.path.join(output_dir, 'heatmap_' + sample_id + '.eps'), dpi=100, format='eps', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'heatmap_' + sample_id + '.png'), dpi=100, format='png', bbox_inches='tight')
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


def plot_boxplot(pd_bins, metric_name, output_dir):
    metric_all = []
    binning_labels = []
    for tool, pd_bins_tool in pd_bins.groupby(utils_labels.TOOL):
        binning_labels.append(tool)
        metric_all.append(pd_bins_tool[metric_name][pd_bins_tool[metric_name].notnull()].tolist())

    fig, axs = plt.subplots(figsize=(6, 5))

    medianprops = dict(linewidth=2.5, color='gold')
    bplot = axs.boxplot(metric_all, notch=0, vert=0, patch_artist=True, labels=binning_labels, medianprops=medianprops, sym='k.')
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
        axs.set_xlabel('Purity per bin (%)', fontsize=14)
    else:
        axs.set_xlabel('Completeness per genome (%)', fontsize=14)

    fig.savefig(os.path.join(output_dir, 'genome', 'boxplot_' + metric_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'genome', 'boxplot_' + metric_name + '.png'), dpi=100, format='png', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'genome', 'boxplot_' + metric_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')

    # remove labels but keep grid
    axs.get_yaxis().set_ticklabels([])
    for tic in axs.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False
    fig.savefig(os.path.join(output_dir, 'genome', 'boxplot_' + metric_name + '_wo_legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


def plot_summary(df_results, output_dir, rank, plot_type, file_name, xlabel, ylabel):
    colors_list = create_colors_list()
    df_groups = df_results.groupby(utils_labels.TOOL)
    binning_type = df_results[utils_labels.BINNING_TYPE].iloc[0]

    if len(df_groups) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axes to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    if plot_type == 'e':
        for i, (tool, pd_summary) in enumerate(df_groups):
            df_mean = pd_summary[[utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP, utils_labels.AVG_PRECISION_BP_SEM, utils_labels.AVG_RECALL_BP_SEM]].mean()
            axs.errorbar(df_mean[utils_labels.AVG_PRECISION_BP], df_mean[utils_labels.AVG_RECALL_BP], xerr=df_mean[utils_labels.AVG_PRECISION_BP_SEM], yerr=df_mean[utils_labels.AVG_RECALL_BP_SEM],
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
    if plot_type == 'f':
        for i, (tool, pd_summary) in enumerate(df_groups):
            df_mean = pd_summary[[utils_labels.AVG_PRECISION_SEQ, utils_labels.AVG_RECALL_SEQ, utils_labels.AVG_PRECISION_SEQ_SEM, utils_labels.AVG_RECALL_SEQ_SEM]].mean()
            axs.errorbar(df_mean[utils_labels.AVG_PRECISION_SEQ], df_mean[utils_labels.AVG_RECALL_SEQ], xerr=df_mean[utils_labels.AVG_PRECISION_SEQ_SEM], yerr=df_mean[utils_labels.AVG_RECALL_SEQ_SEM],
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
    if plot_type == 'w':
        for i, (tool, pd_summary) in enumerate(df_groups):
            df_mean = pd_summary[[utils_labels.PRECISION_PER_BP, utils_labels.RECALL_PER_BP]].mean()
            axs.plot(df_mean[utils_labels.PRECISION_PER_BP], df_mean[utils_labels.RECALL_PER_BP], marker='o', color=colors_list[i], markersize=10)
    elif plot_type == 'p':
        for i, (tool, pd_summary) in enumerate(df_groups):
            df_mean = pd_summary[[utils_labels.ARI_BY_BP, utils_labels.PERCENTAGE_ASSIGNED_BPS]].mean()
            axs.plot(df_mean[utils_labels.ARI_BY_BP], df_mean[utils_labels.PERCENTAGE_ASSIGNED_BPS], marker='o', color=colors_list[i], markersize=10)

    # turn on grid
    axs.minorticks_on()
    axs.grid(which='major', linestyle='-', linewidth='0.5')
    axs.grid(which='minor', linestyle=':', linewidth='0.5')

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}'.format(x * 100) for x in vals])
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}'.format(x * 100) for x in vals])

    if rank:
        file_name = rank + '_' + file_name
        plt.title(rank)
        ylabel = ylabel.replace('genome', 'taxon')

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')

    colors_iter = iter(colors_list)
    circles = []
    for x in range(len(df_groups)):
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=11, markerfacecolor=next(colors_iter)))

    lgd = plt.legend(circles, list(df_groups.groups.keys()), bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False, fontsize=12)

    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def plot_avg_precision_recall(df_results, output_dir, rank=None):
    plot_summary(df_results,
                 output_dir,
                 rank,
                 'e',
                 'avg_purity_completeness_bp',
                 'Average purity per bin (%)',
                 'Average completeness per genome (%)')
    plot_summary(df_results,
                 output_dir,
                 rank,
                 'f',
                 'avg_purity_completeness_seq',
                 'Average purity per bin (%)',
                 'Average completeness per genome (%)')


def plot_precision_recall(summary_per_query, output_dir, rank=None):
    plot_summary(summary_per_query,
                 output_dir,
                 rank,
                 'w',
                 'purity_completeness_bp',
                 'Purity per base pair (%)',
                 'Completeness per base pair (%)')
    plot_summary(summary_per_query,
                 output_dir,
                 rank,
                 'x',
                 'purity_completeness_seq',
                 'Purity per sequence (%)',
                 'Completeness per sequence (%)')


def plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, output_dir, rank=None):
    plot_summary(summary_per_query,
                 output_dir,
                 rank,
                 'p',
                 'ari_vs_assigned_bps',
                 'Adjusted Rand Index (%)',
                 'Percentage of assigned base pairs (%)')


def plot_taxonomic_results(df_summary_t, metrics_list, errors_list, file_name, output_dir):
    colors_list = ["#006cba", "#008000", "#ba9e00", "red"]

    for tool, pd_results in df_summary_t.groupby(utils_labels.TOOL):
        dict_metric_list = []
        for metric in metrics_list:
            rank_to_metric = OrderedDict([(k, .0) for k in load_ncbi_taxinfo.RANKS])
            dict_metric_list.append(rank_to_metric)
        dict_error_list = []
        for error in errors_list:
            rank_to_metric_error = OrderedDict([(k, .0) for k in load_ncbi_taxinfo.RANKS])
            dict_error_list.append(rank_to_metric_error)

        for index, row in pd_results.iterrows():
            for rank_to_metric, metric in zip(dict_metric_list, metrics_list):
                rank_to_metric[row[utils_labels.RANK]] = .0 if np.isnan(row[metric]) else row[metric]
            for rank_to_metric_error, error in zip(dict_error_list, errors_list):
                rank_to_metric_error[row[utils_labels.RANK]] = .0 if np.isnan(row[error]) else row[error]

        fig, axs = plt.subplots(figsize=(6, 5))

        # force axes to be from 0 to 100%
        axs.set_xlim([0, 7])
        axs.set_ylim([0.0, 1.0])
        x_values = range(len(load_ncbi_taxinfo.RANKS))

        y_values_list = []
        for rank_to_metric, color in zip(dict_metric_list, colors_list):
            y_values = list(rank_to_metric.values())
            axs.plot(x_values, y_values, color=color)
            y_values_list.append(y_values)

        for rank_to_metric_error, y_values, color in zip(dict_error_list, y_values_list, colors_list):
            sem = list(rank_to_metric_error.values())
            plt.fill_between(x_values, np.subtract(y_values, sem).tolist(), np.add(y_values, sem).tolist(), color=color, alpha=0.5)

        plt.xticks(x_values, load_ncbi_taxinfo.RANKS, rotation='vertical')

        vals = axs.get_yticks()
        axs.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

        lgd = plt.legend(metrics_list, loc=1, borderaxespad=0., handlelength=2, frameon=False)

        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, 'taxonomic', tool, file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
        fig.savefig(os.path.join(output_dir, 'taxonomic', tool, file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close(fig)


def create_contamination_column(pd_tool_bins):
    pd_tool_bins['newcolumn'] = 1 - pd_tool_bins['purity_bp']


def create_completeness_minus_contamination_column(pd_tool_bins):
    pd_tool_bins['newcolumn'] = pd_tool_bins['completeness_bp'] + pd_tool_bins['purity_bp'] - 1


def plot_contamination(pd_bins, binning_type, title, xlabel, ylabel, create_column_function, output_dir):
    if len(pd_bins) == 0:
        return

    pd_bins_copy = pd_bins[[utils_labels.TOOL, 'purity_bp', 'completeness_bp']].copy().dropna(subset=['purity_bp'])
    create_column_function(pd_bins_copy)

    colors_list = create_colors_list()

    fig, axs = plt.subplots(figsize=(6, 5))

    tools = pd_bins_copy[utils_labels.TOOL].unique().tolist()

    for color, tool in zip(colors_list, tools):
        pd_tool_bins = pd_bins_copy[pd_bins_copy[utils_labels.TOOL] == tool]
        pd_tool_bins = pd_tool_bins.sort_values(by='newcolumn', ascending=False).reset_index()
        pd_tool_bins = pd_tool_bins.drop(['index'], axis=1)

        axs.plot(list(range(1, len(pd_tool_bins) + 1)), pd_tool_bins['newcolumn'], color=color)

    axs.set_ylim(pd_bins_copy['newcolumn'].min(), 1.0)
    axs.set_xlim(1, None)
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    # transform plot_labels to percentages
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}'.format(y * 100) for y in vals])

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel + ' (%)', fontsize=14)

    lgd = plt.legend(tools, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=1, frameon=False, fontsize=12)

    plt.tight_layout()

    file_name = title.lower().replace(' ', '_').replace('-', 'minus').replace('|', '')
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)
