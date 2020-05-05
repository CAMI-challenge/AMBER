# Copyright 2020 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
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

from src.utils import labels as utils_labels
from src.utils import load_ncbi_taxinfo
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import numpy as np
import os, sys, inspect
import pandas as pd
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


def create_legend(color_indices, available_tools, output_dir):
    colors_list = create_colors_list()
    if color_indices:
        colors_list = [colors_list[i] for i in color_indices]

    colors_iter = iter(colors_list)
    circles = [Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=10, markerfacecolor=next(colors_iter)) for label in available_tools]

    fig = plt.figure(figsize=(0.5, 0.5))
    fig.legend(circles, available_tools, loc='center', frameon=False, ncol=5, handletextpad=0.1)
    fig.savefig(os.path.join(output_dir, 'genome', 'legend.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    plt.close(fig)


def plot_precision_vs_bin_size(pd_bins, output_dir):
    pd_plot = pd_bins[pd_bins[utils_labels.TOOL] != utils_labels.GS]

    for tool_label, pd_tool in pd_plot.groupby(utils_labels.TOOL):
        fig, axs = plt.subplots(figsize=(5, 4.5))
        axs.scatter(np.log(pd_tool['total_length']), pd_tool['precision_bp'], marker='o')

        axs.set_xlim([None, np.log(pd_tool['total_length'].max())])
        axs.set_ylim([0.0, 1.0])
        axs.set_title(tool_label, fontsize=12)

        plt.ylabel('Purity per bin (%)', fontsize=12)
        plt.xlabel('Bin size [log(# bp)]', fontsize=12)

        fig.savefig(os.path.join(output_dir, 'genome', tool_label, 'purity_vs_bin_size.png'), dpi=200, format='png', bbox_inches='tight')
        plt.close(fig)


def plot_by_genome_coverage(pd_bins, pd_target_column, available_tools, output_dir):
    colors_list = create_colors_list()
    if len(available_tools) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(5, 4.5))

    for i, (color, tool) in enumerate(zip(colors_list, available_tools)):
        pd_tool = pd_bins[pd_bins[utils_labels.TOOL] == tool].sort_values(by=['genome_index'])
        axs.scatter(pd_tool['genome_coverage'], pd_tool[pd_target_column], marker='o', color=colors_list[i], s=[2] * pd_tool.shape[0])
        window = 50
        rolling_mean = pd_tool[pd_target_column].rolling(window=window, min_periods=10).mean()
        axs.plot(pd_tool['genome_coverage'], rolling_mean, color=colors_list[i])

    axs.set_xlim([0.0, pd_tool['genome_coverage'].max()])
    axs.set_ylim([0.0, 1.0])

    # transform plot_labels to percentages
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}'.format(x * 100) for x in vals], fontsize=12)

    axs.tick_params(axis='x', labelsize=12)

    if pd_target_column == 'precision_bp':
        ylabel = 'Purity per bin (%)'
        file_name = 'purity_by_genome_coverage'
    else:
        ylabel = 'Completeness per genome (%)'
        file_name = 'completeness_by_genome_coverage'

    plt.ylabel(ylabel, fontsize=15)
    plt.xlabel('log$_{10}$(average genome coverage)', fontsize=15)

    colors_iter = iter(colors_list)
    circles = []
    for x in range(len(available_tools)):
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=11, markerfacecolor=next(colors_iter)))
    lgd = plt.legend(circles, available_tools, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False, fontsize=14)

    fig.savefig(os.path.join(output_dir, 'genome', file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


# def get_genome_recall_all_samples_df(sample_id_to_queries_list):
#     pd_query_all = pd.DataFrame()
#     for sample_id in sample_id_to_queries_list:
#         for query in sample_id_to_queries_list[sample_id]:
#             pd_query = pd.DataFrame.from_dict(query.genome_to_recall_bp, orient='index', columns=['recall_bp'])
#             pd_query.index.name = 'genome_id'
#             pd_query['sample_id'] = sample_id
#             pd_query[utils_labels.TOOL] = query.label
#             pd_query = pd_query.reset_index().set_index(['sample_id', utils_labels.TOOL])
#             pd_query_all = pd.concat([pd_query_all, pd_query])
#     return pd_query_all


def get_pd_genomes_recall(sample_id_to_queries_list):
    pd_genomes_recall = pd.DataFrame()
    for sample_id in sample_id_to_queries_list:
        for query in sample_id_to_queries_list[sample_id]:
            pd_query = pd.DataFrame.from_dict(query.genome_to_recall_bp, orient='index', columns=['recall_bp'])
            pd_query.index.name = 'genome_id'
            pd_query['sample_id'] = sample_id
            pd_query[utils_labels.TOOL] = query.label
            pd_query = pd_query.reset_index().set_index(['sample_id', utils_labels.TOOL])
            pd_genomes_recall = pd.concat([pd_genomes_recall, pd_query])
    return pd_genomes_recall


def plot_precision_recall_by_coverage(sample_id_to_queries_list, pd_bins_g, coverages_pd, available_tools, output_dir):
    # compute average genome coverage if coverages for multiple samples were provided
    coverages_pd['mean_coverage'] = coverages_pd.mean(axis=1)
    coverages_pd = coverages_pd.sort_values(by=['mean_coverage'])
    coverages_pd['rank'] = coverages_pd['mean_coverage'].rank()

    pd_genomes_recall = get_pd_genomes_recall(sample_id_to_queries_list)
    pd_genomes_recall['genome_index'] = pd_genomes_recall['genome_id'].map(coverages_pd['rank'].to_dict())
    pd_genomes_recall = pd_genomes_recall.groupby([utils_labels.TOOL, 'genome_id']).mean().reset_index()

    pd_genomes_recall['genome_coverage'] = np.log10(pd_genomes_recall['genome_id'].map(coverages_pd['mean_coverage'].to_dict()))

    # print(pd_genomes_recall)
    # pd_genomes_recall['genome_coverage_log10'] = np.log10(pd_genomes_recall['genome_id'].map(coverages_pd['mean_coverage'].to_dict()))
    # pd_genomes_recall['genome_coverage'] = pd_genomes_recall['genome_id'].map(coverages_pd['mean_coverage'].to_dict())
    # print(pd_genomes_recall)
    # pd_genomes_recall.to_csv('/home/fmeyer/tutorial/amber_mouse_gut/coverages_recall_pd.tsv', sep='\t')
    # exit()

    plot_by_genome_coverage(pd_genomes_recall, 'recall_bp', available_tools, output_dir)

    pd_bins_precision = pd_bins_g[[utils_labels.TOOL, 'precision_bp', 'most_abundant_genome']].copy().dropna(subset=['precision_bp'])
    pd_bins_precision['genome_index'] = pd_bins_precision['most_abundant_genome'].map(coverages_pd['rank'].to_dict())
    pd_bins_precision['genome_coverage'] = np.log10(pd_bins_precision['most_abundant_genome'].map(coverages_pd['mean_coverage'].to_dict()))
    plot_by_genome_coverage(pd_bins_precision, 'precision_bp', available_tools, output_dir)


def plot_heatmap(df_confusion, sample_id, output_dir, label, separate_bar=False, log_scale=False):
    if log_scale:
        df_confusion = df_confusion.apply(np.log10, inplace=True).replace(-np.inf, 0)
    fig, axs = plt.subplots(figsize=(10, 8))
    fontsize = 20

    # replace columns and rows labels by numbers
    d = {value: key for (key, value) in enumerate(df_confusion.columns.tolist(), 1)}
    df_confusion = df_confusion.rename(index=str, columns=d)
    df_confusion.index = range(1, len(df_confusion) + 1)

    xticklabels = int(round(df_confusion.shape[1] / 10, -1))
    yticklabels = int(round(df_confusion.shape[0] / 10, -1))
    sns_plot = sns.heatmap(df_confusion, ax=axs, annot=False, cmap="YlGnBu_r", xticklabels=xticklabels, yticklabels=yticklabels, cbar=False, rasterized=True)

    # sns_plot = sns.heatmap(df_confusion, ax=axs, annot=False, cmap="YlGnBu_r", xticklabels=False, yticklabels=False, cbar=True, rasterized=True)
    sns_plot.set_xlabel("Genomes", fontsize=fontsize)
    sns_plot.set_ylabel("Predicted bins", fontsize=fontsize)
    plt.yticks(fontsize=12, rotation=0)
    plt.xticks(fontsize=12)

    mappable = sns_plot.get_children()[0]

    cbar_ax = fig.add_axes([.915, .11, .017, .77])
    cbar = plt.colorbar(mappable, ax=axs, cax=cbar_ax, orientation='vertical')
    if log_scale:
        cbar.set_label(fontsize=fontsize, label='log$_{10}$(# bp)')
    else:
        fmt = lambda x, pos: '{:.0f}'.format(x / 1000000)
        cbar = plt.colorbar(mappable, ax=axs, cax=cbar_ax, orientation='vertical', label='Millions of base pairs', format=ticker.FuncFormatter(fmt))

    cbar.set_label(fontsize=fontsize, label='# bp')
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.outline.set_edgecolor(None)

    axs.set_title(label, fontsize=fontsize, pad=10)

    axs.set_ylim([len(df_confusion), 0])

    # plt.yticks(fontsize=14, rotation=0)
    # plt.xticks(fontsize=14)

    output_dir = os.path.join(output_dir, 'genome', label)

    fig.savefig(os.path.join(output_dir, 'heatmap_' + sample_id + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'heatmap_' + sample_id + '.png'), dpi=200, format='png', bbox_inches='tight')
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
    fig.savefig(os.path.join(output_dir, 'heatmap_bar.pdf'), dpi=100, format='pdf', bbox_inches='tight')

    plt.close(fig)


def plot_boxplot(sample_id_to_queries_list, metric_name, output_dir, available_tools):
    pd_bins = pd.DataFrame()
    for sample_id in sample_id_to_queries_list:
        for query in sample_id_to_queries_list[sample_id]:
            metric_df = getattr(query, metric_name.replace('_bp', '_df')).copy()
            metric_df[utils_labels.TOOL] = query.label
            metric_df['sample_id'] = sample_id
            metric_df = metric_df.reset_index().set_index(['sample_id', utils_labels.TOOL])
            pd_bins = pd.concat([pd_bins, metric_df])

    metric_all = []

    for tool in available_tools:
        pd_tool = pd_bins.iloc[pd_bins.index.get_level_values(utils_labels.TOOL) == tool]
        metric_all.append(pd_tool[metric_name][pd_tool[metric_name].notnull()].tolist())

    fig, axs = plt.subplots(figsize=(6, 5))

    medianprops = dict(linewidth=2.5, color='gold')
    bplot = axs.boxplot(metric_all, notch=0, vert=0, patch_artist=True, labels=available_tools, medianprops=medianprops, sym='k.')
    colors_iter = iter(create_colors_list())

    # turn on grid
    axs.grid(which='major', linestyle=':', linewidth='0.5', color='lightgrey')

    # force axes to be from 0 to 100%
    axs.set_xlim([-0.01, 1.01])

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}'.format(x * 100) for x in vals])

    # enable code to rotate labels
    tick_labels = axs.get_yticklabels()
    plt.setp(tick_labels, fontsize=13) ## rotation=55

    for box in bplot['boxes']:
        box.set(facecolor=next(colors_iter), linewidth=0.1)
    plt.ylim(plt.ylim()[::-1])

    if metric_name == 'precision_bp':
        axs.set_xlabel('Purity per bin (%)', fontsize=13)
        metric_name = 'purity_bp'
    else:
        axs.set_xlabel('Completeness per genome (%)', fontsize=13)
        metric_name = 'completeness_bp'

    fig.savefig(os.path.join(output_dir, 'genome', 'boxplot_' + metric_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'genome', 'boxplot_' + metric_name + '.png'), dpi=200, format='png', bbox_inches='tight')

    # remove labels but keep grid
    # axs.get_yaxis().set_ticklabels([])
    # for tic in axs.yaxis.get_major_ticks():
    #     tic.tick1line.set_visible(False)
    #     tic.tick2line.set_visible(False)
    #     tic.label1.set_visible(False)
    #     tic.label2.set_visible(False)
    # fig.savefig(os.path.join(output_dir, 'genome', 'boxplot_' + metric_name + '_wo_legend.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    plt.close(fig)


def plot_summary(color_indices, df_results, labels, output_dir, rank, plot_type, file_name, xlabel, ylabel):
    available_tools = df_results[utils_labels.TOOL].unique()
    tools = [tool for tool in labels if tool in available_tools]

    colors_list = create_colors_list()
    if color_indices:
        colors_list = [colors_list[i] for i in color_indices]
    df_mean = df_results.groupby(utils_labels.TOOL).mean().reindex(tools)

    binning_type = df_results[utils_labels.BINNING_TYPE].iloc[0]

    if len(df_mean) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(5, 4.5))

    # force axes to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    if plot_type == 'e':
        for i, (tool, df_row) in enumerate(df_mean.iterrows()):
            axs.errorbar(df_row[utils_labels.AVG_PRECISION_BP], df_row[utils_labels.AVG_RECALL_BP], xerr=df_row['avg_precision_bp_var'], yerr=df_row['avg_recall_bp_var'],
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
    if plot_type == 'f':
        for i, (tool, df_row) in enumerate(df_mean.iterrows()):
            axs.errorbar(df_row[utils_labels.AVG_PRECISION_SEQ], df_row[utils_labels.AVG_RECALL_SEQ], xerr=df_row[utils_labels.AVG_PRECISION_SEQ_SEM], yerr=df_row[utils_labels.AVG_RECALL_SEQ_SEM],
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
    if plot_type == 'w':
        for i, (tool, df_row) in enumerate(df_mean.iterrows()):
            axs.plot(df_row[utils_labels.PRECISION_PER_BP], df_row[utils_labels.RECALL_PER_BP], marker='o', color=colors_list[i], markersize=10)
    if plot_type == 'x':
        for i, (tool, df_row) in enumerate(df_mean.iterrows()):
            axs.plot(df_row[utils_labels.PRECISION_PER_SEQ], df_row[utils_labels.RECALL_PER_SEQ], marker='o', color=colors_list[i], markersize=10)
    elif plot_type == 'p':
        for i, (tool, df_row) in enumerate(df_mean.iterrows()):
            axs.plot(df_row[utils_labels.ARI_BY_BP], df_row[utils_labels.PERCENTAGE_ASSIGNED_BPS], marker='o', color=colors_list[i], markersize=10)

    # turn on grid
    # axs.minorticks_on()
    axs.grid(which='major', linestyle=':', linewidth='0.5')
    # axs.grid(which='minor', linestyle=':', linewidth='0.5')

    # transform plot_labels to percentages
    if plot_type != 'p':
        vals = axs.get_xticks()
        axs.set_xticklabels(['{:3.0f}'.format(x * 100) for x in vals], fontsize=11)
    else:
        axs.tick_params(axis='x', labelsize=12)
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}'.format(x * 100) for x in vals], fontsize=11)

    if rank:
        file_name = rank + '_' + file_name
        plt.title(rank)
        ylabel = ylabel.replace('genome', 'taxon')

    plt.xlabel(xlabel, fontsize=13)
    plt.ylabel(ylabel, fontsize=13)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')

    colors_iter = iter(colors_list)
    circles = []
    for x in range(len(df_mean)):
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=11, markerfacecolor=next(colors_iter)))
    lgd = plt.legend(circles, tools, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False, fontsize=12)

    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.png'), dpi=200, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def plot_avg_precision_recall(colors, df_results, labels, output_dir, rank=None):
    plot_summary(colors,
                 df_results,
                 labels,
                 output_dir,
                 rank,
                 'e',
                 'avg_purity_recall_bp',
                 'Average purity per bin (%)',
                 'Average completeness per genome (%)')
    plot_summary(colors,
                 df_results,
                 labels,
                 output_dir,
                 rank,
                 'f',
                 'avg_purity_completeness_seq',
                 'Average purity per bin (%)',
                 'Average completeness per genome (%)')


def plot_precision_recall(colors, summary_per_query, labels, output_dir, rank=None):
    plot_summary(colors,
                 summary_per_query,
                 labels,
                 output_dir,
                 rank,
                 'w',
                 'purity_recall_bp',
                 'Purity for sample (%)',
                 'Completeness for sample (%)')
    plot_summary(colors,
                 summary_per_query,
                 labels,
                 output_dir,
                 rank,
                 'x',
                 'purity_completeness_seq',
                 'Purity for sample (%)',
                 'Completeness for sample (%)')


def plot_adjusted_rand_index_vs_assigned_bps(colors, summary_per_query, labels, output_dir, rank=None):
    plot_summary(colors,
                 summary_per_query,
                 labels,
                 output_dir,
                 rank,
                 'p',
                 'ari_vs_assigned_bps',
                 'Adjusted Rand index',
                 'Percentage of binned base pairs')


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
    pd_tool_bins['newcolumn'] = 1 - pd_tool_bins['precision_bp']


def create_completeness_minus_contamination_column(pd_tool_bins):
    pd_tool_bins['newcolumn'] = pd_tool_bins['recall_bp'] + pd_tool_bins['precision_bp'] - 1


def plot_contamination(pd_bins, binning_type, title, xlabel, ylabel, create_column_function, output_dir):
    if len(pd_bins) == 0:
        return

    pd_bins_copy = pd_bins[[utils_labels.TOOL, 'precision_bp', 'recall_bp']].copy().dropna(subset=['precision_bp'])
    create_column_function(pd_bins_copy)

    colors_list = create_colors_list()

    fig, axs = plt.subplots(figsize=(6, 5))

    tools = pd_bins_copy[utils_labels.TOOL].unique().tolist()

    for color, tool in zip(colors_list, tools):
        pd_tool_bins = pd_bins_copy[pd_bins_copy[utils_labels.TOOL] == tool]
        pd_tool_bins = pd_tool_bins.sort_values(by='newcolumn', ascending=False).reset_index()
        pd_tool_bins = pd_tool_bins.drop(['index'], axis=1)

        axs.plot(list(range(1, len(pd_tool_bins) + 1)), pd_tool_bins['newcolumn'], color=color)

    min_value = pd_bins_copy['newcolumn'].min()
    axs.set_ylim(min_value if min_value < 1.0 else .9, 1.0)
    axs.set_xlim(1, None)
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    # transform plot_labels to percentages
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}'.format(y * 100) for y in vals])

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel + ' [%]', fontsize=14)

    lgd = plt.legend(tools, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=1, frameon=False, fontsize=12)

    plt.tight_layout()

    file_name = title.lower().replace(' ', '_').replace('-', 'minus').replace('|', '')
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, binning_type, file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def get_number_of_hq_bins(tools, pd_bins):
    pd_counts = pd.DataFrame()
    pd_bins_copy = pd_bins[[utils_labels.TOOL, 'precision_bp', 'recall_bp']].copy().dropna(subset=['precision_bp'])
    for tool in tools:
        pd_tool_bins = pd_bins_copy[pd_bins_copy[utils_labels.TOOL] == tool]
        x50 = pd_tool_bins[(pd_tool_bins['recall_bp'] > .5) & (pd_tool_bins['precision_bp'] > .9)].shape[0]
        x70 = pd_tool_bins[(pd_tool_bins['recall_bp'] > .7) & (pd_tool_bins['precision_bp'] > .9)].shape[0]
        x90 = pd_tool_bins[(pd_tool_bins['recall_bp'] > .9) & (pd_tool_bins['precision_bp'] > .9)].shape[0]
        pd_tool_counts = pd.DataFrame([[x90, x70, x50]], columns=['>90%', '>70%', '>50%'], index=[tool])
        pd_counts = pd_counts.append(pd_tool_counts)
    return pd_counts


def get_number_of_hq_bins_by_score(tools, pd_bins):
    pd_counts = pd.DataFrame()
    pd_bins_copy = pd_bins[[utils_labels.TOOL, 'precision_bp', 'recall_bp']].copy().dropna(subset=['precision_bp'])
    pd_bins_copy['newcolumn'] = pd_bins_copy['recall_bp'] + 5 * (pd_bins_copy['precision_bp'] - 1)
    for tool in tools:
        pd_tool_bins = pd_bins_copy[pd_bins_copy[utils_labels.TOOL] == tool]
        x50 = pd_tool_bins[pd_tool_bins['newcolumn'] > .5].shape[0]
        x70 = pd_tool_bins[pd_tool_bins['newcolumn'] > .7].shape[0]
        x90 = pd_tool_bins[pd_tool_bins['newcolumn'] > .9].shape[0]
        x50 -= x70
        x70 -= x90
        pd_tool_counts = pd.DataFrame([[x90, x70, x50]], columns=['>90', '>70', '>50'], index=[tool])
        pd_counts = pd_counts.append(pd_tool_counts)
    return pd_counts


def plot_counts(pd_bins, tools, output_dir, output_file, get_bin_counts_function):
    pd_counts = get_bin_counts_function(tools, pd_bins)
    fig, axs = plt.subplots(figsize=(11, 5))
    if output_file == 'bin_counts':
        fig = pd_counts.plot.bar(ax=axs, stacked=False, color=['#28334AFF', '#FBDE44FF', '#F65058FF'], width=.8, legend=None).get_figure()
    else:
        fig = pd_counts.plot.bar(ax=axs, stacked=True, color=['#9B4A97FF', '#FC766AFF', '#F9A12EFF'], width=.8, legend=None).get_figure()

    axs.tick_params(axis='x', labelrotation=45, length=0)
    axs.set_xticklabels(tools, horizontalalignment='right', fontsize=14)
    axs.set_xlabel(None)

    # axs.yaxis.set_major_locator(MaxNLocator(integer=True))

    h, l = axs.get_legend_handles_labels()
    axs.set_ylabel('#genome bins', fontsize=14)

    # axs.grid(which='major', linestyle=':', linewidth='0.5')
    # axs.grid(which='minor', linestyle=':', linewidth='0.5')

    ph = [plt.plot([], marker='', ls='')[0]]
    handles = ph + h

    if output_file == 'bin_counts':
        labels = ['Contamination < 10%           Completeness  '] + l
        bbox_to_anchor = (0.49, 1.02)
    else:
        labels = ['Score  '] + l
        y_values = (pd_counts['>90'] + pd_counts['>70'] + pd_counts['>50']).tolist()
        for i, v in enumerate(y_values):
            axs.text(i - .25, v + 5, str(v), color='black', fontweight='bold')
        bbox_to_anchor = (0.47, 1.02)

    lgd = plt.legend(handles, labels, bbox_to_anchor=bbox_to_anchor, columnspacing=.5, loc=8, borderaxespad=0., handlelength=1, frameon=False, fontsize=14, ncol=5)

    # plt.subplots_adjust(hspace=0.6, wspace=0.2)

    fig.savefig(os.path.join(output_dir, 'genome', output_file + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'genome', output_file + '.png'), dpi=200, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)
