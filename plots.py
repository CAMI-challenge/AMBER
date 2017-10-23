#!/usr/bin/env python

import argparse
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import re
from utils import load_data


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
        f = open(path + '/precision_recall.tsv', 'r')
        data_list.append(load_data.load_tsv_table(f))

        # load label and order
        f = open(path + '/label.txt', 'r')
        line = f.readline().rstrip('\n')
        match = p_order.match(line)
        match_string = match.group()
        order.append(int(match_string[2:match.end() - 1]))
        binning_labels.append(line[match.end():])
        f.close()

    return data_list, binning_labels, order


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
    axs.set_xticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    for box in bplot['boxes']:
        box.set(facecolor=next(colors_iter), linewidth=0.1)
    plt.ylim(plt.ylim()[::-1])

    if metric_name == 'precision':
        axs.set_xlabel(metric_name.title() + ' per bin')
    else:
        axs.set_xlabel(metric_name.title() + ' per genome')

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


def plot_summary(summary_per_query, output_dir, plot_type, file_name, xlabel, ylabel):
    colors_list = create_colors_list()
    if len(summary_per_query) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axes to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    i = 0
    plot_labels = []
    if plot_type == 'e':
        for summary in summary_per_query:
            axs.errorbar(float(summary['avg_precision']), float(summary['avg_recall']), xerr=float(summary['sem_precision']), yerr=float(summary['sem_recall']),
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
            plot_labels.append(summary['binning_label'])
            i += 1
    if plot_type == 'w':
        for summary in summary_per_query:
            axs.plot(float(summary['avg_precision_per_bp']), float(summary['avg_recall_per_bp']), marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(summary['binning_label'])
            i += 1
    elif plot_type == 'p':
        for summary in summary_per_query:
            axs.plot(float(summary['a_rand_index_by_bp']), float(summary['percent_assigned_bps']), marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(summary['binning_label'])
            i += 1

    # turn on grid
    axs.minorticks_on()
    axs.grid(which='major', linestyle='-', linewidth='0.5')
    axs.grid(which='minor', linestyle=':', linewidth='0.5')

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}%'.format(x * 100) for x in vals])
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')
    lgd = plt.legend(plot_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)

    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def plot_avg_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'e',
                 'avg_precision_recall',
                 'Average precision per bin', #'Truncated average precision per bin $\overline{p}_{99}$',
                 'Average recall per genome') #'Average recall per genome $\overline{r}$')


def plot_weighed_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'w',
                 'avg_precision_recall_per_bp',
                 'Average precision per base pair', # 'Average precision per base pair $\overline{p}_{bp}$',
                 'Average recall per base pair') # 'Average recall per base pair $\overline{r}_{bp}$')


def plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'p',
                 'ari_vs_assigned_bps',
                 'Adjusted Rand Index',
                 'Percentage of assigned base pairs')


def main():
    parser = argparse.ArgumentParser(description="Create plots from one or more tables of results")
    parser.add_argument("files", nargs='+', help="File(s) including system path")
    parser.add_argument('-o', '--output_dir', help="Directory to save the plots in", required=True)
    args = parser.parse_args()
    load_data.make_sure_path_exists(args.output_dir)
    results = load_results(args.files)
    plot_avg_precision_recall(results, args.output_dir)
    plot_weighed_precision_recall(results, args.output_dir)
    plot_adjusted_rand_index_vs_assigned_bps(results, args.output_dir)


if __name__ == "__main__":
    main()
