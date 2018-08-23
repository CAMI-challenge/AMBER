#!/usr/bin/env python

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src import plots
from src.utils import labels


def plot_by_genome(data, out_file=None, sort_by='completeness'):
    not_sort_by = list(set(['purity','completeness']) - set([sort_by]))[0]  # get the metric not sorted by
    data = sorted(data, key=lambda x: x[sort_by])
    genomes = []
    precision = []
    recall = []
    for genome in data:
        genomes.append(genome['mapped_genome'])
        precision.append(genome['purity'])
        recall.append(genome['completeness'])
    sort = {'purity': precision, 'completeness': recall}

    fig, ax1 = plt.subplots(figsize=(len(genomes) * 0.15, 5))

    ax1.plot(np.arange(len(genomes)), sort[sort_by], color='black')
    plt.xticks(np.arange(len(genomes)), genomes, rotation='vertical', fontsize="smaller")
    ax1.plot(np.arange(len(genomes)), sort[not_sort_by], '.', color='red')

    # transform y labels to percentages
    vals = ax1.get_yticks()
    ax1.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    plt.legend((sort_by.title(), not_sort_by.title()))
    plt.grid(True)
    plt.tight_layout()
    if out_file is None:
        plt.show()
    else:
        plt.savefig(os.path.normpath(out_file + '.png'), dpi=100, format='png', bbox_inches='tight')
        plt.savefig(os.path.normpath(out_file + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    plt.close(fig)


def plot_precision_recall_per_bin(pd_bins, output_dir):
    colors_list = plots.create_colors_list()
    df_groups = pd_bins.groupby(labels.TOOL)
    if len(df_groups) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axis to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    # for query_metrics in bin_metrics_per_query:
    for i, (tool, pd_summary) in enumerate(df_groups):
        precision = pd_summary['purity'].tolist()
        recall = pd_summary['completeness'].tolist()
        axs.scatter(precision, recall, marker='o', color=colors_list[i], s=[8] * len(precision))

    # turn on grid
    axs.minorticks_on()
    axs.grid(which='major', linestyle='-', linewidth='0.5')
    axs.grid(which='minor', linestyle=':', linewidth='0.5')

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}%'.format(x * 100) for x in vals])
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    plt.xlabel('Purity per bin')
    plt.ylabel('Completeness per genome')
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'genome', 'purity_completeness_per_bin.eps'), dpi=100, format='eps', bbox_inches='tight')

    lgd = plt.legend(list(df_groups.groups.keys()), bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)
    for handle in lgd.legendHandles:
        handle.set_sizes([100.0])

    fig.savefig(os.path.join(output_dir, 'genome', 'purity_completeness_per_bin.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'genome', 'purity_completeness_per_bin.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)
