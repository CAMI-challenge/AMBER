#!/usr/bin/env python

import sys
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
from utils import load_data
from utils import argparse_parents


def create_colors_list():
    colors_list = []
    for color in plt.cm.Set1(np.linspace(0, 1, 9)):
        colors_list.append(tuple(color))
    for color in plt.cm.Set2(np.linspace(0, 1, 8)):
        colors_list.append(tuple(color))
    for color in plt.cm.Set3(np.linspace(0, 1, 12)):
        colors_list.append(tuple(color))
    return colors_list


def plot_by_genome(data, out_file=None, sort_by='recall'):
    not_sort_by = list(set(['precision','recall']) - set([sort_by]))[0]  # get the metric not sorted by
    data = sorted(data, key=lambda x: x[sort_by])
    genomes = []
    precision = []
    recall = []
    for genome in data:
        genomes.append(genome['mapped_genome'])
        precision.append(genome['precision'])
        recall.append(genome['recall'])
    sort = {'precision': precision, 'recall': recall}

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


def plot_by_genome2(summary_per_query, output_dir):
    colors_list = create_colors_list()
    if len(summary_per_query) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axis to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    i = 0
    plot_labels = []
    for summary in summary_per_query:
        precision = []
        recall = []
        for metrics in summary[1]:
            precision.append(metrics['precision'])
            recall.append(metrics['recall'])
        axs.scatter(precision, recall, marker='o', color=colors_list[i])
        plot_labels.append(summary[0]['binning_label'])
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

    lgd = plt.legend(plot_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)
    plt.xlabel('Precision')
    plt.ylabel('Recall')
    plt.tight_layout()
    fig.savefig(os.path.normpath(output_dir + '/precision_recall_per_bin.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/precision_recall_per_bin.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot precision and recall per genome. Genomes can be sorted by recall (default) or precision")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help=argparse_parents.HELP_FILE)
    parser.add_argument('-s','--sort_by', help='Sort by either precision or recall (default: recall)', choices=set(['precision','recall']))
    parser.add_argument('-o','--out_file', help='Path to store image (default: only show image)')
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    metrics = load_data.load_tsv_table(sys.stdin if not sys.stdin.isatty() else args.file)
    if args.sort_by is not None:
        plot_by_genome(metrics, args.out_file, args.sort_by)
    else:
        plot_by_genome(metrics, args.out_file)
        

if __name__ == "__main__":
    main()
