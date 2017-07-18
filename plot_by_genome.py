#!/usr/bin/env python

import sys
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
from utils import load_data
from utils import argparse_parents


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

    fig, ax1 = plt.subplots(figsize=(9, 5))

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
