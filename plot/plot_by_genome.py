#!/usr/bin/env python

# so we can use the utils stuff in this folder
import sys
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
import argparse
import numpy as np
from utils import load_data
import matplotlib.pyplot as plt

def plot_by_genome(by_genome_file, out_file = None, sort_by='recall'):
    not_sort_by = list(set(['precision','recall']) - set([sort_by]))[0] # get the metric not sorted by
    data = load_data.load_tsv_table(by_genome_file)
    data = sorted(data,key=lambda x:x[sort_by])
    genomes = []
    precision = []
    recall = []
    for genome in data:
        genomes.append(genome['mapped_genome'])
        precision.append(genome['precision'])
        recall.append(genome['recall'])
    sort = {'genomes':genomes,'precision':precision,'recall':recall}
    fig, ax1 = plt.subplots()
    ax1.plot(np.arange(len(sort['genomes'])), sort[sort_by], color='black')
    plt.xticks(np.arange(len(sort['genomes'])), sort['genomes'], rotation='vertical')
    ax1.plot(np.arange(len(sort['genomes'])), sort[not_sort_by], '.', color='red')
    plt.title("%s by genome" % sort_by)
    plt.legend((sort_by, not_sort_by))
    plt.ylabel("Percentage")
    plt.tight_layout()
    if out_file is None:
        plt.show()
    else:
        plt.savefig(out_file)

def main():
    parser = argparse.ArgumentParser(description="Plot sorted precision and recall per genome")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help="File containing precision and recall for each genome")
    parser.add_argument('-s','--sort_by', help='sort by either precision or recall',choices=set(['precision','recall']))
    parser.add_argument('-o','--out_file', help='Path to store image (default: show)')
    args = parser.parse_args()
    if not args.file:
        parser.print_help()
        parser.exit(1)
    if args.sort_by is not None:
        plot_by_genome(args.file, args.out_file, args.sort_by)
    else:
        plot_by_genome(args.file, args.out_file)
        

if __name__ == "__main__":
    main()
