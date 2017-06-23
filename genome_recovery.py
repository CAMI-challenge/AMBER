#!/usr/bin/python

import sys
import precision_recall_average
import numpy as np
import argparse


def calc_table(metrics):
    genome_recovery = [0, 0, 0, 0, 0, 0]  # >0.9/<0.05 | >0.9/<0.1 | >0.7/<0.05 | >0.7/<0.1 | >0.5/<0.05 | >0.5/<0.1
    for metric in metrics:
        precision = float(metric['precision'])
        recall = float(metric['recall'])
        if np.isnan(precision):
            continue
        if recall > 0.9:
            if precision > 0.95:
                genome_recovery[0] += 1
            if precision > 0.9:
                genome_recovery[1] += 1
        if recall > 0.7:
            if precision > 0.95:
                genome_recovery[2] += 1
            if precision > 0.9:
                genome_recovery[3] += 1
        if recall > 0.5:
            if precision > 0.95:
                genome_recovery[4] += 1
            if precision > 0.9:
                genome_recovery[5] += 1

    return genome_recovery


def print_table(genome_recovery, label):
    if not label:
        label = ""
    line = "%s\t>50%% complete\t>70%% complete\t>90%% complete" % label
    print(line)
    line = "<10%% contamination\t%s\t%s\t%s" % (genome_recovery[5], genome_recovery[3], genome_recovery[1])
    print(line)
    line = "<5%% contamination\t%s\t%s\t%s" % (genome_recovery[4], genome_recovery[2], genome_recovery[0])
    print(line)


def main():
    parser = argparse.ArgumentParser(description="Compute number of genomes in ranges of completeness and contamination")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help="File containing precision and recall for each genome")
    parser.add_argument('-s', '--filter', action="store_true", help="Filter out 1%% smallest bins - default is false")
    parser.add_argument('-l', '--label', help="Binning name", required=False)
    args = parser.parse_args()
    filter_tail_bool = args.filter
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    metrics = precision_recall_average.load_tsv_table(sys.stdin if not sys.stdin.isatty() else args.file)
    if filter_tail_bool:
        metrics = precision_recall_average.load_tsv_table(sys.stdin)
    metrics = precision_recall_average.filter_tail(metrics)
    print_table(calc_table(metrics), args.label)


if __name__ == "__main__":
    main()
