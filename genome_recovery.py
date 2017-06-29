#!/usr/bin/env python

import sys
import numpy as np
import argparse
from utils import filter_tail
from utils import load_data
from utils import argparse_parents


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


def print_table(genome_recovery, label, stream=sys.stdout):
    if not label:
        label = ""
    stream.write("%s\t>50%% complete\t>70%% complete\t>90%% complete\n" % label)
    stream.write("<10%% contamination\t%s\t%s\t%s\n" % (genome_recovery[5], genome_recovery[3], genome_recovery[1]))
    stream.write("<5%% contamination\t%s\t%s\t%s\n" % (genome_recovery[4], genome_recovery[2], genome_recovery[0]))


def main():
    parser = argparse.ArgumentParser(description="Compute number of genomes in ranges of completeness and contamination",
                                     parents=[argparse_parents.PARSER_MULTI])
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    metrics = load_data.load_tsv_table(sys.stdin if not sys.stdin.isatty() else args.file)
    if args.filter:
        metrics = filter_tail.filter_tail(metrics, args.filter)
    print_table(calc_table(metrics), args.label)


if __name__ == "__main__":
    main()
