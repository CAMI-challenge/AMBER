#!/usr/bin/env python

import sys
import numpy as np
import argparse
from utils import filter_tail
from utils import load_data
from utils import argparse_parents


def get_defaults():
    return [.5, .7, .9], [.1, .05]


def calc_table(metrics, min_completeness, max_contamination):
    if not min_completeness:
        min_completeness = get_defaults()[0]
    if not max_contamination:
        max_contamination = get_defaults()[1]

    genome_recovery = np.zeros((len(max_contamination), len(min_completeness)), dtype=int)
    for metric in metrics:
        precision = float(metric['purity'])
        recall = float(metric['completeness'])
        if np.isnan(precision):
            continue
        contamination = 1 - precision
        for i in range(len(max_contamination)):
            for j in range(len(min_completeness)):
                if recall > min_completeness[j] and contamination < max_contamination[i]:
                    genome_recovery[i][j] += 1
    return genome_recovery


def print_table(genome_recovery, label, min_completeness, max_contamination, stream=sys.stdout):
    if not label:
        label = ""
    if not min_completeness:
        min_completeness = get_defaults()[0]
    if not max_contamination:
        max_contamination = get_defaults()[1]

    completeness = "".join(">{}% complete\t".format(str(int(x * 100))) for x in min_completeness).rstrip("\t")
    contamination = ["<{}% contamination".format(str(int(x * 100))) for x in max_contamination]

    stream.write("%s\t%s\n" % (label, completeness))

    for y_label, num_genomes in zip(contamination, genome_recovery):
        contamination_values = "".join("{}\t".format(str(values)) for values in num_genomes).rstrip("\t")
        stream.write("%s\t%s\n" % (y_label, contamination_values))


def main():
    parser = argparse.ArgumentParser(description="Calculate number of genome bins recovered with more than the specified thresholds of completeness and contamination. Default: >50%, >70%, >90% completeness vs. <10%, <5% contamination",
                                     parents=[argparse_parents.PARSER_MULTI])
    parser.add_argument('-x', '--min_completeness', help="Comma-separated list of min. completeness thresholds (default: 50,70,90)", required=False)
    parser.add_argument('-y', '--max_contamination', help="Comma-separated list of max. contamination thresholds (default: 10,5)", required=False)
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    metrics = load_data.load_tsv_table(sys.stdin if not sys.stdin.isatty() else args.file)

    min_completeness = None
    max_contamination = None
    if args.min_completeness:
        min_completeness = [int(x.strip())/100.0 for x in args.min_completeness.split(',')]
    if args.max_contamination:
        max_contamination = [int(x.strip())/100.0 for x in args.max_contamination.split(',')]

    if args.filter:
        metrics = filter_tail.filter_tail(metrics, args.filter)
    results = calc_table(metrics, min_completeness, max_contamination)
    print_table(results, args.label, min_completeness, max_contamination)


if __name__ == "__main__":
    main()
