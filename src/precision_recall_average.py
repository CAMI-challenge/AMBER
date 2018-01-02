#!/usr/bin/env python

import argparse
import math
import sys

import numpy as np
from src.utils import argparse_parents
from src.utils import filter_tail
from src.utils import load_data

from src.utils import labels


def compute_precision_and_recall(data, filter_tail_percentage):
    avg_precision = .0
    avg_recall = .0
    count_p = 0
    count_r = 0

    if filter_tail_percentage:
        filter_tail.filter_tail(data, filter_tail_percentage)

    # for each predicted bin (row of the table)
    for bin in data:
        # compute the average recall over all bins if the mapped genome size > 0
        real_size = float(bin['real_size'])
        if real_size > 0:
            recall = float(bin['completeness'])
            count_r += 1
            current_avg = avg_recall
            avg_recall = (recall - current_avg) / count_r + current_avg

        # compute the average precision over all bins
        if not np.isnan(bin['purity']):
            precision = bin['purity']
            count_p += 1
            current_avg = avg_precision
            avg_precision = (precision - current_avg) / count_p + current_avg

    sum_diffs_precision = .0
    sum_diffs_recall = .0
    for bin in data:
        real_size = float(bin['real_size'])
        if real_size > 0:
            recall = float(bin['completeness'])
            sum_diffs_recall += math.pow(recall - avg_recall, 2)
        if not np.isnan(bin['purity']):
            precision = bin['purity']
            sum_diffs_precision += math.pow(precision - avg_precision, 2)

    std_deviation_precision = math.sqrt(sum_diffs_precision / count_p)
    std_error_precision = std_deviation_precision / math.sqrt(count_p)

    std_deviation_recall = math.sqrt(sum_diffs_recall / count_r)
    std_error_recall = std_deviation_recall / math.sqrt(count_r)

    return avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall


def print_precision_recall_table_header(stream=sys.stdout):
    stream.write("%s\n" % "\t".join((labels.TOOL, labels.AVG_PRECISION, labels.STD_DEV_PRECISION,
                                     labels.SEM_PRECISION, labels.AVG_RECALL, labels.STD_DEV_RECALL, labels.
                                     SEM_RECALL)))


def print_precision_recall(label, avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall,
                           stream=sys.stdout):
    if not label:
        label = ""
    stream.write("%s\n" % "\t".join((label,
                                    format(avg_precision, '.3f'),
                                    format(std_deviation_precision, '.3f'),
                                    format(std_error_precision, '.3f'),
                                    format(avg_recall, '.3f'),
                                    format(std_deviation_recall, '.3f'),
                                    format(std_error_recall, '.3f'))))


def main():
    parser = argparse.ArgumentParser(description="Compute precision and recall, including standard deviation and standard error of the mean, from table of precision and recall per genome. The table can be provided as file or via the standard input",
                                     parents=[argparse_parents.PARSER_MULTI])
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    metrics = load_data.load_tsv_table(sys.stdin if not sys.stdin.isatty() else args.file)
    avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall =\
        compute_precision_and_recall(metrics, args.filter)
    print_precision_recall_table_header()
    print_precision_recall(args.label, avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall)

if __name__ == "__main__":
    main()
