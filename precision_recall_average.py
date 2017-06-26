#!/usr/bin/python

import sys
import argparse
import numpy as np
import math


def load_tsv_table(stream):
    data = []
    next(stream)
    for line in stream:
        line = line.strip()
        if len(line) == 0 or line.startswith("@"):
            continue
        row_data = line.split('\t')

        mapped_genome = row_data[0]
        real_size = int(float(row_data[5]))
        predicted_size = int(float(row_data[3]))
        correctly_predicted = int(float(row_data[4]))

        if row_data[1] != "NA" and predicted_size > 0:
            precision = float(row_data[1])
        else:
            precision = np.nan
        data.append({'mapped_genome': mapped_genome, 'precision': precision, 'recall': row_data[2],
                     'predicted_size': predicted_size, 'correctly_predicted': correctly_predicted, 'real_size': real_size})
    return data


def filter_tail(data, percentage):
    value = float(percentage) / 100.0
    # sort bins by increasing predicted size
    data = sorted(data, key=lambda x: x['predicted_size'])
    sum_of_predicted_sizes = sum([int(float(bin['predicted_size'])) for bin in data])
    cumsum_of_predicted_sizes = 0
    for bin in data:
        predicted_size = int(float(bin['predicted_size']))
        cumsum_of_predicted_sizes += predicted_size
        if float(cumsum_of_predicted_sizes) / float(sum_of_predicted_sizes) <= value:
            bin['precision'] = np.nan
        else:
            break
    return data


def compute_precision_and_recall(data, filter_tail_percentage):
    avg_precision = .0
    avg_recall = .0
    count_p = 0
    count_r = 0

    if filter_tail_percentage:
        data = filter_tail(data, filter_tail_percentage)

    # for each predicted bin (row of the table)
    for bin in data:
        # compute the average recall over all bins if the mapped genome size > 0
        real_size = float(bin['real_size'])
        if real_size > 0:
            recall = float(bin['recall'])
            count_r += 1
            current_avg = avg_recall
            avg_recall = (recall - current_avg) / count_r + current_avg

        # compute the average precision over all bins
        if not np.isnan(bin['precision']):
            precision = bin['precision']
            count_p += 1
            current_avg = avg_precision
            avg_precision = (precision - current_avg) / count_p + current_avg

    sum_diffs_precision = .0
    sum_diffs_recall = .0
    for bin in data:
        real_size = float(bin['real_size'])
        if real_size > 0:
            recall = float(bin['recall'])
            sum_diffs_recall += math.pow(recall - avg_recall, 2)
        if not np.isnan(bin['precision']):
            precision = bin['precision']
            sum_diffs_precision += math.pow(precision - avg_precision, 2)

    std_deviation_precision = math.sqrt(sum_diffs_precision / count_p)
    std_error_precision = std_deviation_precision / math.sqrt(count_p)

    std_deviation_recall = math.sqrt(sum_diffs_recall / count_r)
    std_error_recall = std_deviation_recall / math.sqrt(count_r)

    return avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall


def print_precision_recall_table_header():
    print "tool\tprecision\tstd_dev_precision\tsem_precision\trecall\tstd_dev_recall\tsem_recall"


def print_precision_recall(label, avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall):
    if not label:
        label = ""
    print "%s\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f" % (label,
                                                            avg_precision,
                                                            std_deviation_precision,
                                                            std_error_precision,
                                                            avg_recall,
                                                            std_deviation_recall,
                                                            std_error_recall)


def main():
    parser = argparse.ArgumentParser(description="Compute precision and recall from table file of precision and recall or standard input")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help="File containing precision and recall for each genome")
    parser.add_argument('-p', '--filter', help="Filter out [FILTER]%% smallest bins - default is 0")
    parser.add_argument('-l', '--label', help="Binning name", required=False)
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    metrics = load_tsv_table(sys.stdin if not sys.stdin.isatty() else args.file)
    avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall =\
        compute_precision_and_recall(metrics, args.filter)
    print_precision_recall_table_header()
    print_precision_recall(args.label, avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall)

if __name__ == "__main__":
    main()
