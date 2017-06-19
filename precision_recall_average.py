#!/usr/bin/python

import sys
import argparse
import numpy as np


def load_data(stream):
    data = []
    for line in stream:
        line = line.strip()
        if len(line) == 0 or line.startswith("@"):
            continue
        row_data = line.split('\t')
        bin = row_data[0]
        real_size = int(float(row_data[5]))
        predicted_size = int(float(row_data[3]))

        if row_data[1] != "NA" and predicted_size > 0:
            precision = float(row_data[1])
        else:
            precision = np.nan
        data.append({'bin': bin, 'precision': precision, 'recall': row_data[2],
                     'predicted_size': predicted_size, 'real_size': real_size})
    return data


def filter_tail(data):
    # sort bins by increasing predicted size
    data = sorted(data, key=lambda x: x['predicted_size'])
    sum_of_predicted_sizes = sum([int(float(bin['predicted_size'])) for bin in data])
    cumsum_of_predicted_sizes = 0
    for bin in data:
        predicted_size = int(float(bin['predicted_size']))
        cumsum_of_predicted_sizes += predicted_size
        if float(cumsum_of_predicted_sizes) / float(sum_of_predicted_sizes) <= .01:
            bin['precision'] = np.nan
        else:
            break
    return data


def compute_precision_and_recall(data, filter_tail_bool):
    avg_precision = .0
    avg_recall = .0
    count_p = 0
    count_r = 0

    if filter_tail_bool:
        data = filter_tail(data)

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
    return avg_precision, avg_recall


def main():
    parser = argparse.ArgumentParser(description="Compute precision and recall from file or standard input")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help="File containing precision and recall for each bin")
    parser.add_argument('-f', '--filter', action="store_true", help="Filter out 1%% smallest bins - default is false")
    args = parser.parse_args()
    filter_tail_bool = args.filter
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    data = load_data(sys.stdin if not sys.stdin.isatty() else args.file)
    avgprecision, avgrecall = compute_precision_and_recall(data, filter_tail_bool)
    print "%1.3f" % avgprecision
    print "%1.3f" % avgrecall


if __name__ == "__main__":
    main()



# main("../metadata/ANI/unique_common.tsv",
#      "../binning/tables/low_unsupervised_by_genome_all.tsv",
#      "../binning/tables/medium_unsupervised_by_genome_all.tsv",
#      "../binning/tables/high_unsupervised_by_genome_all.tsv",
#      "../binning/tables/low_supervised_by_bin_all.tsv",
#      "../binning/tables/medium_supervised_by_bin_all.tsv",
#      "../binning/tables/high_supervised_by_bin_all.tsv")
