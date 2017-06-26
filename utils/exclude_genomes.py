#!/usr/bin/python

import sys
import argparse


def load_unique_common(unique_common_file_path):
    genome_to_unique_common = {}
    with open(unique_common_file_path) as read_handler:
        for line in read_handler:
            genome_to_unique_common[line.split("\t")[0]] = line.split("\t")[1].strip('\n')
    return genome_to_unique_common


def print_filtered_data(stream, unique_common_file_path, keyword):
    genome_to_unique_common = load_unique_common(unique_common_file_path)
    for line in stream:
        line = line.strip()
        if len(line) == 0 or line.startswith("@"):
            print line
            continue
        bin = line.split('\t')[0]
        if bin in genome_to_unique_common and (keyword is None or genome_to_unique_common[bin] == keyword):
            continue
        print line


def filter_data(bin_metrics, unique_common_file_path, keyword):
    genome_to_unique_common = load_unique_common(unique_common_file_path)
    filtered_bin_metrics = []
    for bin in bin_metrics:
        bin_id = bin['mapped_genome']
        if bin_id not in genome_to_unique_common or genome_to_unique_common[bin_id] != keyword:
            filtered_bin_metrics.append(bin)
    return filtered_bin_metrics


def main():
    parser = argparse.ArgumentParser(description="Exclude genome bins from table file of precision and recall or standard input")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help="File containing precision and recall for each genome")
    parser.add_argument('-r', '--genomes_file', help="File with list of genomes to be removed", required=True)
    parser.add_argument('-k', '--keyword', help="Keyword in second column of input for bins to be removed (no keyword=remove all in list)", required=False)
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    print_filtered_data(sys.stdin if not sys.stdin.isatty() else args.file,
                        args.genomes_file,
                        args.keyword)


if __name__ == "__main__":
    main()