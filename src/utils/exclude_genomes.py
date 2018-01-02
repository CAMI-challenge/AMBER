#!/usr/bin/env python

import sys
import os
import argparse

try:
    import argparse_parents
    import load_data
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import argparse_parents
        import load_data
    finally:
        sys.path.remove(os.path.dirname(__file__))


def print_filtered_data(stream, unique_common_file_path, keyword):
    genome_to_unique_common = load_data.load_unique_common(unique_common_file_path)
    for line in stream:
        line = line.strip()
        if len(line) == 0 or line.startswith("@"):
            print(line)
            continue
        bin = line.split('\t')[0]
        if bin in genome_to_unique_common and (keyword is None or genome_to_unique_common[bin] == keyword):
            continue
        print(line)


def filter_genomes(gold_standard, unique_common_file_path, keyword):
    if not unique_common_file_path:
        return gold_standard.genome_id_to_list_of_contigs.keys()
    genome_to_unique_common = load_data.load_unique_common(unique_common_file_path)
    genome_ids = []
    for genome_id in gold_standard.genome_id_to_list_of_contigs:
        if genome_id not in genome_to_unique_common or (keyword and genome_to_unique_common[genome_id] != keyword):
            genome_ids.append(genome_id)
    return genome_ids


def filter_data(bin_metrics, unique_common_file_path, keyword):
    genome_to_unique_common = load_data.load_unique_common(unique_common_file_path)
    filtered_bin_metrics = []
    for bin in bin_metrics:
        bin_id = bin['mapped_genome']
        if bin_id not in genome_to_unique_common or (keyword and genome_to_unique_common[bin_id] != keyword):
            filtered_bin_metrics.append(bin)
    if len(filtered_bin_metrics) == 0:
        sys.exit('All bins have been filtered out due to option --remove_genomes.')
    return filtered_bin_metrics


def main():
    parser = argparse.ArgumentParser(description="Exclude genome bins from table of precision and recall per genome. The table can be provided as file or via the standard input")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help=argparse_parents.HELP_FILE)
    parser.add_argument('-r', '--remove_genomes', help=argparse_parents.HELP_GENOMES_FILE, required=True)
    parser.add_argument('-k', '--keyword', help=argparse_parents.HELP_KEYWORD, required=False)
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    print_filtered_data(sys.stdin if not sys.stdin.isatty() else args.file,
                        args.remove_genomes,
                        args.keyword)


if __name__ == "__main__":
    main()