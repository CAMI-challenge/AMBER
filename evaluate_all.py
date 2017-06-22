#!/usr/bin/python

import argparse
import precision_recall_per_genome
import precision_recall_average
from utils import exclude_genomes
from utils import load_data


def evaluate_all(gold_standard, query_files, filter_tail_bool, genomes_file, keyword):
    precision_recall_average.print_precision_recall_table_header()
    for query_file in query_files:
        bin_metrics = precision_recall_per_genome.compute_metrics(query_file, gold_standard)

        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)

        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_bool)

        precision_recall_average.print_precision_recall(query_file.split('/')[-1],
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        std_error_precision,
                                                        std_error_recall)


def main():
    parser = argparse.ArgumentParser(description="Compute metrics for multiple binning files at once")
    parser.add_argument("query_files", nargs='+', help="Query files")
    parser.add_argument("-g", "--gold_standard_file", help="gold standard - ground truth - file", required=True)
    parser.add_argument("-f", "--fasta_file",
                        help="FASTA or FASTQ file w/ sequences of gold standard - required if gold standard file misses column _LENGTH")
    parser.add_argument('-s', '--filter', action="store_true", help="Filter out 1%% smallest bins - default is false")
    parser.add_argument('-r', '--genomes_file', help="File with list of genomes to be removed", required=False)
    parser.add_argument('-k', '--keyword', help="Keyword in second column of input for bins to be removed (no keyword=remove all in list)", required=False)
    args = parser.parse_args()
    gold_standard = load_data.get_genome_mapping(args.gold_standard_file, args.fasta_file)
    evaluate_all(gold_standard,
                 args.query_files,
                 args.filter,
                 args.genomes_file,
                 args.keyword)


if __name__ == "__main__":
    main()
