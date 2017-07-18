#!/usr/bin/env python

import argparse
import precision_recall_per_genome
import precision_recall_average
from utils import exclude_genomes
from utils import load_data
from utils import argparse_parents


def evaluate_all(gold_standard, queries, labels, filter_tail_percentage, genomes_file, keyword):
    precision_recall_average.print_precision_recall_table_header()
    labels_iterator = iter(labels)
    for query in queries:

        bin_metrics = precision_recall_per_genome.compute_metrics(query, gold_standard)

        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)

        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_percentage)

        precision_recall_average.print_precision_recall(next(labels_iterator) if len(labels) > 0 else query.path.split('/')[-1],
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        std_error_precision,
                                                        std_error_recall)


def main():
    parser = argparse.ArgumentParser(description="Compute precision and recall, including standard deviation and standard error of the mean, for binning files",
                                     parents=[argparse_parents.PARSER_MULTI2])
    args = parser.parse_args()
    labels = []
    if args.labels:
        labels = [x.strip() for x in args.labels.split(',')]
        if len(labels) != len(args.bin_files):
            parser.error('number of labels does not match the number of query files')
    gold_standard = load_data.get_genome_mapping(args.gold_standard_file, args.fasta_file)
    queries = load_data.open_queries(args.bin_files)
    evaluate_all(gold_standard,
                 queries,
                 labels,
                 args.filter,
                 args.genomes_file,
                 args.keyword)


if __name__ == "__main__":
    main()
