#!/usr/bin/env python

import argparse
import os
import errno
import precision_recall_per_genome
import precision_recall_average
import ari
from utils import exclude_genomes
from utils import load_data
from utils import argparse_parents


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def evaluate_all(gold_standard_file, fasta_file, query_files, labels, filter_tail_percentage, genomes_file, keyword, output_dir):
    gold_standard = load_data.get_genome_mapping(gold_standard_file, fasta_file)
    labels_iterator = iter(labels)
    for query_file in query_files:
        tools_id = query_file.split('/')[-1]
        path = os.path.join(output_dir, tools_id)
        make_sure_path_exists(path)

        query = load_data.open_query(query_file)

        # PRECISION RECALL PER GENOME
        f_precision_recall = open(path + "/precision_recall.tsv", 'w')
        bin_metrics = precision_recall_per_genome.compute_metrics(query, gold_standard)
        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)
        precision_recall_per_genome.print_metrics(bin_metrics, f_precision_recall)
        f_precision_recall.close()

        # PRECISION RECALL AVG
        f_precision_recall_avg = open(path + "/precision_recall_avg.tsv", 'w')
        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_percentage)
        precision_recall_average.print_precision_recall_table_header(f_precision_recall_avg)
        precision_recall_average.print_precision_recall(next(labels_iterator) if len(labels) > 0 else tools_id,
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        std_error_precision,
                                                        std_error_recall,
                                                        f_precision_recall_avg)
        f_precision_recall_avg.close()

        # ARI
        f_ari = open(path + "/ari.tsv", 'w')
        ari_by_bp, ari_unweighed = ari.compute_metrics(query, gold_standard)
        ari.print_ari(ari_by_bp, ari_unweighed, f_ari)
        f_ari.close()




def main():
    parser = argparse.ArgumentParser(description="Compute all metrics for binning files",
                                     parents=[argparse_parents.PARSER_MULTI2])
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    args = parser.parse_args()
    labels = []
    if args.labels:
        labels = [x.strip() for x in args.labels.split(',')]
        if len(labels) != len(args.query_files):
            parser.error('number of labels does not match the number of query files')
    evaluate_all(args.gold_standard_file,
                 args.fasta_file,
                 args.query_files,
                 labels,
                 args.filter,
                 args.genomes_file,
                 args.keyword,
                 args.output_dir)


if __name__ == "__main__":
    main()
