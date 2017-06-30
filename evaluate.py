#!/usr/bin/env python

import argparse
import os
import sys
import errno
import precision_recall_per_genome
import precision_recall_average
import precision_recall_by_bpcount
import ari
import genome_recovery
from utils import exclude_genomes
from utils import load_data
from utils import argparse_parents
from utils import labels


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def evaluate_all(gold_standard_file, fasta_file, query_files, labels, filter_tail_percentage, genomes_file, keyword, output_dir):
    gold_standard = load_data.get_genome_mapping(gold_standard_file, fasta_file)
    labels_iterator = iter(labels)
    summary_per_query = []
    for query_file in query_files:
        tool_id = query_file.split('/')[-1]
        binning_label = next(labels_iterator) if len(labels) > 0 else tool_id
        path = os.path.join(output_dir, tool_id)
        make_sure_path_exists(path)

        query = load_data.open_query(query_file)

        # PRECISION RECALL PER GENOME
        bin_metrics = precision_recall_per_genome.compute_metrics(query, gold_standard)
        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)
        f = open(path + "/precision_recall.tsv", 'w')
        precision_recall_per_genome.print_metrics(bin_metrics, f)
        f.close()

        # PRECISION RECALL AVG
        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_percentage)
        f = open(path + "/precision_recall_avg.tsv", 'w')
        precision_recall_average.print_precision_recall_table_header(f)
        precision_recall_average.print_precision_recall(binning_label,
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        std_error_precision,
                                                        std_error_recall,
                                                        f)
        f.close()

        # PRECISION RECALL BY BP COUNTS
        precision, recall = precision_recall_by_bpcount.compute_metrics(query, gold_standard)
        f = open(path + "/precision_recall_by_bpcount.tsv", 'w')
        precision_recall_by_bpcount.print_precision_recall_by_bpcount(precision, recall, f)
        f.close()

        # ARI
        ari_by_bp, ari_unweighed = ari.compute_metrics(query, gold_standard)
        f = open(path + "/ari.tsv", 'w')
        ari.print_ari(ari_by_bp, ari_unweighed, f)
        f.close()

        # GENOME RECOVERY
        genome_recovery_val = genome_recovery.calc_table(bin_metrics)

        summary_per_query.append((binning_label,
                                  format(avg_precision, '.3f'),
                                  format(std_deviation_precision, '.3f'),
                                  format(std_error_precision, '.3f'),
                                  format(avg_recall, '.3f'),
                                  format(std_deviation_recall, '.3f'),
                                  format(std_error_recall, '.3f'),
                                  format(precision, '.3f'),
                                  format(recall, '.3f'),
                                  format(ari_by_bp, '.3f'),
                                  format(ari_unweighed, '.3f'),
                                  str(genome_recovery_val[5]),
                                  str(genome_recovery_val[3]),
                                  str(genome_recovery_val[1]),
                                  str(genome_recovery_val[4]),
                                  str(genome_recovery_val[2]),
                                  str(genome_recovery_val[0])))
    return summary_per_query


def print_summary(summary_per_query, stream=sys.stdout):
    stream.write("%s\n" % "\t".join((labels.TOOL, labels.AVG_PRECISION, labels.STD_DEV_PRECISION,
                                     labels.SEM_PRECISION, labels.AVG_RECALL, labels.STD_DEV_RECALL, labels.
                                     SEM_RECALL, labels.PRECISION, labels.RECALL, labels.ARI_BY_BP,
                                     labels.ARI_UNWEIGHED,
                                     ">0.5compl<0.1cont",
                                     ">0.7compl<0.1cont",
                                     ">0.9compl<0.1cont",
                                     ">0.5compl<0.05cont",
                                     ">0.7compl<0.05cont",
                                     ">0.9compl<0.05cont")))
    for summary in summary_per_query:
        stream.write("%s\n" % "\t".join(summary))


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics for binning files",
                                     parents=[argparse_parents.PARSER_MULTI2])
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    args = parser.parse_args()
    binning_labels = []
    if args.labels:
        binning_labels = [x.strip() for x in args.labels.split(',')]
        if len(binning_labels) != len(args.query_files):
            parser.error('number of labels does not match the number of query files')
    summary_per_query = evaluate_all(args.gold_standard_file,
                                     args.fasta_file,
                                     args.query_files,
                                     binning_labels,
                                     args.filter,
                                     args.genomes_file,
                                     args.keyword,
                                     args.output_dir)
    print_summary(summary_per_query)


if __name__ == "__main__":
    main()
