#!/usr/bin/env python

import argparse
import sys
from collections import Counter
from collections import defaultdict

from utils import argparse_parents
from utils import load_data


def calc_precision_recall(bin_id_to_genome_id_to_total_length, bin_id_to_total_lengths, genome_id_to_total_length):
    # now calculate precision looping bins
    true_positives_p = 0
    denominator_p = sum(bin_id_to_total_lengths.values())

    for predicted_bin, genome_assigns in list(bin_id_to_genome_id_to_total_length.items()):
        # get maximal genome assignment
        if len(genome_assigns) > 0:
            maxAssign = max(genome_assigns.values())
        else:
            maxAssign = 0.0

        true_positives_p += maxAssign

    # now calculate precision as TP/(FP + TP)
    precision = float(true_positives_p) / float(denominator_p)

    # now calculate recall looping genomes
    true_positives_r = 0.0
    for genome_id in genome_id_to_total_length:
        # now loop bins
        bin_assigns = []
        for bin_id in bin_id_to_total_lengths:
            if genome_id in bin_id_to_genome_id_to_total_length[bin_id]:
                bin_assigns.append(bin_id_to_genome_id_to_total_length[bin_id][genome_id])
        if len(bin_assigns) > 0:
            maxAssign = max(bin_assigns)
        else:
            maxAssign = 0.0
        true_positives_r += maxAssign

    denominator_r = sum(genome_id_to_total_length.values())
    recall = true_positives_r / denominator_r

    return precision, recall


def compute_metrics(query, gold_standard):
    """
    @param query
    @param gold_standard
    @return: query, gold_standard
    """
    bin_id_to_total_length = {}
    bin_id_to_genome_id_to_total_length = defaultdict(Counter)
    for predicted_bin in query.bin_id_to_list_of_sequence_id:
        if predicted_bin not in bin_id_to_total_length:
            bin_id_to_total_length[predicted_bin] = 0
        for sequence_id in query.bin_id_to_list_of_sequence_id[predicted_bin]:
            bin_id_to_total_length[predicted_bin] += gold_standard.sequence_id_to_lengths[sequence_id]
            genome_id = gold_standard.sequence_id_to_genome_id[sequence_id]
            bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += gold_standard.sequence_id_to_lengths[sequence_id]
    precision, recall = calc_precision_recall(bin_id_to_genome_id_to_total_length, bin_id_to_total_length,
                                              gold_standard.genome_id_to_total_length)
    return precision, recall


def print_precision_recall_by_bpcount(precision, recall, stream=sys.stdout):
    stream.write("precision\trecall\n%1.3f\t%1.3f\n" % (precision, recall))


def main():
    parser = argparse.ArgumentParser(description="Compute precision and recall weighed by base pair counts (not averaged over genome bins) from binning file",
                                     parents=[argparse_parents.PARSER_GS])
    args = parser.parse_args()
    if not args.bin_file:
        parser.print_help()
        parser.exit(1)
    gold_standard = load_data.get_genome_mapping(args.gold_standard_file, args.fasta_file)
    query = load_data.open_query(args.bin_file)
    precision, recall = compute_metrics(query, gold_standard)
    print_precision_recall_by_bpcount(precision, recall)


if __name__ == "__main__":
    main()
