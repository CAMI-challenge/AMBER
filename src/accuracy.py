#!/usr/bin/env python

import argparse
from collections import Counter
from collections import defaultdict

from src.utils import argparse_parents
from src.utils import load_data


def calc_accuracy(bin_id_to_genome_id_to_total_length, gold_standard):
    # query_sequence_ids = set(query.sequence_id_to_bin_id.keys())
    # gs_sequence_ids = set(gold_standard.sequence_id_to_lengths.keys())
    # unassigned = 0
    # for unassigned_seq_id in gs_sequence_ids - query_sequence_ids:
    #     unassigned += gold_standard.sequence_id_to_lengths[unassigned_seq_id]
    # denominator_p = sum(bin_id_to_total_lengths.values()) + unassigned

    # summing up assigned and unassigned bps is equivalent to summing up the length of all genomes
    denominator_p = 0
    for genome_id in gold_standard.genome_id_to_total_length.keys():
        denominator_p += gold_standard.genome_id_to_total_length[genome_id]

    true_positives_p = 0
    for predicted_bin, genome_assigns in list(bin_id_to_genome_id_to_total_length.items()):
        # get maximal genome assignment
        if len(genome_assigns) > 0:
            max_assign = max(genome_assigns.values())
        else:
            max_assign = 0.0
        true_positives_p += max_assign

    # calculate accuracy as TP/(FP + TP + U)
    accuracy = float(true_positives_p) / float(denominator_p)
    return accuracy


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
    accuracy = calc_accuracy(bin_id_to_genome_id_to_total_length, gold_standard)
    return accuracy


def main():
    parser = argparse.ArgumentParser(description="Compute accuracy from binning file",
                                     parents=[argparse_parents.PARSER_GS])
    args = parser.parse_args()
    if not args.bin_file:
        parser.print_help()
        parser.exit(1)
    gold_standard = load_data.get_genome_mapping(args.gold_standard_file, args.fasta_file)
    query = load_data.open_query(args.bin_file)
    accuracy = compute_metrics(query, gold_standard)
    print(accuracy)

if __name__ == "__main__":
    main()
