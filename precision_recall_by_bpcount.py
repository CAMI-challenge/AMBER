#!/usr/bin/python

import argparse
from utils import load_data
from collections import Counter
from collections import defaultdict


def calc_precision_recall(bin_id_to_genome_id_to_total_length, bin_id_to_total_lengths, genome_id_to_total_length):
    # now calculate precision looping bins
    true_positives_p = 0
    denominator_p = sum(bin_id_to_total_lengths.values())
    false_positives_p = 0.0

    for predicted_bin, genome_assigns in bin_id_to_genome_id_to_total_length.iteritems():
        # get maximal genome assignment
        if len(genome_assigns) > 0:
            maxAssign = max(genome_assigns.values())
        else:
            maxAssign = 0.0

        true_positives_p += maxAssign

    # now calculate precision as TP/FP + TP
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


def compute_metrics(file_path_mapping, file_path_query, file_fasta):
    """
    This script calculates precision and recall for binned contigs

    @param file_path_mapping:
    @param file_path_query:
    @return:
    """
    gold_standard = load_data.get_genome_mapping(file_path_mapping, file_fasta)
    bin_id_to_list_of_sequence_id = {}
    bin_id_to_total_length = {}
    with open(file_path_query) as read_handler:
        for sequence_id, predicted_bin, length in load_data.read_binning_file(read_handler):
            if predicted_bin not in bin_id_to_total_length:
                bin_id_to_list_of_sequence_id[predicted_bin] = []
                bin_id_to_total_length[predicted_bin] = 0
            bin_id_to_list_of_sequence_id[predicted_bin].append(sequence_id)
            bin_id_to_total_length[predicted_bin] += gold_standard.contig_id_to_lengths[sequence_id]

    bin_id_to_genome_id_to_total_length = defaultdict(Counter)
    for predicted_bin in bin_id_to_list_of_sequence_id:
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            genome_id = gold_standard.contig_id_to_genome_id[sequence_id]
            bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += gold_standard.contig_id_to_lengths[sequence_id]

    precision, recall = calc_precision_recall(bin_id_to_genome_id_to_total_length, bin_id_to_total_length,
                                              gold_standard.genome_id_to_total_length)
    print "precision\trecall\n%1.3f\t%1.3f" % (precision, recall)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gold_standard_file", help="gold standard - ground truth - file", required=True)
    parser.add_argument("-q", "--query_file", help="query file", required=True)
    parser.add_argument("-f", "--fasta_file",
                        help="FASTA or FASTQ file w/ sequences of gold standard - required if gold standard file misses column _LENGTH")
    args = parser.parse_args()
    if not args.gold_standard_file or not args.query_file:
        parser.print_help()
        parser.exit(1)
    compute_metrics(file_path_mapping=args.gold_standard_file,
                    file_path_query=args.query_file,
                    file_fasta=args.fasta_file)


if __name__ == "__main__":
    main()
