#!/usr/bin/env python

import sys
import argparse
from utils import load_data
from utils import argparse_parents
from utils import labels


def choose2(n):
    return (n * (n - 1)) / 2.0


def compute_percentage_of_assigned_bps(query, gold_standard):
    assigned_bps = 0
    for bin_id in query.bin_id_to_list_of_sequence_id:
        for sequence_id in query.bin_id_to_list_of_sequence_id[bin_id]:
            assigned_bps += gold_standard.sequence_id_to_lengths[sequence_id]
    gold_size_in_bps = sum(gold_standard.genome_id_to_total_length.values())
    return float(assigned_bps) / float(gold_size_in_bps)


def preprocess_by_bp_counts(query, gold_standard):
    bin_id_to_genome_id_to_length = {}
    for bin_id in query.bin_id_to_list_of_sequence_id:
        bin_id_to_genome_id_to_length[bin_id] = {}
        for sequence_id in query.bin_id_to_list_of_sequence_id[bin_id]:
            genome_id = gold_standard.sequence_id_to_genome_id[sequence_id]
            if genome_id not in bin_id_to_genome_id_to_length[bin_id]:
                bin_id_to_genome_id_to_length[bin_id][genome_id] = 0
            bin_id_to_genome_id_to_length[bin_id][genome_id] += gold_standard.sequence_id_to_lengths[sequence_id]
    genome_id_to_bin_id_to_length = {}
    for sequence_id in query.sequence_id_to_bin_id:
        genome_id = gold_standard.sequence_id_to_genome_id[sequence_id]
        if genome_id not in genome_id_to_bin_id_to_length:
            genome_id_to_bin_id_to_length[genome_id] = {}
        bin_id = query.sequence_id_to_bin_id[sequence_id]
        if bin_id not in genome_id_to_bin_id_to_length[genome_id]:
            genome_id_to_bin_id_to_length[genome_id][bin_id] = 0
        genome_id_to_bin_id_to_length[genome_id][bin_id] += gold_standard.sequence_id_to_lengths[sequence_id]
    return bin_id_to_genome_id_to_length, genome_id_to_bin_id_to_length


def preprocess_by_sequence_counts(query, gold_standard):
    # metrics = precision_recall_average.filter_tail(metrics, 1)
    bin_id_to_genome_id_to_total_sequences = {}
    for bin_id in query.bin_id_to_list_of_sequence_id:
        bin_id_to_genome_id_to_total_sequences[bin_id] = {}
        for sequence_id in query.bin_id_to_list_of_sequence_id[bin_id]:
            genome_id = gold_standard.sequence_id_to_genome_id[sequence_id]
            # metric_entry = filter(lambda x: x['mapped_genome'] == genome_id, metrics)
            # if len(metric_entry) > 0 and np.isnan(metric_entry[0]['precision']):
            #     continue
            if genome_id not in bin_id_to_genome_id_to_total_sequences[bin_id]:
                bin_id_to_genome_id_to_total_sequences[bin_id][genome_id] = 0
            bin_id_to_genome_id_to_total_sequences[bin_id][genome_id] += 1
    genome_id_to_bin_id_to_total_sequences = {}
    for sequence_id in query.sequence_id_to_bin_id:
        genome_id = gold_standard.sequence_id_to_genome_id[sequence_id]
        # metric_entry = filter(lambda x: x['mapped_genome'] == genome_id, metrics)
        # if len(metric_entry) > 0 and np.isnan(metric_entry[0]['precision']):
        #     continue
        if genome_id not in genome_id_to_bin_id_to_total_sequences:
            genome_id_to_bin_id_to_total_sequences[genome_id] = {}
        bin_id = query.sequence_id_to_bin_id[sequence_id]
        if bin_id not in genome_id_to_bin_id_to_total_sequences[genome_id]:
            genome_id_to_bin_id_to_total_sequences[genome_id][bin_id] = 0
        genome_id_to_bin_id_to_total_sequences[genome_id][bin_id] += 1
    return bin_id_to_genome_id_to_total_sequences, genome_id_to_bin_id_to_total_sequences


def combinations(bin_id_to_genome_id_to_counts, genome_id_to_bin_id_to_counts):
    bin_genome_comb = 0.0
    bin_comb = 0.0
    genome_comb = 0.0
    num_bp = 0
    for bin_id in bin_id_to_genome_id_to_counts:
        for genome_id in bin_id_to_genome_id_to_counts[bin_id]:
            bin_genome_comb += choose2(bin_id_to_genome_id_to_counts[bin_id][genome_id])
    for bin_id in bin_id_to_genome_id_to_counts:
        bin_totals = 0
        for genome_id in bin_id_to_genome_id_to_counts[bin_id]:
            bin_totals += bin_id_to_genome_id_to_counts[bin_id][genome_id]
        num_bp += bin_totals
        bin_comb += choose2(bin_totals)
    for genome_id in genome_id_to_bin_id_to_counts:
        genome_totals = 0
        for bin_id in genome_id_to_bin_id_to_counts[genome_id]:
            genome_totals += genome_id_to_bin_id_to_counts[genome_id][bin_id]
        genome_comb += choose2(genome_totals)
    return bin_comb, genome_comb, bin_genome_comb, choose2(num_bp)


def compute_adjusted_rand_index(bin_id_to_genome_id_to_counts, genome_id_to_bin_id_to_counts):
    bin_comb, genome_comb, bin_genome_comb, num_bp_comb = combinations(bin_id_to_genome_id_to_counts, genome_id_to_bin_id_to_counts)
    temp = bin_comb * genome_comb / num_bp_comb
    ret = bin_genome_comb - temp
    return ret / (((bin_comb + genome_comb) / 2.0) - temp)


def compute_rand_index(bin_id_to_genome_id_to_counts, genome_id_to_bin_id_to_counts):
    bin_comb, genome_comb, bin_genome_comb, num_bp_comb = combinations(bin_id_to_genome_id_to_counts, genome_id_to_bin_id_to_counts)
    return (num_bp_comb - bin_comb - genome_comb + 2 * bin_genome_comb) / num_bp_comb


def print_rand_indices(ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq, percentage_of_assigned_bps, stream=sys.stdout):
    stream.write("%s\n%s\n" % ("\t".join((labels.RI_BY_BP, labels.RI_BY_SEQ, labels.ARI_BY_BP, labels.ARI_BY_SEQ, labels.PERCENTAGE_ASSIGNED_BPS)),
                 "\t".join((format(ri_by_bp, '.3f'),
                            format(ri_by_seq, '.3f'),
                            format(ari_by_bp, '.3f'),
                            format(ari_by_seq, '.3f'),
                            format(percentage_of_assigned_bps, '.3f')))))


def compute_metrics(query, gold_standard):
    # metrics = precision_recall_average.load_tsv_table(sys.stdin)
    bin_id_to_genome_id_to_total_sequences, genome_id_to_bin_id_to_total_sequences = preprocess_by_sequence_counts(query, gold_standard)
    bin_id_to_genome_id_to_length, genome_id_to_bin_id_to_length = preprocess_by_bp_counts(query, gold_standard)

    ri_by_seq = compute_rand_index(bin_id_to_genome_id_to_total_sequences, genome_id_to_bin_id_to_total_sequences)
    ri_by_bp = compute_rand_index(bin_id_to_genome_id_to_length, genome_id_to_bin_id_to_length)
    ari_by_seq = compute_adjusted_rand_index(bin_id_to_genome_id_to_total_sequences, genome_id_to_bin_id_to_total_sequences)
    ari_by_bp = compute_adjusted_rand_index(bin_id_to_genome_id_to_length, genome_id_to_bin_id_to_length)

    percentage_of_assigned_bps = compute_percentage_of_assigned_bps(query, gold_standard)

    return ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq, percentage_of_assigned_bps


def main():
    parser = argparse.ArgumentParser(description="Compute (adjusted) rand index from binning file, unweighed and weighed by base pairs, and percentage of binned base pairs",
                                     parents=[argparse_parents.PARSER_GS])
    args = parser.parse_args()
    if not args.query_file:
        parser.print_help()
        parser.exit(1)
    gold_standard = load_data.get_genome_mapping(args.gold_standard_file, args.fasta_file)
    query = load_data.open_query(args.query_file)
    ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq, percentage_of_assigned_bps = compute_metrics(query, gold_standard)
    print_rand_indices(ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq, percentage_of_assigned_bps)


if __name__ == "__main__":
    main()
