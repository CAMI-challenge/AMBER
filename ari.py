#!/usr/bin/env python

import sys
import argparse
from utils import load_data


def choose2(n):
    return (n * (n - 1)) / 2.0


def compute_ari(bin_id_to_genome_id_to_counts, genome_id_to_bin_id_to_counts):
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

    temp = bin_comb * genome_comb / choose2(num_bp)
    ret = bin_genome_comb - temp
    return ret / (((bin_comb + genome_comb) / 2.0) - temp)


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


def print_ari(ari_by_bp, ari_unweighed, stream=sys.stdout):
    stream.write("ari_weighed_by_bp\tari_unweighed\n%1.3f\t%1.3f\n" % (ari_by_bp, ari_unweighed))


def compute_metrics(file_path_mapping, file_path_query, file_fasta):
    # metrics = precision_recall_average.load_tsv_table(sys.stdin)
    gold_standard = load_data.get_genome_mapping(file_path_mapping, file_fasta)
    query = load_data.open_query(file_path_query)

    bin_id_to_genome_id_to_total_sequences, genome_id_to_bin_id_to_total_sequences = preprocess_by_sequence_counts(query, gold_standard)
    ari_unweighted = compute_ari(bin_id_to_genome_id_to_total_sequences, genome_id_to_bin_id_to_total_sequences)

    bin_id_to_genome_id_to_length, genome_id_to_bin_id_to_length = preprocess_by_bp_counts(query, gold_standard)
    ari_weighted = compute_ari(bin_id_to_genome_id_to_length, genome_id_to_bin_id_to_length)
    return ari_unweighted, ari_weighted


def main():
    parser = argparse.ArgumentParser(description="Compute adjusted rand index from binning file")
    parser.add_argument("-g", "--gold_standard_file", help="gold standard - ground truth - file", required=True)
    parser.add_argument("query_file", help="Query file")
    parser.add_argument("-f", "--fasta_file",
                        help="FASTA or FASTQ file w/ sequences of gold standard - required if gold standard file misses column _LENGTH")
    args = parser.parse_args()
    if not args.query_file:
        parser.print_help()
        parser.exit(1)
    ari_unweighed, ari_by_bp = compute_metrics(file_path_mapping=args.gold_standard_file,
                                               file_path_query=args.query_file,
                                               file_fasta=args.fasta_file)
    print_ari(ari_by_bp, ari_unweighed)


if __name__ == "__main__":
    main()
