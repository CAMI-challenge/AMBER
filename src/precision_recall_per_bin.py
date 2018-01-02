#!/usr/bin/env python

# Script for computing precision and recall. It takes as input:
# - a gold standard file in bioboxes format
# (https://github.com/bioboxes/rfc/blob/4bb19a633a6a969c2332f1f298852114c5f89b1b/data-format/binning.mkd)
# with optional column _LENGTH
# - a (compressed) fasta or fastq file, required if _LENGTH is not present in the gold standard file
# - the bins to be evaluated in the same format as above
# It writes to standard output a table containing precision and recall for each bin.

import argparse
import collections
import sys

import numpy as np
import pandas as pd

from utils import argparse_parents
from utils import load_data


def map_genomes_by_recall(gold_standard, bin_id_to_list_of_sequence_id):
    """
    This script maps a predicted bin to the genome that maximizes recall

    @attention: In case of reads, read ids might not be paired read id and cause error: ReadID/1 ReadID/2

    @param sequence_id_to_genome_id:
    @param anonymous_sequence_id_to_lengths:
    @param bin_id_to_list_of_sequence_id
    @return:
    """
    bin_id_to_genome_id_to_total_length = {}
    mapped_genomes = set()
    bin_id_to_mapped_genome = {}
    for predicted_bin in bin_id_to_list_of_sequence_id:
        bin_id_to_genome_id_to_total_length[predicted_bin] = {}
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            genome_id = gold_standard.sequence_id_to_genome_id[sequence_id]
            if genome_id not in bin_id_to_genome_id_to_total_length[predicted_bin]:
                bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] = 0
            bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += gold_standard.sequence_id_to_lengths[sequence_id]
        max_genome_percentage = .0
        best_genome_id = ""
        for genome_id in bin_id_to_genome_id_to_total_length[predicted_bin]:
            genome_percentage = bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] / gold_standard.genome_id_to_total_length[genome_id]
            if max_genome_percentage < genome_percentage:
                max_genome_percentage = genome_percentage
                best_genome_id = genome_id
            elif max_genome_percentage == genome_percentage and gold_standard.genome_id_to_total_length[genome_id] > gold_standard.genome_id_to_total_length[best_genome_id]:
                best_genome_id = genome_id
        mapped_genomes.add(best_genome_id)
        bin_id_to_mapped_genome[predicted_bin] = best_genome_id
    return bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes


def map_genomes(gold_standard, bin_id_to_list_of_sequence_id):
    """
    This script maps a predicted bin to the genome that maximizes precision

    @attention: In case of reads, read ids might not be paired read id and cause error: ReadID/1 ReadID/2

    @param sequence_id_to_genome_id:
    @param anonymous_sequence_id_to_lengths:
    @param bin_id_to_list_of_sequence_id
    @return:
    """
    bin_id_to_genome_id_to_total_length = {}
    mapped = set()
    bin_id_to_mapped_genome = {}
    for predicted_bin in bin_id_to_list_of_sequence_id:
        bin_id_to_genome_id_to_total_length[predicted_bin] = {}
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            genome_id = gold_standard.sequence_id_to_genome_id[sequence_id]
            if genome_id not in bin_id_to_genome_id_to_total_length[predicted_bin]:
                bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] = 0
            bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += gold_standard.sequence_id_to_lengths[sequence_id]
        max_length = 0
        best_genome_id = ""
        for genome_id in bin_id_to_genome_id_to_total_length[predicted_bin]:
            if max_length < bin_id_to_genome_id_to_total_length[predicted_bin][genome_id]:
                max_length = bin_id_to_genome_id_to_total_length[predicted_bin][genome_id]
                best_genome_id = genome_id
        mapped.add(best_genome_id)
        bin_id_to_mapped_genome[predicted_bin] = best_genome_id
    return bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped


def compute_confusion_matrix(bin_id_to_mapped_genome,
                             bin_id_to_genome_id_to_total_length,
                             gold_standard,
                             query):
    df_confusion = pd.DataFrame(bin_id_to_genome_id_to_total_length).T

    query_sequence_ids = set(query.sequence_id_to_bin_id.keys())
    gs_sequence_ids = set(gold_standard.sequence_id_to_lengths.keys())
    genome_id_to_unassigned_bps = collections.Counter()
    for unassigned_seq_id in gs_sequence_ids - query_sequence_ids:
        genome_id = gold_standard.sequence_id_to_genome_id[unassigned_seq_id]
        genome_id_to_unassigned_bps[genome_id] += gold_standard.sequence_id_to_lengths[unassigned_seq_id]

    df_unassigned = pd.DataFrame.from_dict(genome_id_to_unassigned_bps, orient='index').rename(columns={0: 'unassigned'}).T
    table = df_confusion.append(df_unassigned)
    table.fillna(0, inplace=True)
    # use log scale
    #table = table.applymap(np.log).fillna(0)

    # sort bins by the number of true positives (length of mapped genome within the bin)
    bin_id_to_mapped_genome_by_length = collections.OrderedDict(sorted(bin_id_to_mapped_genome.items(), key=lambda t: bin_id_to_genome_id_to_total_length[t[0]][t[1]], reverse=True))

    # sort genomes
    genome_order = []
    for bin_id in bin_id_to_mapped_genome_by_length:
        mapped_genome = bin_id_to_mapped_genome_by_length[bin_id]
        if mapped_genome not in genome_order:
            genome_order.append(mapped_genome)
    genome_order += list(set(table.columns.values.tolist()) - set(genome_order))
    for genome_id in genome_id_to_unassigned_bps.keys():
        if genome_id not in genome_order:
            genome_order.append(genome_id)

    table = table.loc[list(bin_id_to_mapped_genome_by_length.keys()) + ['unassigned'], genome_order]

    for genome_id in gold_standard.genome_id_to_total_length.keys():
        if genome_id not in table.columns.values.tolist():
            table[genome_id] = 0

    return table


def compute_precision_recall(gold_standard,
                             bin_id_to_list_of_sequence_id,
                             bin_id_to_mapped_genome,
                             bin_id_to_genome_id_to_total_length,
                             mapped):
    bin_id_to_total_lengths = {}
    for predicted_bin in bin_id_to_list_of_sequence_id:
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            if predicted_bin not in bin_id_to_total_lengths:
                bin_id_to_total_lengths[predicted_bin] = 0
            bin_id_to_total_lengths[predicted_bin] += gold_standard.sequence_id_to_lengths[sequence_id]

    bin_metrics = []
    for predicted_bin in bin_id_to_list_of_sequence_id:
        best_genome_id = bin_id_to_mapped_genome[predicted_bin]
        # length of genome in bin divided by bin size
        precision = float(bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id]) / float(bin_id_to_total_lengths[predicted_bin])
        recall = float(bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id]) / float(gold_standard.genome_id_to_total_length[best_genome_id])
        bin_metrics.append({'mapped_genome': best_genome_id,
                            'purity': precision,
                            'completeness': recall,
                            'predicted_size': bin_id_to_total_lengths[predicted_bin],
                            'correctly_predicted': bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id],
                            'real_size': gold_standard.genome_id_to_total_length[best_genome_id]})
    for genome_id in gold_standard.genome_id_to_list_of_contigs:
        if genome_id not in mapped:
            bin_metrics.append({'mapped_genome': genome_id,
                                'purity': np.nan,
                                'completeness': .0,
                                'predicted_size': 0,
                                'correctly_predicted': 0,
                                'real_size': gold_standard.genome_id_to_total_length[genome_id]})
    # sort bins by completeness
    return sorted(bin_metrics, key=lambda t: t['completeness'], reverse=True)


def compute_metrics(query, gold_standard, bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes):
    bin_metrics = compute_precision_recall(gold_standard,
                                           query.bin_id_to_list_of_sequence_id,
                                           bin_id_to_mapped_genome,
                                           bin_id_to_genome_id_to_total_length,
                                           mapped_genomes)
    return bin_metrics


def print_metrics(bin_metrics, stream=sys.stdout):
    stream.write("genome\tpurity\tcompleteness\tpredicted_size\tcorrectly_predicted\treal_size\n")
    for bin in bin_metrics:
        stream.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
            bin['mapped_genome'],
            bin['purity'] if not np.isnan(bin['purity']) else 'NA',
            bin['completeness'],
            bin['predicted_size'],
            bin['correctly_predicted'],
            bin['real_size']))


def main():
    parser = argparse.ArgumentParser(description="Compute table of purity and completeness per genome bin",
                                     parents=[argparse_parents.PARSER_GS])
    parser.add_argument('-m', '--map_by_completeness', help=argparse_parents.HELP_MAP_BY_RECALL, action='store_true')
    args = parser.parse_args()
    if not args.gold_standard_file or not args.bin_file:
        parser.print_help()
        parser.exit(1)
    gold_standard = load_data.get_genome_mapping(args.gold_standard_file, args.fasta_file)
    query = load_data.open_query(args.bin_file)

    if args.map_by_completeness:
            bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes = map_genomes_by_recall(gold_standard, query.bin_id_to_list_of_sequence_id)
    else:
        bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes = map_genomes(gold_standard, query.bin_id_to_list_of_sequence_id)

    bin_metrics = compute_metrics(query, gold_standard, bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes)
    print_metrics(bin_metrics)


if __name__ == "__main__":
    main()
