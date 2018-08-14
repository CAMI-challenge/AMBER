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
import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import numpy as np
import pandas as pd
from src.utils import load_data


def g_compute_confusion_matrix_all(gold_standard, queries_list, sequence_id_to_length):
    for query in queries_list:
        g_compute_confusion_matrix(gold_standard, query, sequence_id_to_length)


def g_compute_confusion_matrix(gold_standard, query, sequence_id_to_length):
    query.bin_id_to_genome_id_to_length = {}
    query.bin_id_to_mapped_genome = {}
    for bin_id in query.bin_id_to_list_of_sequence_ids:
        query.bin_id_to_genome_id_to_length[bin_id] = {}
        for sequence_id in query.bin_id_to_list_of_sequence_ids[bin_id]:
            genome_id = gold_standard.sequence_id_to_bin_id[sequence_id]
            if genome_id not in query.bin_id_to_genome_id_to_length[bin_id]:
                query.bin_id_to_genome_id_to_length[bin_id][genome_id] = 0
            query.bin_id_to_genome_id_to_length[bin_id][genome_id] += sequence_id_to_length[sequence_id]


def map_genomes_by_recall_all(gold_standard, queries_list):
    for query in queries_list:
        map_genomes_by_recall(gold_standard, query)


def map_genomes_by_recall(gold_standard, query):
    """
    This script maps a predicted bin to the genome that maximizes recall
    @return:
    """
    query.bin_id_to_mapped_genome = {}
    query.bin_id_to_true_positives = {}
    query.mapped_genomes = set()
    for bin_id in query.bin_id_to_list_of_sequence_ids:
        max_genome_percentage = .0
        best_genome_id = ""
        for genome_id in query.bin_id_to_genome_id_to_length[bin_id]:
            genome_percentage = query.bin_id_to_genome_id_to_length[bin_id][genome_id] / gold_standard.bin_id_to_length[genome_id]
            if max_genome_percentage < genome_percentage:
                max_genome_percentage = genome_percentage
                best_genome_id = genome_id
            elif max_genome_percentage == genome_percentage and gold_standard.bin_id_to_length[genome_id] > gold_standard.bin_id_to_length[best_genome_id]:
                best_genome_id = genome_id
        query.bin_id_to_true_positives[bin_id] = gold_standard.bin_id_to_length[best_genome_id]
        query.mapped_genomes.add(best_genome_id)
        query.bin_id_to_mapped_genome[bin_id] = best_genome_id


def map_genomes_by_precision_all(queries_list):
    for query in queries_list:
        map_genomes_by_precision(query)


def map_genomes_by_precision(query):
    """
    This script maps a predicted bin to the genome that maximizes precision
    @return:
    """
    query.bin_id_to_mapped_genome = {}
    query.bin_id_to_true_positives = {}
    query.mapped = set()
    for bin_id in query.bin_id_to_list_of_sequence_ids:
        max_length = 0
        best_genome_id = ""
        for genome_id in query.bin_id_to_genome_id_to_length[bin_id]:
            if max_length < query.bin_id_to_genome_id_to_length[bin_id][genome_id]:
                max_length = query.bin_id_to_genome_id_to_length[bin_id][genome_id]
                best_genome_id = genome_id
        query.bin_id_to_true_positives[bin_id] = query.bin_id_to_genome_id_to_length[bin_id][best_genome_id]
        query.mapped.add(best_genome_id)
        query.bin_id_to_mapped_genome[bin_id] = best_genome_id


def transform_confusion_matrix_all(gold_standard, queries_list):
    for query in queries_list:
        transform_confusion_matrix(gold_standard, query)


def transform_confusion_matrix(gold_standard,
                             query):
    bin_id_to_mapped_genome = query.bin_id_to_mapped_genome
    bin_id_to_genome_id_to_length = query.bin_id_to_genome_id_to_length

    df_confusion = pd.DataFrame(bin_id_to_genome_id_to_length).T

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
    bin_id_to_mapped_genome_by_length = collections.OrderedDict(sorted(bin_id_to_mapped_genome.items(), key=lambda t: bin_id_to_genome_id_to_length[t[0]][t[1]], reverse=True))

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
                             query):
    # bin_metrics = []
    bin_id_to_precision = {}
    bin_id_to_recall = {}
    for bin in query.bins:
        bin_id_to_precision[bin.id] = float(bin.true_positives) / float(bin.length)

        if isinstance(bin, load_data.GenomeBin):
            gold_standard_query = gold_standard.genome_query
        else:
            gold_standard_query = gold_standard.taxonomic_query

        if bin.mapping_id in gold_standard_query.get_bin_ids():
            bin_id_to_recall[bin.id] = float(bin.true_positives) / float(gold_standard_query.get_bin_by_id(bin.mapping_id).length)
        else:
            bin_id_to_recall[bin.id] = .0

    return bin_id_to_precision, bin_id_to_recall


    #     bin_metrics.append({'mapped_genome': mapping_id,
    #                         'purity': precision,
    #                         'completeness': recall,
    #                         'predicted_size': query.bin_id_to_total_length[bin_id],
    #                         'correctly_predicted': query.bin_id_to_genome_id_to_length[bin_id][mapping_id],
    #                         'real_size': gold_standard.bin_id_to_total_length[mapping_id]})
    # for genome_id in gold_standard.bin_id_to_list_of_sequence_ids:
    #     if genome_id not in query.mapped:
    #         bin_metrics.append({'mapped_genome': genome_id,
    #                             'purity': np.nan,
    #                             'completeness': .0,
    #                             'predicted_size': 0,
    #                             'correctly_predicted': 0,
    #                             'real_size': gold_standard.bin_id_to_total_length[genome_id]})
    # # sort bins by completeness
    # return sorted(bin_metrics, key=lambda t: t['completeness'], reverse=True)


def compute_metrics(gold_standard, query):
    bin_metrics = compute_precision_recall(gold_standard,
                                           query)
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
