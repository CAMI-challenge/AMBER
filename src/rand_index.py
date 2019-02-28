#!/usr/bin/env python

import os, sys, inspect
from collections import defaultdict
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src.utils import labels
from src import binning_classes


def choose2(n):
    return (n * (n - 1)) / 2.0


def preprocess_counts(bin_ids, query, by_bp_counts):
    if not bin_ids:
        return
    bin_id_to_mapping_id_to_length = defaultdict(lambda: defaultdict(int))
    mapping_id_to_bin_id_to_length = defaultdict(lambda: defaultdict(int))

    if isinstance(query, binning_classes.GenomeQuery):
        gs_sequence_id_to_mapping_id = query.gold_standard.sequence_id_to_bin_id
        sequence_id_to_bin_id = query.sequence_id_to_bin_id
    else:
        rank = query.get_bin_by_id(bin_ids[0]).rank
        gs_sequence_id_to_mapping_id = query.gold_standard.rank_to_sequence_id_to_bin_id[rank]
        sequence_id_to_bin_id = query.rank_to_sequence_id_to_bin_id[rank]

    for bin in query.get_bins_by_id(bin_ids):
        for sequence_id in bin.sequence_ids:
            if sequence_id in gs_sequence_id_to_mapping_id:
                mapping_id = gs_sequence_id_to_mapping_id[sequence_id]
                bin_id_to_mapping_id_to_length[bin.id][mapping_id] += binning_classes.Bin.sequence_id_to_length[sequence_id] if by_bp_counts else 1

    for sequence_id in sequence_id_to_bin_id:
        if sequence_id in gs_sequence_id_to_mapping_id:
            mapping_id = gs_sequence_id_to_mapping_id[sequence_id]
            bin_id = sequence_id_to_bin_id[sequence_id]
            mapping_id_to_bin_id_to_length[mapping_id][bin_id] += binning_classes.Bin.sequence_id_to_length[sequence_id] if by_bp_counts else 1

    return bin_id_to_mapping_id_to_length, mapping_id_to_bin_id_to_length


def combinations(bin_id_to_mapping_id_to_counts, mapping_id_to_bin_id_to_counts):
    bin_mapping_comb = 0.0
    bin_comb = 0.0
    mapping_comb = 0.0
    total_counts = 0
    for bin_id in bin_id_to_mapping_id_to_counts:
        for genome_id in bin_id_to_mapping_id_to_counts[bin_id]:
            bin_mapping_comb += choose2(bin_id_to_mapping_id_to_counts[bin_id][genome_id])
    for bin_id in bin_id_to_mapping_id_to_counts:
        bin_totals = 0
        for genome_id in bin_id_to_mapping_id_to_counts[bin_id]:
            bin_totals += bin_id_to_mapping_id_to_counts[bin_id][genome_id]
        total_counts += bin_totals
        bin_comb += choose2(bin_totals)
    for genome_id in mapping_id_to_bin_id_to_counts:
        genome_totals = 0
        for bin_id in mapping_id_to_bin_id_to_counts[genome_id]:
            genome_totals += mapping_id_to_bin_id_to_counts[genome_id][bin_id]
        mapping_comb += choose2(genome_totals)
    return bin_comb, mapping_comb, bin_mapping_comb, choose2(total_counts)


def compute_adjusted_rand_index(bin_id_to_mapping_id_to_counts, mapping_id_to_bin_id_to_counts):
    bin_comb, mapping_comb, bin_mapping_comb, num_bp_comb = combinations(bin_id_to_mapping_id_to_counts, mapping_id_to_bin_id_to_counts)
    temp = (bin_comb * mapping_comb / num_bp_comb) if num_bp_comb != 0 else .0
    ret = bin_mapping_comb - temp
    denominator = (((bin_comb + mapping_comb) / 2.0) - temp)
    return (ret / denominator) if denominator != 0 else .0


def compute_rand_index(bin_id_to_mapping_id_to_counts, mapping_id_to_bin_id_to_counts):
    bin_comb, mapping_comb, bin_mapping_comb, num_bp_comb = combinations(bin_id_to_mapping_id_to_counts, mapping_id_to_bin_id_to_counts)
    return ((num_bp_comb - bin_comb - mapping_comb + 2 * bin_mapping_comb) / num_bp_comb) if num_bp_comb != 0 else .0


def print_rand_indices(ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq, percentage_of_assigned_bps, stream=sys.stdout):
    stream.write("%s\n%s\n" % ("\t".join((labels.RI_BY_BP, labels.RI_BY_SEQ, labels.ARI_BY_BP, labels.ARI_BY_SEQ, labels.PERCENTAGE_ASSIGNED_BPS)),
                 "\t".join((format(ri_by_bp, '.3f'),
                            format(ri_by_seq, '.3f'),
                            format(ari_by_bp, '.3f'),
                            format(ari_by_seq, '.3f'),
                            format(percentage_of_assigned_bps, '.3f')))))


def compute_metrics(bin_ids, query):
    if not bin_ids:
        return np.nan, np.nan, np.nan, np.nan
    bin_id_to_mapping_id_to_total_sequences, mapping_id_to_bin_id_to_total_sequences = preprocess_counts(bin_ids, query, False)
    bin_id_to_mapping_id_to_length, mapping_id_to_bin_id_to_length = preprocess_counts(bin_ids, query, True)

    ri_by_seq = compute_rand_index(bin_id_to_mapping_id_to_total_sequences, mapping_id_to_bin_id_to_total_sequences)
    ri_by_bp = compute_rand_index(bin_id_to_mapping_id_to_length, mapping_id_to_bin_id_to_length)
    ari_by_seq = compute_adjusted_rand_index(bin_id_to_mapping_id_to_total_sequences, mapping_id_to_bin_id_to_total_sequences)
    ari_by_bp = compute_adjusted_rand_index(bin_id_to_mapping_id_to_length, mapping_id_to_bin_id_to_length)

    return ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq
