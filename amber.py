#!/usr/bin/env python

import argparse
import os
import errno
import matplotlib
import numpy as np
from collections import OrderedDict
from collections import defaultdict
from version import __version__
from src import genome_recovery
from src import plot_by_genome
from src import plots
from src import precision_recall_per_bin
from src import rand_index
from src import amber_html
matplotlib.use('Agg')
import pandas as pd
from src.utils import load_data
from src.utils import argparse_parents
from src.utils import labels as utils_labels
from src import binning_classes


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def create_output_directories(output_dir, queries_list):
    for query in queries_list:
        make_sure_path_exists(os.path.join(output_dir, query.binning_type, query.label))


def get_labels(labels, bin_files):
    if labels:
        labels_list = [x.strip() for x in labels.split(',')]
        if len(set(labels_list)) != len(bin_files):
            exit('Number of different labels does not match the number of binning files. Please check parameter -l, --labels.')
        return labels_list
    tool_id = []
    for bin_file in bin_files:
        tool_id.append(bin_file.split('/')[-1])
    return tool_id


def compute_metrics_over_bins(rank, gs_pd_bins_rank, pd_bins_rank, query):
    true_positive_bps_all_bins = pd_bins_rank['true_positive_bps'].sum()
    true_positive_seqs_all_bins = pd_bins_rank['true_positive_seqs'].sum()

    all_bins_length = pd_bins_rank['predicted_size'].sum()
    all_bins_num_seqs = pd_bins_rank['predicted_num_seqs'].sum()
    all_bins_false_positive_bps = all_bins_length - true_positive_bps_all_bins
    all_bins_false_positive_seqs = all_bins_num_seqs - true_positive_seqs_all_bins

    precision_bp = float(true_positive_bps_all_bins) / float(all_bins_length)
    precision_seq = float(true_positive_seqs_all_bins) / float(all_bins_num_seqs)
    misclassification_rate_bp = all_bins_false_positive_bps / float(all_bins_length)
    misclassification_rate_seq = all_bins_false_positive_seqs / float(all_bins_num_seqs)

    true_positives_recall_bp = 0
    true_positives_recall_seq = 0
    for i in gs_pd_bins_rank.index:
        mapping_id = gs_pd_bins_rank.at[i, 'mapping_id']
        bin_assigns = []
        bin_assigns_seqs = []
        for i2 in pd_bins_rank.index:
            bin_id = pd_bins_rank.at[i2, 'id']
            if bin_id:
                bin = query.get_bin_by_id(bin_id)
                if mapping_id in bin.mapping_id_to_length:
                    bin_assigns.append(bin.mapping_id_to_length[mapping_id])
                    bin_assigns_seqs.append(bin.mapping_id_to_num_seqs[mapping_id])
        if len(bin_assigns) > 0:
            true_positives_recall_bp += max(bin_assigns)
            true_positives_recall_seq += max(bin_assigns_seqs)

    gs_length = gs_pd_bins_rank['true_size'].sum()
    gs_num_seqs = gs_pd_bins_rank['true_size'].sum()
    recall_bp = float(true_positives_recall_bp) / float(gs_length)
    recall_seq = float(true_positives_recall_seq) / float(gs_num_seqs)
    accuracy_bp = float(true_positive_bps_all_bins) / float(gs_length)
    accuracy_seq = float(true_positive_seqs_all_bins) / float(gs_num_seqs)
    percentage_of_assigned_bps = float(all_bins_length) / float(gs_length)

    if isinstance(query, binning_classes.TaxonomicQuery):
        if rank in query.rank_to_overbinned_seqs:
            length_overbinned_seqs = sum([binning_classes.Bin.sequence_id_to_length[sequence_id] for sequence_id in query.rank_to_overbinned_seqs[rank]])
            percentage_of_overbinned_bps = length_overbinned_seqs / float(gs_length)
        else:
            percentage_of_overbinned_bps = .0
    else:
        percentage_of_overbinned_bps = np.nan

    return precision_bp, recall_bp, accuracy_bp, misclassification_rate_bp,\
           precision_seq, recall_seq, accuracy_seq, misclassification_rate_seq,\
           percentage_of_assigned_bps, percentage_of_overbinned_bps


def compute_percentage_of_assigned_seqs(gold_standard, query):
    percentage_of_assigned_seqs = {}
    percentage_of_overbinned_seqs = {}
    if isinstance(query, binning_classes.GenomeQuery):
        num_seqs = len(query.get_sequence_ids())
        gs_num_seqs = len(gold_standard.genome_query.get_sequence_ids())
        percentage_of_assigned_seqs['NA'] = num_seqs / gs_num_seqs
        percentage_of_overbinned_seqs['NA'] = np.nan
    else:
        num_seqs = defaultdict(int)
        gs_num_seqs = defaultdict(int)
        for bin in gold_standard.taxonomic_query.bins:
            gs_num_seqs[bin.rank] += len(bin.sequence_ids)
        for bin in query.bins:
            num_seqs[bin.rank] += len(bin.sequence_ids)
        for rank in num_seqs.keys():
            percentage_of_assigned_seqs[rank] = float(num_seqs[rank]) / float(gs_num_seqs[rank])
        for rank in query.rank_to_overbinned_seqs.keys():
            if rank in num_seqs:
                percentage_of_overbinned_seqs[rank] = float(len(query.rank_to_overbinned_seqs[rank])) / float(gs_num_seqs[rank])
    return percentage_of_assigned_seqs, percentage_of_overbinned_seqs


def evaluate_all(gold_standard,
                 queries_list,
                 min_completeness, max_contamination):
    if gold_standard.genome_query:
        gs_genome_bins_metrics = gold_standard.genome_query.get_bins_metrics(gold_standard)
        gs_pd_genome_bins = pd.DataFrame.from_dict(gs_genome_bins_metrics)
        gs_pd_genome_bins['rank'] = 'NA'
    if gold_standard.taxonomic_query:
        gs_tax_bins_metrics = gold_standard.taxonomic_query.get_bins_metrics(gold_standard)
        gs_pd_tax_bins = pd.DataFrame.from_dict(gs_tax_bins_metrics)

    pd_bins_all = pd.DataFrame()
    df_summary = pd.DataFrame()
    for query in queries_list:
        percentage_of_assigned_seqs, percentage_of_overbinned_seqs = compute_percentage_of_assigned_seqs(gold_standard, query)

        # Compute metrics per bin
        query.compute_true_positives(gold_standard)
        precision_recall_per_bin.compute_precision_recall(gold_standard, query)
        bins_metrics = query.get_bins_metrics(gold_standard)

        pd_bins = pd.DataFrame.from_dict(bins_metrics)
        pd_bins[utils_labels.TOOL] = query.label

        if isinstance(query, binning_classes.GenomeQuery):
            pd_bins['rank'] = 'NA'
            gs_pd_bins = gs_pd_genome_bins
        else:
            gs_pd_bins = gs_pd_tax_bins
        pd_bins_all = pd.concat([pd_bins_all, pd_bins], ignore_index=True, sort=False)

        # Compute metrics over bins
        for rank, pd_bins_rank in pd_bins.groupby('rank'):
            gs_pd_bins_rank = gs_pd_bins[gs_pd_bins['rank'] == rank]
            precision_bp_rows = pd_bins_rank[pd_bins_rank['purity_bp'].notnull()]['purity_bp']
            precision_seq_rows = pd_bins_rank[pd_bins_rank['purity_seq'].notnull()]['purity_seq']
            recall_bp_rows = pd_bins_rank[pd_bins_rank['true_size'] > 0]['completeness_bp']
            recall_seq_rows = pd_bins_rank[pd_bins_rank['true_size'] > 0]['completeness_seq']

            avg_precision_bp = precision_bp_rows.mean()
            sem_precision_bp = precision_bp_rows.sem()
            avg_recall_bp = recall_bp_rows.mean()
            sem_recall_bp = recall_bp_rows.sem()

            avg_precision_seq = precision_seq_rows.mean()
            sem_precision_seq = precision_seq_rows.sem()
            avg_recall_seq = recall_seq_rows.mean()
            sem_recall_seq = recall_seq_rows.sem()

            precision_bp, recall_bp, accuracy_bp, misclassification_rate_bp,\
                precision_seq, recall_seq, accuracy_seq, misclassification_rate_seq,\
                percentage_of_assigned_bps, percentage_of_overbinned_bps = compute_metrics_over_bins(
                rank, gs_pd_bins_rank, pd_bins_rank, query)

            bin_ids = pd_bins_rank['id'][pd_bins_rank['id'].notnull()].tolist()
            ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq = rand_index.compute_metrics(bin_ids, query, gold_standard)

            genome_recovery_val = genome_recovery.calc_dict(pd_bins_rank, min_completeness, max_contamination)

            df = pd.DataFrame(OrderedDict([(utils_labels.TOOL, query.label),
                                           (utils_labels.BINNING_TYPE, query.binning_type),
                                           (utils_labels.RANK, rank),
                                           (utils_labels.AVG_PRECISION_BP, [avg_precision_bp]),
                                           (utils_labels.AVG_PRECISION_BP_SEM, [sem_precision_bp]),
                                           (utils_labels.AVG_RECALL_BP, [avg_recall_bp]),
                                           (utils_labels.AVG_RECALL_BP_SEM, [sem_recall_bp]),

                                           (utils_labels.AVG_PRECISION_SEQ, [avg_precision_seq]),
                                           (utils_labels.AVG_PRECISION_SEQ_SEM, [sem_precision_seq]),
                                           (utils_labels.AVG_RECALL_SEQ, [avg_recall_seq]),
                                           (utils_labels.AVG_RECALL_SEQ_SEM, [sem_recall_seq]),

                                           (utils_labels.AVG_PRECISION_PER_BP, [precision_bp]),
                                           (utils_labels.AVG_PRECISION_PER_SEQ, [precision_seq]),

                                           (utils_labels.AVG_RECALL_PER_BP, [recall_bp]),
                                           (utils_labels.AVG_RECALL_PER_SEQ, [recall_seq]),

                                           (utils_labels.ACCURACY_PER_BP, [accuracy_bp]),
                                           (utils_labels.ACCURACY_PER_SEQ, [accuracy_seq]),

                                           (utils_labels.PERCENTAGE_ASSIGNED_BPS, [percentage_of_assigned_bps]),
                                           (utils_labels.PERCENTAGE_ASSIGNED_SEQS, [percentage_of_assigned_seqs[rank]]),
                                           (utils_labels.PERCENTAGE_ASSIGNED_BPS_UNKNOWN, [percentage_of_overbinned_bps]),
                                           (utils_labels.PERCENTAGE_ASSIGNED_SEQS_UNKNOWN, [percentage_of_overbinned_seqs[rank] if rank in percentage_of_overbinned_seqs else .0]),
                                           (utils_labels.RI_BY_BP, [ri_by_bp]),
                                           (utils_labels.RI_BY_SEQ, [ri_by_seq]),
                                           (utils_labels.ARI_BY_BP, [ari_by_bp]),
                                           (utils_labels.ARI_BY_SEQ, [ari_by_seq]),

                                           (utils_labels.MISCLASSIFICATION_PER_BP, [misclassification_rate_bp]),
                                           (utils_labels.MISCLASSIFICATION_PER_SEQ, [misclassification_rate_seq])]))

            df_genome_recovery = pd.DataFrame.from_dict(genome_recovery_val, orient='index').T
            df = df.join(df_genome_recovery)
            df_summary = pd.concat([df_summary, df], ignore_index=True)
    return df_summary, pd_bins_all


def plot_heat_maps(gold_standard, queries_list, output_dir):
    for query in queries_list:
        if isinstance(query, binning_classes.GenomeQuery):
            df_confusion = precision_recall_per_bin.transform_confusion_matrix(gold_standard, query)
            plots.plot_heatmap(df_confusion, os.path.join(output_dir, query.label))


def plot_genome_binning(gold_standard, queries_list, df_summary, pd_bins, plot_heatmaps, output_dir):
    df_summary_g = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'genome']
    if len(df_summary_g) == 0:
        return

    if plot_heatmaps:
        plot_heat_maps(gold_standard, queries_list, output_dir)

    plots.create_legend(df_summary_g, output_dir)
    plots.plot_avg_precision_recall(df_summary_g, output_dir)
    plots.plot_weighed_precision_recall(df_summary_g, output_dir)
    plots.plot_adjusted_rand_index_vs_assigned_bps(df_summary_g, output_dir)

    pd_bins_g = pd_bins[pd_bins['rank'] == 'NA']
    plots.plot_boxplot(pd_bins_g, 'purity', output_dir)
    plots.plot_boxplot(pd_bins_g, 'completeness', output_dir)

    plot_by_genome.plot_precision_recall_per_bin(pd_bins_g, output_dir)


def plot_taxonomic_binning(df_summary, output_dir):
    df_summary_t = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'taxonomic']
    for rank, pd_group in df_summary_t.groupby('rank'):
        plots.plot_avg_precision_recall(pd_group, output_dir, rank)
        plots.plot_weighed_precision_recall(pd_group, output_dir, rank)
        plots.plot_adjusted_rand_index_vs_assigned_bps(pd_group, output_dir, rank)


def main(args=None, tax_id_to_parent=None, tax_id_to_rank=None, tax_id_to_name=None):
    parser = argparse.ArgumentParser(description="AMBER: Assessment of Metagenome BinnERs",
                                     parents=[argparse_parents.PARSER_MULTI2], prog='AMBER')
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    parser.add_argument('--stdout', help="Print summary to stdout", action='store_true')
    parser.add_argument('-d', '--desc', help="Description for HTML page", required=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    group_g = parser.add_argument_group('genome binning-specific arguments')
    group_g.add_argument('-n', '--min_length', help="Minimum length of sequences", type=int, required=False)
    group_g.add_argument('-m', '--map_by_completeness', help=argparse_parents.HELP_MAP_BY_RECALL, action='store_true')
    group_g.add_argument('-x', '--min_completeness', help=argparse_parents.HELP_THRESHOLDS_COMPLETENESS, required=False)
    group_g.add_argument('-y', '--max_contamination', help=argparse_parents.HELP_THRESHOLDS_CONTAMINATION, required=False)
    group_g.add_argument('-c', '--plot_heatmaps', help="Plot heatmaps of confusion matrices (can take some minutes)", action='store_true')
    group_g.add_argument('-p', '--filter', help=argparse_parents.HELP_FILTER)
    group_g.add_argument('-r', '--remove_genomes', help=argparse_parents.HELP_GENOMES_FILE, required=False)
    group_g.add_argument('-k', '--keyword', help=argparse_parents.HELP_KEYWORD, required=False)

    group_t = parser.add_argument_group('taxonomic binning-specific arguments')
    group_t.add_argument('--ncbi_nodes_file', help="NCBI nodes file", required=False)
    group_t.add_argument('--ncbi_names_file', help="NCBI names file", required=False)

    args = parser.parse_args(args)

    min_completeness = None
    max_contamination = None
    if args.min_completeness:
        min_completeness = [int(x.strip())/100.0 for x in args.min_completeness.split(',')]
    if args.max_contamination:
        max_contamination = [int(x.strip())/100.0 for x in args.max_contamination.split(',')]

    labels = get_labels(args.labels, args.bin_files)

    if not tax_id_to_parent:
        tax_id_to_parent, tax_id_to_rank, tax_id_to_name = load_data.load_ncbi_info(args.ncbi_nodes_file,
                                                                                    args.ncbi_names_file,)
    gold_standard, queries_list = load_data.load_queries(args.gold_standard_file,
                                                         args.fasta_file,
                                                         args.bin_files,
                                                         args.map_by_completeness,
                                                         args.filter,
                                                         args.remove_genomes,
                                                         args.keyword,
                                                         tax_id_to_parent,
                                                         tax_id_to_rank,
                                                         tax_id_to_name,
                                                         args.min_length,
                                                         labels)

    output_dir = os.path.abspath(args.output_dir)
    create_output_directories(output_dir, queries_list)

    df_summary, pd_bins = evaluate_all(gold_standard,
                                       queries_list,
                                       min_completeness, max_contamination)
    df_summary.to_csv(os.path.join(output_dir, 'results.tsv'), sep='\t', index=False, float_format='%.3f')
    if args.stdout:
        print(df_summary.to_string(index=False))

    plot_genome_binning(gold_standard, queries_list, df_summary, pd_bins, args.plot_heatmaps, output_dir)
    plot_taxonomic_binning(df_summary, output_dir)
    plots.plot_taxonomic_results(df_summary, output_dir)

    pd_bins_g = pd_bins[pd_bins['rank'] == 'NA']
    for tool, pd_group in pd_bins_g.groupby(utils_labels.TOOL):
        columns = ['id', 'mapping_id', 'purity', 'completeness', 'predicted_size', 'true_positives', 'true_size']
        table = pd_group[columns].rename(columns={'id': 'bin_id', 'mapping_id': 'mapped_genome'})
        table.to_csv(os.path.join(output_dir, 'genome', tool, 'precision_recall_per_bin.tsv'), sep='\t', index=False)

    pd_bins_t = pd_bins[pd_bins['rank'] != 'NA']
    for tool, pd_group in pd_bins_t.groupby(utils_labels.TOOL):
        columns = ['id', 'rank', 'purity_bp', 'completeness_bp', 'predicted_size', 'true_positive_bps', 'true_size']
        table = pd_group[columns].rename(columns={'id': 'tax_id'})
        table.to_csv(os.path.join(output_dir, 'taxonomic', tool, 'precision_recall_per_bin.tsv'), sep='\t', index=False)

    amber_html.create_html(gold_standard, df_summary, pd_bins, args.output_dir, args.desc)


if __name__ == "__main__":
    main()
