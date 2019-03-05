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
from src.utils import load_ncbi_taxinfo


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

    precision_bp = (true_positive_bps_all_bins / all_bins_length) if all_bins_length > 0 else np.nan
    precision_seq = (true_positive_seqs_all_bins / all_bins_num_seqs) if all_bins_num_seqs > 0 else np.nan
    misclassification_rate_bp = (all_bins_false_positive_bps / all_bins_length) if all_bins_length > 0 else np.nan
    misclassification_rate_seq = (all_bins_false_positive_seqs / all_bins_num_seqs) if all_bins_num_seqs > 0 else np.nan

    gs_length = gs_pd_bins_rank['true_size'].sum()
    gs_num_seqs = gs_pd_bins_rank['true_num_seqs'].sum()

    if rank == 'NA':
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
        recall_bp = true_positives_recall_bp / gs_length
        recall_seq = true_positives_recall_seq / gs_num_seqs
    else:
        recall_bp = true_positive_bps_all_bins / gs_length
        recall_seq = true_positive_seqs_all_bins / gs_num_seqs

    accuracy_bp = true_positive_bps_all_bins / gs_length
    accuracy_seq = true_positive_seqs_all_bins / gs_num_seqs
    percentage_of_assigned_bps = all_bins_length / gs_length

    return precision_bp, recall_bp, accuracy_bp, misclassification_rate_bp,\
           precision_seq, recall_seq, accuracy_seq, misclassification_rate_seq,\
           percentage_of_assigned_bps


def compute_percentage_of_assigned_seqs(query):
    percentage_of_assigned_seqs = {}
    if isinstance(query, binning_classes.GenomeQuery):
        num_seqs = len(query.get_sequence_ids())
        gs_num_seqs = len(query.gold_standard.get_sequence_ids())
        percentage_of_assigned_seqs['NA'] = num_seqs / gs_num_seqs
    else:
        num_seqs = defaultdict(int)
        gs_num_seqs = defaultdict(int)
        for bin in query.gold_standard.bins:
            gs_num_seqs[bin.rank] += len(bin.sequence_ids)
        for bin in query.bins:
            num_seqs[bin.rank] += len(bin.sequence_ids)
        for rank in num_seqs.keys():
            percentage_of_assigned_seqs[rank] = num_seqs[rank] / gs_num_seqs[rank]
        for rank in load_ncbi_taxinfo.RANKS:
            if rank not in percentage_of_assigned_seqs:
                percentage_of_assigned_seqs[rank] = .0
    return percentage_of_assigned_seqs


def evaluate_all(queries_list, min_completeness, max_contamination):
    pd_bins_all = pd.DataFrame()
    df_summary = pd.DataFrame()
    for query in queries_list:
        percentage_of_assigned_seqs = compute_percentage_of_assigned_seqs(query)

        # Compute metrics per bin
        query.compute_true_positives()
        query.compute_precision_recall()
        bins_metrics = query.get_bins_metrics()

        pd_bins = pd.DataFrame.from_dict(bins_metrics)
        pd_bins[utils_labels.TOOL] = query.label

        gs_pd_bins = pd.DataFrame.from_dict(query.gold_standard.get_bins_metrics())

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
                percentage_of_assigned_bps = compute_metrics_over_bins(rank, gs_pd_bins_rank, pd_bins_rank, query)

            bin_ids = pd_bins_rank['id'][pd_bins_rank['id'].notnull()].tolist()
            ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq = rand_index.compute_metrics(bin_ids, query)

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


def plot_heat_maps(queries_list, output_dir):
    for query in queries_list:
        if isinstance(query, binning_classes.GenomeQuery):
            df_confusion = precision_recall_per_bin.transform_confusion_matrix(query)
            plots.plot_heatmap(df_confusion, os.path.join(output_dir, query.label))


def plot_genome_binning(queries_list, df_summary, pd_bins, plot_heatmaps, output_dir):
    df_summary_g = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'genome']
    if len(df_summary_g) == 0:
        return

    if plot_heatmaps:
        plot_heat_maps(queries_list, output_dir)

    plots.create_legend(df_summary_g, output_dir)
    plots.plot_avg_precision_recall(df_summary_g, output_dir)
    plots.plot_weighed_precision_recall(df_summary_g, output_dir)
    plots.plot_adjusted_rand_index_vs_assigned_bps(df_summary_g, output_dir)

    pd_bins_g = pd_bins[pd_bins['rank'] == 'NA']
    plots.plot_boxplot(pd_bins_g, 'purity_bp', output_dir)
    plots.plot_boxplot(pd_bins_g, 'completeness_bp', output_dir)

    plot_by_genome.plot_precision_recall_per_bin(pd_bins_g, output_dir)


def plot_taxonomic_binning(df_summary, output_dir):
    df_summary_t = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'taxonomic']
    for rank, pd_group in df_summary_t.groupby('rank'):
        plots.plot_avg_precision_recall(pd_group, output_dir, rank)
        plots.plot_weighed_precision_recall(pd_group, output_dir, rank)
        plots.plot_adjusted_rand_index_vs_assigned_bps(pd_group, output_dir, rank)


def main(args=None):
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

    options = binning_classes.Options(filter_tail_percentage=args.filter,
                                      filter_genomes_file=args.remove_genomes,
                                      filter_keyword=args.keyword,
                                      map_by_completeness=args.map_by_completeness,
                                      min_length=args.min_length)

    load_data.load_ncbi_info(args.ncbi_nodes_file, args.ncbi_names_file)

    queries_list = load_data.load_queries(args.gold_standard_file,
                                          args.fasta_file,
                                          args.bin_files,
                                          options,
                                          labels)

    output_dir = os.path.abspath(args.output_dir)
    create_output_directories(output_dir, queries_list)

    df_summary, pd_bins = evaluate_all(queries_list,
                                       min_completeness, max_contamination)
    df_summary.to_csv(os.path.join(output_dir, 'results.tsv'), sep='\t', index=False, float_format='%.3f')
    if args.stdout:
        print(df_summary.to_string(index=False))

    plot_genome_binning(queries_list, df_summary, pd_bins, args.plot_heatmaps, output_dir)
    plot_taxonomic_binning(df_summary, output_dir)
    plots.plot_taxonomic_results(df_summary, output_dir)

    pd_bins_g = pd_bins[pd_bins['rank'] == 'NA']
    for tool, pd_group in pd_bins_g.groupby(utils_labels.TOOL):
        columns = ['id', 'mapping_id', 'purity_bp', 'completeness_bp', 'predicted_size', 'true_positive_bps', 'true_size']
        table = pd_group[columns].rename(columns={'id': 'bin_id', 'mapping_id': 'mapped_genome'})
        table.to_csv(os.path.join(output_dir, 'genome', tool, 'precision_recall_per_bin.tsv'), sep='\t', index=False)

    pd_bins_t = pd_bins[pd_bins['rank'] != 'NA']
    for tool, pd_group in pd_bins_t.groupby(utils_labels.TOOL):
        columns = ['id', 'rank', 'purity_bp', 'completeness_bp', 'predicted_size', 'true_positive_bps', 'true_size']
        table = pd_group[columns].rename(columns={'id': 'tax_id'})
        table.to_csv(os.path.join(output_dir, 'taxonomic', tool, 'precision_recall_per_bin.tsv'), sep='\t', index=False)

    # TODO: use num_genomes of query gs
    num_genomes = len(queries_list[0].gold_standard.get_bins_metrics())
    amber_html.create_html(num_genomes, df_summary, pd_bins, args.output_dir, args.desc)


if __name__ == "__main__":
    main()
