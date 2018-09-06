#!/usr/bin/env python

import argparse
import os
import errno
import matplotlib
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


def compute_metrics_per_bp(gs_pd_bins_rank, pd_bins_rank, query):
    true_positives_all_bins = pd_bins_rank['true_positives'].sum()
    all_bins_length = pd_bins_rank['predicted_size'].sum()
    all_bins_false_positives = all_bins_length - true_positives_all_bins
    precision_by_bp = float(true_positives_all_bins) / float(all_bins_length)
    misclassification_rate = all_bins_false_positives / float(all_bins_length)

    true_positives_recall = 0
    for gs_index, gs_row in gs_pd_bins_rank.iterrows():
        bin_assigns = []
        for index, row in pd_bins_rank.iterrows():
            if row['id']:
                bin = query.get_bin_by_id(row['id'])
                if gs_row['mapping_id'] in bin.mapping_id_to_length:
                    bin_assigns.append(bin.mapping_id_to_length[gs_row['mapping_id']])
        if len(bin_assigns) > 0:
            true_positives_recall += max(bin_assigns)

    gs_length = gs_pd_bins_rank['real_size'].sum()
    recall_by_bp = float(true_positives_recall) / float(gs_length)
    accuracy = float(true_positives_all_bins) / float(gs_length)
    percentage_of_assigned_bps = float(all_bins_length) / float(gs_length)

    return precision_by_bp, recall_by_bp, accuracy, percentage_of_assigned_bps, misclassification_rate


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
            precision_rows = pd_bins_rank[pd_bins_rank['purity'].notnull()]['purity']
            recall_rows = pd_bins_rank[pd_bins_rank['real_size'] > 0]['completeness']

            avg_precision = precision_rows.mean()
            sem_precision = precision_rows.sem()
            std_precision = precision_rows.std()
            avg_recall = recall_rows.mean()
            sem_recall = recall_rows.sem()
            std_recall = recall_rows.std()

            precision_by_bp, recall_by_bp, accuracy, percentage_of_assigned_bps, misclassification_rate = compute_metrics_per_bp(
                gs_pd_bins_rank, pd_bins_rank, query)

            bin_ids = pd_bins_rank['id'][pd_bins_rank['id'].notnull()].tolist()
            ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq = rand_index.compute_metrics(bin_ids, query, gold_standard)

            genome_recovery_val = genome_recovery.calc_dict(pd_bins_rank, min_completeness, max_contamination)

            df = pd.DataFrame({utils_labels.TOOL: query.label,
                               utils_labels.BINNING_TYPE: query.binning_type,
                               utils_labels.RANK: rank,
                               utils_labels.AVG_PRECISION: [avg_precision],
                               utils_labels.AVG_PRECISION_STD: [std_precision],
                               utils_labels.AVG_PRECISION_SEM: [sem_precision],
                               utils_labels.AVG_RECALL: [avg_recall],
                               utils_labels.AVG_RECALL_STD: [std_recall],
                               utils_labels.AVG_RECALL_SEM: [sem_recall],
                               utils_labels.AVG_PRECISION_PER_BP: [precision_by_bp],
                               utils_labels.AVG_RECALL_PER_BP: [recall_by_bp],
                               utils_labels.ACCURACY: [accuracy],
                               utils_labels.PERCENTAGE_ASSIGNED_BPS: [percentage_of_assigned_bps],
                               utils_labels.RI_BY_BP: [ri_by_bp],
                               utils_labels.RI_BY_SEQ: [ri_by_seq],
                               utils_labels.ARI_BY_BP: [ari_by_bp],
                               utils_labels.ARI_BY_SEQ: [ari_by_seq],
                               utils_labels.MISCLASSIFICATION: [misclassification_rate]})
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


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics and figures for one or more binning files; output summary to screen and results per binning file to chosen directory",
                                     parents=[argparse_parents.PARSER_MULTI2], prog='AMBER')
    parser.add_argument('-n', '--min_length', help="Minimum length of sequences", type=int, required=False)
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)

    parser.add_argument('--ncbi_nodes_file', help="NCBI nodes file", required=False)
    parser.add_argument('-m', '--map_by_completeness', help=argparse_parents.HELP_MAP_BY_RECALL, action='store_true')
    parser.add_argument('-x', '--min_completeness', help=argparse_parents.HELP_THRESHOLDS_COMPLETENESS, required=False)
    parser.add_argument('-y', '--max_contamination', help=argparse_parents.HELP_THRESHOLDS_CONTAMINATION, required=False)
    parser.add_argument('-c', '--plot_heatmaps', help="Plot heatmaps of confusion matrices (can take some minutes)", action='store_true')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()

    min_completeness = None
    max_contamination = None
    if args.min_completeness:
        min_completeness = [int(x.strip())/100.0 for x in args.min_completeness.split(',')]
    if args.max_contamination:
        max_contamination = [int(x.strip())/100.0 for x in args.max_contamination.split(',')]

    labels = get_labels(args.labels, args.bin_files)

    gold_standard, queries_list = load_data.load_queries(args.gold_standard_file,
                                                         args.fasta_file,
                                                         args.bin_files,
                                                         args.map_by_completeness,
                                                         args.filter,
                                                         args.remove_genomes,
                                                         args.keyword,
                                                         args.ncbi_nodes_file,
                                                         args.min_length,
                                                         labels)

    output_dir = os.path.abspath(args.output_dir)
    create_output_directories(output_dir, queries_list)

    df_summary, pd_bins = evaluate_all(gold_standard,
                                       queries_list,
                                       min_completeness, max_contamination)
    df_summary.to_csv(os.path.join(output_dir, 'results.tsv'), sep='\t', index=False, float_format='%.3f')

    plot_genome_binning(gold_standard, queries_list, df_summary, pd_bins, args.plot_heatmaps, output_dir)
    plots.plot_taxonomic_results(df_summary, output_dir)

    pd_bins_g = pd_bins[pd_bins['rank'] == 'NA']
    for tool, pd_group in pd_bins_g.groupby(utils_labels.TOOL):
        columns = ['id', 'mapping_id', 'purity', 'completeness', 'predicted_size', 'true_positives', 'real_size']
        table = pd_group[columns].rename(columns={'id': 'bin_id', 'mapping_id': 'mapped_genome'})
        table.to_csv(os.path.join(output_dir, 'genome', tool, 'precision_recall_per_bin.tsv'), sep='\t', index=False)

    pd_bins_t = pd_bins[pd_bins['rank'] != 'NA']
    for tool, pd_group in pd_bins_t.groupby(utils_labels.TOOL):
        columns = ['id', 'rank', 'purity', 'completeness', 'predicted_size', 'true_positives', 'real_size']
        table = pd_group[columns].rename(columns={'id': 'tax_id'})
        table.to_csv(os.path.join(output_dir, 'taxonomic', tool, 'precision_recall_per_bin.tsv'), sep='\t', index=False)

    amber_html.create_html(gold_standard, df_summary, pd_bins, args.output_dir, "my desc")


if __name__ == "__main__":
    main()
