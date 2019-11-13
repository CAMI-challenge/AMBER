#!/usr/bin/env python

import argparse
import os
import sys
import errno
import matplotlib
import numpy as np
import logging
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


def get_logger(output_dir, silent):
    make_sure_path_exists(output_dir)
    logger = logging.getLogger('amber')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    logging_fh = logging.FileHandler(os.path.join(output_dir, 'log.txt'))
    logging_fh.setFormatter(formatter)
    logger.addHandler(logging_fh)

    if not silent:
        logging_stdout = logging.StreamHandler(sys.stdout)
        logging_stdout.setFormatter(formatter)
        logger.addHandler(logging_stdout)
    return logger


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def create_output_directories(output_dir, sample_id_to_queries_list):
    logging.getLogger('amber').info('Creating output directories...')
    for sample_id in sample_id_to_queries_list:
        for query in sample_id_to_queries_list[sample_id]:
            make_sure_path_exists(os.path.join(output_dir, query.binning_type, query.label))
    logging.getLogger('amber').info('done')


def get_labels(labels, bin_files):
    if labels:
        labels_list = [x.strip() for x in labels.split(',')]
        if len(set(labels_list)) != len(bin_files):
            logging.getLogger('amber').critical('Number of different labels does not match the number of binning files. Please check parameter -l, --labels.')
            exit(1)
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
            if gs_num_seqs[rank] == 0:
                percentage_of_assigned_seqs[rank] = 0
            else:
                percentage_of_assigned_seqs[rank] = num_seqs[rank] / gs_num_seqs[rank]
        for rank in load_ncbi_taxinfo.RANKS:
            if rank not in percentage_of_assigned_seqs:
                percentage_of_assigned_seqs[rank] = .0
    return percentage_of_assigned_seqs


def evaluate_samples_queries(sample_id_to_queries_list, min_completeness, max_contamination):
    pd_bins_all = pd.DataFrame()
    df_summary_all = pd.DataFrame()
    for sample_id in sample_id_to_queries_list:
        df_summary, pd_bins = evaluate_all(sample_id_to_queries_list[sample_id], sample_id, min_completeness, max_contamination)
        pd_bins_all = pd.concat([pd_bins_all, pd_bins], ignore_index=True)
        df_summary_all = pd.concat([df_summary_all, df_summary], ignore_index=True)
    return df_summary_all, pd_bins_all


def evaluate_all(queries_list, sample_id, min_completeness, max_contamination):
    pd_bins_all = pd.DataFrame()
    df_summary = pd.DataFrame()
    for query in queries_list:
        logging.getLogger('amber').info('Computing metrics for ' + query.label + ' - ' + query.binning_type + ' binning, ' + sample_id + '...')
        percentage_of_assigned_seqs = compute_percentage_of_assigned_seqs(query)

        # Compute metrics per bin
        query.compute_true_positives()
        query.compute_precision_recall()
        bins_metrics = query.get_bins_metrics()

        pd_bins = pd.DataFrame.from_dict(bins_metrics)
        pd_bins[utils_labels.TOOL] = query.label
        gs_pd_bins = pd.DataFrame.from_dict(query.gold_standard.get_bins_metrics())
        pd_bins_all = pd.concat([pd_bins_all, pd_bins], ignore_index=True, sort=False)

        unifrac_bp, unifrac_seq = query.compute_unifrac()

        # Compute metrics over bins
        for rank, pd_bins_rank in pd_bins.groupby('rank'):
            gs_pd_bins_rank = gs_pd_bins[gs_pd_bins['rank'] == rank]
            if gs_pd_bins_rank.empty:
                continue

            precision_bp_rows = pd_bins_rank[pd_bins_rank['purity_bp'].notnull()]['purity_bp']
            precision_seq_rows = pd_bins_rank[pd_bins_rank['purity_seq'].notnull()]['purity_seq']
            recall_bp_rows = pd_bins_rank[pd_bins_rank['true_size'] > 0]['completeness_bp']
            recall_seq_rows = pd_bins_rank[pd_bins_rank['true_size'] > 0]['completeness_seq']

            avg_precision_bp = precision_bp_rows.mean()
            avg_recall_bp = recall_bp_rows.mean()
            f1_score_bp = 2 * avg_precision_bp * avg_recall_bp / (avg_precision_bp + avg_recall_bp)

            avg_precision_seq = precision_seq_rows.mean()
            avg_recall_seq = recall_seq_rows.mean()
            f1_score_seq = 2 * avg_precision_seq * avg_recall_seq / (avg_precision_seq + avg_recall_seq)

            precision_bp, recall_bp, accuracy_bp, misclassification_rate_bp,\
                precision_seq, recall_seq, accuracy_seq, misclassification_rate_seq,\
                percentage_of_assigned_bps = compute_metrics_over_bins(rank, gs_pd_bins_rank, pd_bins_rank, query)

            bin_ids = pd_bins_rank['id'][pd_bins_rank['id'].notnull()].tolist()
            ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq = rand_index.compute_metrics(bin_ids, query)

            genome_recovery_val = genome_recovery.calc_dict(pd_bins_rank, min_completeness, max_contamination)

            df = pd.DataFrame(OrderedDict([(utils_labels.TOOL, query.label),
                                           (utils_labels.BINNING_TYPE, query.binning_type),
                                           (utils_labels.SAMPLE, sample_id),
                                           (utils_labels.RANK, rank),
                                           (utils_labels.AVG_PRECISION_BP, [avg_precision_bp]),
                                           (utils_labels.AVG_PRECISION_BP_SEM, [precision_bp_rows.sem()]),
                                           (utils_labels.AVG_RECALL_BP, [avg_recall_bp]),
                                           (utils_labels.AVG_RECALL_BP_SEM, [recall_bp_rows.sem()]),
                                           (utils_labels.F1_SCORE_BP, [f1_score_bp]),

                                           (utils_labels.AVG_PRECISION_SEQ, [avg_precision_seq]),
                                           (utils_labels.AVG_PRECISION_SEQ_SEM, [precision_seq_rows.sem()]),
                                           (utils_labels.AVG_RECALL_SEQ, [avg_recall_seq]),
                                           (utils_labels.AVG_RECALL_SEQ_SEM, [recall_seq_rows.sem()]),
                                           (utils_labels.F1_SCORE_SEQ, [f1_score_seq]),

                                           (utils_labels.PRECISION_PER_BP, [precision_bp]),
                                           (utils_labels.PRECISION_PER_SEQ, [precision_seq]),
                                           (utils_labels.RECALL_PER_BP, [recall_bp]),
                                           (utils_labels.RECALL_PER_SEQ, [recall_seq]),
                                           (utils_labels.F1_SCORE_PER_BP, [2 * precision_bp * recall_bp / (precision_bp + recall_bp)]),
                                           (utils_labels.F1_SCORE_PER_SEQ, [2 * precision_seq * recall_seq / (precision_seq + recall_seq)]),

                                           (utils_labels.ACCURACY_PER_BP, [accuracy_bp]),
                                           (utils_labels.ACCURACY_PER_SEQ, [accuracy_seq]),

                                           (utils_labels.PERCENTAGE_ASSIGNED_BPS, [percentage_of_assigned_bps]),
                                           (utils_labels.PERCENTAGE_ASSIGNED_SEQS, [percentage_of_assigned_seqs[rank]]),
                                           (utils_labels.RI_BY_BP, [ri_by_bp]),
                                           (utils_labels.RI_BY_SEQ, [ri_by_seq]),
                                           (utils_labels.ARI_BY_BP, [ari_by_bp]),
                                           (utils_labels.ARI_BY_SEQ, [ari_by_seq]),

                                           (utils_labels.UNIFRAC_BP, [unifrac_bp]),
                                           (utils_labels.UNIFRAC_SEQ, [unifrac_seq]),

                                           (utils_labels.MISCLASSIFICATION_PER_BP, [misclassification_rate_bp]),
                                           (utils_labels.MISCLASSIFICATION_PER_SEQ, [misclassification_rate_seq])]))

            df_genome_recovery = pd.DataFrame.from_dict(genome_recovery_val, orient='index').T
            df = df.join(df_genome_recovery)
            df_summary = pd.concat([df_summary, df], ignore_index=True, sort=True)
            pd_bins_all['sample_id'] = sample_id
        logging.getLogger('amber').info('done')
    return df_summary, pd_bins_all


def plot_heat_maps(sample_id_to_queries_list, output_dir):
    for sample_id in sample_id_to_queries_list:
        for query in sample_id_to_queries_list[sample_id]:
            if isinstance(query, binning_classes.GenomeQuery):
                df_confusion = precision_recall_per_bin.transform_confusion_matrix(query)
                plots.plot_heatmap(df_confusion, sample_id, os.path.join(output_dir, 'genome', query.label))


def plot_genome_binning(sample_id_to_queries_list, df_summary, pd_bins, plot_heatmaps, output_dir):
    df_summary_g = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'genome']
    if len(df_summary_g) == 0:
        return

    logging.getLogger('amber').info('Creating genome binning plots...')

    if plot_heatmaps:
        plot_heat_maps(sample_id_to_queries_list, output_dir)

    plots.create_legend(df_summary_g, output_dir)
    plots.plot_avg_precision_recall(df_summary_g, output_dir)
    plots.plot_precision_recall(df_summary_g, output_dir)
    plots.plot_adjusted_rand_index_vs_assigned_bps(df_summary_g, output_dir)

    pd_bins_g = pd_bins[pd_bins['rank'] == 'NA']
    plots.plot_boxplot(pd_bins_g, 'purity_bp', output_dir)
    plots.plot_boxplot(pd_bins_g, 'completeness_bp', output_dir)

    plot_by_genome.plot_precision_recall_per_bin(pd_bins_g, output_dir)

    plots.plot_contamination(pd_bins[pd_bins['rank'] == 'NA'], 'genome', 'Contamination', 'Index of bin (sorted by contamination (bp))', 'Contamination (bp)', plots.create_contamination_column, output_dir)
    plots.plot_contamination(pd_bins[pd_bins['rank'] == 'NA'], 'genome', 'Completeness - contamination', 'Index of bin (sorted by completeness - contamination (bp))', 'Completeness - contamination (bp)', plots.create_completeness_minus_contamination_column, output_dir)
    logging.getLogger('amber').info('done')


def plot_taxonomic_binning(df_summary, pd_bins, output_dir):
    df_summary_t = df_summary[df_summary[utils_labels.BINNING_TYPE] == 'taxonomic']
    if len(df_summary_t) == 0:
        return

    logging.getLogger('amber').info('Creating taxonomic binning plots...')

    for rank, pd_group in df_summary_t.groupby('rank'):
        plots.plot_avg_precision_recall(pd_group, output_dir, rank)
        plots.plot_precision_recall(pd_group, output_dir, rank)
        plots.plot_adjusted_rand_index_vs_assigned_bps(pd_group, output_dir, rank)

    metrics_list = [utils_labels.AVG_PRECISION_BP, utils_labels.AVG_RECALL_BP]
    errors_list = [utils_labels.AVG_PRECISION_BP_SEM, utils_labels.AVG_RECALL_BP_SEM]
    plots.plot_taxonomic_results(df_summary_t, metrics_list, errors_list, 'avg_precision_recall_bp', output_dir)

    metrics_list = [utils_labels.AVG_PRECISION_SEQ, utils_labels.AVG_RECALL_SEQ]
    errors_list = [utils_labels.AVG_PRECISION_SEQ_SEM, utils_labels.AVG_RECALL_SEQ_SEM]
    plots.plot_taxonomic_results(df_summary_t, metrics_list, errors_list, 'avg_precision_recall_seq', output_dir)

    metrics_list = [utils_labels.PRECISION_PER_BP, utils_labels.RECALL_PER_BP, utils_labels.PRECISION_PER_SEQ, utils_labels.RECALL_PER_SEQ]
    plots.plot_taxonomic_results(df_summary_t, metrics_list, [], 'precision_recall', output_dir)

    for rank in load_ncbi_taxinfo.RANKS:
        pd_bins_rank = pd_bins[pd_bins['rank'] == rank]
        plots.plot_contamination(pd_bins_rank, 'taxonomic', rank + ' | Contamination', 'Index of bin (sorted by contamination (bp))', 'Contamination (bp)', plots.create_contamination_column, output_dir)
        plots.plot_contamination(pd_bins_rank, 'taxonomic', rank + ' | Completeness - contamination', 'Index of bin (sorted by completeness - contamination (bp))', 'Completeness - contamination (bp)', plots.create_completeness_minus_contamination_column, output_dir)

    logging.getLogger('amber').info('done')


def main(args=None):
    parser = argparse.ArgumentParser(description="AMBER: Assessment of Metagenome BinnERs",
                                     parents=[argparse_parents.PARSER_MULTI2], prog='AMBER')
    parser.add_argument('-p', '--filter', help=argparse_parents.HELP_FILTER)
    parser.add_argument('-n', '--min_length', help="Minimum length of sequences", type=int, required=False)
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    parser.add_argument('--stdout', help="Print summary to stdout", action='store_true')
    parser.add_argument('-d', '--desc', help="Description for HTML page", required=False)
    parser.add_argument('--silent', help='Silent mode', action='store_true')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    group_g = parser.add_argument_group('genome binning-specific arguments')
    group_g.add_argument('-m', '--map_by_completeness', help=argparse_parents.HELP_MAP_BY_RECALL, action='store_true')
    group_g.add_argument('-x', '--min_completeness', help=argparse_parents.HELP_THRESHOLDS_COMPLETENESS, required=False)
    group_g.add_argument('-y', '--max_contamination', help=argparse_parents.HELP_THRESHOLDS_CONTAMINATION, required=False)
    group_g.add_argument('-c', '--plot_heatmaps', help="Plot heatmaps of confusion matrices (can take some minutes)", action='store_true')
    group_g.add_argument('-r', '--remove_genomes', help=argparse_parents.HELP_GENOMES_FILE, required=False)
    group_g.add_argument('-k', '--keyword', help=argparse_parents.HELP_KEYWORD, required=False)

    group_t = parser.add_argument_group('taxonomic binning-specific arguments')
    group_t.add_argument('--ncbi_nodes_file', help="NCBI nodes file", required=False)
    group_t.add_argument('--ncbi_names_file', help="NCBI names file", required=False)
    group_t.add_argument('--ncbi_merged_file', help="NCBI merged file", required=False)
    group_t.add_argument('--rank_as_genome_binning', help="Assess taxonomic binning at a rank also as genome binning. Valid ranks: superkingdom, phylum, class, order, family, genus, species, strain", required=False)

    args = parser.parse_args(args)
    output_dir = os.path.abspath(args.output_dir)
    logger = get_logger(output_dir, args.silent)

    min_completeness = None
    max_contamination = None
    if args.min_completeness:
        min_completeness = [int(x.strip())/100.0 for x in args.min_completeness.split(',')]
    if args.max_contamination:
        max_contamination = [int(x.strip())/100.0 for x in args.max_contamination.split(',')]

    labels = get_labels(args.labels, args.bin_files)

    genome_to_unique_common = load_data.load_unique_common(args.remove_genomes)

    options = binning_classes.Options(filter_tail_percentage=args.filter,
                                      genome_to_unique_common=genome_to_unique_common,
                                      filter_keyword=args.keyword,
                                      map_by_completeness=args.map_by_completeness,
                                      min_length=args.min_length,
                                      rank_as_genome_binning=args.rank_as_genome_binning)

    load_data.load_ncbi_info(args.ncbi_nodes_file, args.ncbi_names_file, args.ncbi_merged_file)

    sample_id_to_gs_list, sample_ids_list, sample_id_to_num_genomes, sample_id_to_queries_list = \
        load_data.load_queries(args.gold_standard_file,
                               args.bin_files,
                               options,
                               labels)

    create_output_directories(output_dir, sample_id_to_queries_list)

    df_summary_gs, pd_bins_gs = evaluate_samples_queries(sample_id_to_gs_list,
                                                         min_completeness, max_contamination)
    df_summary, pd_bins = evaluate_samples_queries(sample_id_to_queries_list,
                                                   min_completeness, max_contamination)

    logger.info('Saving computed metrics...')
    df_summary.to_csv(os.path.join(output_dir, 'results.tsv'), sep='\t', index=False, float_format='%.3f')
    if args.stdout:
        summary_columns = [utils_labels.TOOL] + [col for col in df_summary if col != utils_labels.TOOL]
        print(df_summary[summary_columns].to_string(index=False))
    pd_bins_g = pd_bins[pd_bins['rank'] == 'NA']

    for tool, pd_group in pd_bins_g.groupby(utils_labels.TOOL):
        bins_columns = amber_html.get_genome_bins_columns()
        table = pd_group[['sample_id'] + list(bins_columns.keys())].rename(columns=dict(bins_columns))
        table.to_csv(os.path.join(output_dir, 'genome', tool, 'metrics_per_bin.tsv'), sep='\t', index=False)

    pd_bins_t = pd_bins[pd_bins['rank'] != 'NA']
    for tool, pd_group in pd_bins_t.groupby(utils_labels.TOOL):
        bins_columns = amber_html.get_tax_bins_columns()
        if pd_group['name'].isnull().any():
            del bins_columns['name']
        table = pd_group[['sample_id'] + list(bins_columns.keys())].rename(columns=dict(bins_columns))
        table.to_csv(os.path.join(output_dir, 'taxonomic', tool, 'metrics_per_bin.tsv'), sep='\t', index=False)
    logger.info('done')

    plot_genome_binning(sample_id_to_queries_list, df_summary, pd_bins, args.plot_heatmaps, output_dir)
    plot_taxonomic_binning(df_summary, pd_bins, output_dir)

    amber_html.create_html(sample_id_to_num_genomes,
                           pd.concat([df_summary_gs, df_summary], ignore_index=True),
                           pd_bins,
                           [utils_labels.GS] + labels,
                           sample_ids_list,
                           args.output_dir,
                           args.desc)
    logger.info('AMBER finished successfully. All results have been saved to {}'.format(output_dir))


if __name__ == "__main__":
    main()
