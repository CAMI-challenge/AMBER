#!/usr/bin/env python

import argparse
import collections
import os
import errno
import matplotlib
from version import __version__
from src import accuracy
from src import genome_recovery
from src import html_plots
from src import plot_by_genome
from src import plots
from src import precision_recall_per_bin
from src import rand_index
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from src.utils import exclude_genomes
from src.utils import load_data
from src.utils import argparse_parents
from src.utils import labels as utils_labels
from src.utils import load_ncbi_taxinfo
from src.utils import filter_tail


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def create_output_directories(output_dir, labels):
    make_sure_path_exists(output_dir)
    for label in labels:
        make_sure_path_exists(os.path.join(output_dir, label))


def get_labels(labels, bin_files):
    if labels:
        labels_list = [x.strip() for x in labels.split(',')]
        if len(labels_list) != len(bin_files):
            exit('Number of labels does not match the number of binning files. Please check parameter -l, --labels.')
        return labels_list
    tool_id = []
    for bin_file in bin_files:
        tool_id.append(bin_file.split('/')[-1])
    return tool_id


def load_queries(gold_standard_file, fastx_file, query_files, map_by_completeness, ncbi_nodes_file, min_length, labels):
    if not min_length:
        min_length = 0

    if ncbi_nodes_file:
        tax_id_to_parent, tax_id_to_rank = load_ncbi_taxinfo.load_tax_info(ncbi_nodes_file)
    else:
        tax_id_to_parent = tax_id_to_rank = None

    g_gold_standard, t_gold_standard = load_data.open_query(gold_standard_file,
                                                            True,
                                                            fastx_file,
                                                            tax_id_to_parent, tax_id_to_rank,
                                                            None,
                                                            min_length)
    gold_standard = load_data.GoldStandard(g_gold_standard, t_gold_standard)

    queries_list = []
    for query_file, label in zip(query_files, labels):
        g_query, t_query = load_data.open_query(query_file,
                                                False,
                                                None,
                                                tax_id_to_parent, tax_id_to_rank,
                                                gold_standard,
                                                0)
        if g_query:
            g_query.label = label
            g_query.map_by_completeness = map_by_completeness
            queries_list.append(g_query)
        if t_query:
            t_query.label = label
            queries_list.append(t_query)

    # TODO if there is a g_query (t_query), there must be a g_gold_standard (t_gold_standard)

    return gold_standard, queries_list


def evaluate_all(gold_standard,
                 queries_list,
                 filter_tail_percentage,
                 filter_genomes_file,
                 keyword,
                 min_completeness, max_contamination,
                 output_dir):

    if gold_standard.genome_query:
        gs_genome_bins_metrics = gold_standard.genome_query.get_bins_metrics(gold_standard)
        gs_pd_genome_bins = pd.DataFrame.from_dict(gs_genome_bins_metrics)
        gs_pd_genome_bins['rank'] = 'NA'

    if gold_standard.taxonomic_query:
        gs_tax_bins_metrics = gold_standard.taxonomic_query.get_bins_metrics(gold_standard)
        gs_pd_tax_bins = pd.DataFrame.from_dict(gs_tax_bins_metrics)

    for query in queries_list:
        query.compute_true_positives(gold_standard)

        precision_recall_per_bin.compute_precision_recall(gold_standard, query)

        bins_metrics = query.get_bins_metrics(gold_standard)

        if isinstance(query, load_data.GenomeQuery):
            if filter_tail_percentage:
                filter_tail.filter_tail(bins_metrics, filter_tail_percentage)
            if filter_genomes_file:
                bins_metrics = exclude_genomes.filter_data(bins_metrics, filter_genomes_file, keyword)
            pd_bins = pd.DataFrame.from_dict(bins_metrics)
            pd_bins['rank'] = 'NA'
            gs_pd_bins = gs_pd_genome_bins
        else:
            pd_bins = pd.DataFrame.from_dict(bins_metrics)
            gs_pd_bins = gs_pd_tax_bins

        for rank, pd_bins_rank in pd_bins.groupby('rank'):
            gs_pd_bins_rank = gs_pd_bins[gs_pd_bins['rank'] == rank]
            print(rank)
            precision_rows = pd_bins_rank[pd_bins_rank['purity'].notnull()]['purity']
            recall_rows = pd_bins_rank[pd_bins_rank['real_size'] > 0]['completeness']
            avg_precision = precision_rows.mean()
            sem_precision = precision_rows.sem()
            std_precision = precision_rows.std()
            avg_recall = recall_rows.mean()
            sem_recall = recall_rows.sem()
            std_recall = recall_rows.std()

            true_positives_all_bins = pd_bins_rank['true_positives'].sum()
            all_bins_length = pd_bins_rank['predicted_size'].sum()
            precision_by_bp = float(true_positives_all_bins) / float(all_bins_length)

            # print(gs_pd_bins_rank)

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

            print("true_positives_recall {}".format(true_positives_recall))


            gs_length = gs_pd_bins_rank['real_size'].sum()
            recall_by_bp = float(true_positives_recall) / float(gs_length)

            accuracy = float(true_positives_all_bins) / float(gs_length)

            print("precision by bp:\t{}".format(precision_by_bp))
            print("recall by bp:\t{}".format(recall_by_bp))
            print("accuracy:\t{}".format(accuracy))

            # exit()
        # print(df.to_csv(sep='\t', index=False, float_format='%.3f'))
        # pd_metrics = pd.concat([pd_metrics, df], ignore_index=True)



    exit()


    summary_per_query = []
    bin_metrics_per_query = []
    count = 0

    for query, label in zip(g_queries_list, labels):
        path = os.path.join(output_dir, label)

        # PRECISION RECALL PER BIN
        bin_metrics = precision_recall_per_bin.compute_metrics(g_gold_standard, query)
        if filter_genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, filter_genomes_file, keyword)
        f = open(os.path.join(path, "purity_completeness.tsv"), 'w')
        precision_recall_per_bin.print_metrics(bin_metrics, f)
        # slow code disabled
        # plot_by_genome.plot_by_genome(bin_metrics, path + '/genomes_sorted_by_recall', 'recall')
        # plot_by_genome.plot_by_genome(bin_metrics, path + '/genomes_sorted_by_precision', 'precision')
        f.close()

        # AVG PRECISION RECALL
        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, sem_precision, sem_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_percentage)
        f = open(os.path.join(path, "purity_completeness_avg.tsv"), 'w')
        precision_recall_average.print_precision_recall_table_header(f)
        precision_recall_average.print_precision_recall(label,
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        sem_precision,
                                                        sem_recall,
                                                        f)
        f.close()

        # PRECISION RECALL BY BP COUNTS
        precision, recall = precision_recall_by_bpcount.compute_metrics(query, g_gold_standard)
        f = open(os.path.join(path, "purity_completeness_by_bpcount.tsv"), 'w')
        precision_recall_by_bpcount.print_precision_recall_by_bpcount(precision, recall, f)
        f.close()

        # (ADJUSTED) RAND INDEX
        ri_by_seq, ri_by_bp, a_rand_index_by_bp, a_rand_index_by_seq, percent_assigned_bps = rand_index.compute_metrics(query, g_gold_standard)
        f = open(os.path.join(path, "rand_index.tsv"), 'w')
        rand_index.print_rand_indices(ri_by_seq, ri_by_bp, a_rand_index_by_bp, a_rand_index_by_seq, percent_assigned_bps, f)
        f.close()

        # GENOME RECOVERY
        genome_recovery_val = genome_recovery.calc_dict(bin_metrics, min_completeness, max_contamination)

        # ACCURACY
        acc = accuracy.compute_metrics(query, g_gold_standard)

        summary_per_query.append(collections.OrderedDict([(utils_labels.TOOL, label),
                                                          (utils_labels.AVG_PRECISION, avg_precision),
                                                          (utils_labels.STD_DEV_PRECISION, std_deviation_precision),
                                                          (utils_labels.SEM_PRECISION, sem_precision),
                                                          (utils_labels.AVG_RECALL, avg_recall),
                                                          (utils_labels.STD_DEV_RECALL, std_deviation_recall),
                                                          (utils_labels.SEM_RECALL, sem_recall),
                                                          (utils_labels.PRECISION, precision),
                                                          (utils_labels.RECALL, recall),
                                                          (utils_labels.RI_BY_BP, ri_by_bp),
                                                          (utils_labels.RI_BY_SEQ, ri_by_seq),
                                                          (utils_labels.ARI_BY_BP, a_rand_index_by_bp),
                                                          (utils_labels.ARI_BY_SEQ, a_rand_index_by_seq),
                                                          (utils_labels.PERCENTAGE_ASSIGNED_BPS, percent_assigned_bps),
                                                          (utils_labels.ACCURACY, acc)] +
                                                         [(k, v) for k, v in genome_recovery_val.items()]))
        bin_metrics_per_query.append(bin_metrics)
        count += 1
    return summary_per_query, bin_metrics_per_query


def create_legend(summary_per_query, output_dir):
    colors_iter = iter(plots.create_colors_list())
    labels = []
    circles = []
    for summary in summary_per_query:
        labels.append(summary[utils_labels.TOOL])
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=10, markerfacecolor=next(colors_iter)))

    fig = plt.figure(figsize=(0.5, 0.5))
    fig.legend(circles, labels, loc='center', frameon=False, ncol=5, handletextpad=0.1)
    fig.savefig(os.path.normpath(output_dir + '/legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


def compute_rankings(summary_per_query, output_dir):
    f = open(os.path.normpath(output_dir + '/rankings.txt'), 'w')
    f.write("Tool\tAverage purity\n")
    sorted_by = sorted(summary_per_query, key=lambda x: x[utils_labels.AVG_PRECISION], reverse=True)
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary[utils_labels.TOOL], summary[utils_labels.AVG_PRECISION]))

    sorted_by = sorted(summary_per_query, key=lambda x: x[utils_labels.AVG_RECALL], reverse=True)
    f.write("\nTool\tAverage completeness\n")
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary[utils_labels.TOOL], summary[utils_labels.AVG_RECALL]))

    sorted_by = sorted(summary_per_query, key=lambda x: x[utils_labels.AVG_PRECISION] + x[utils_labels.AVG_RECALL], reverse=True)
    f.write("\nTool\tAverage purity + Average completeness\tAverage purity\tAverage completeness\n")
    for summary in sorted_by:
        f.write("%s\t%1.3f\t%1.3f\t%1.3f\n" % (summary[utils_labels.TOOL],
                                               summary[utils_labels.AVG_PRECISION] + summary[utils_labels.AVG_RECALL],
                                               summary[utils_labels.AVG_PRECISION],
                                               summary[utils_labels.AVG_RECALL]))
    f.close()


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
    output_dir = os.path.abspath(args.output_dir)
    create_output_directories(output_dir, labels)

    gold_standard, queries_list = load_queries(args.gold_standard_file,
                                               args.fasta_file,
                                               args.bin_files,
                                               args.map_by_completeness,
                                               args.ncbi_nodes_file,
                                               args.min_length,
                                               labels)

    # if args.plot_heatmaps and g_queries_list:
    #     df_confusion = precision_recall_per_bin.transform_confusion_matrix_all(g_gold_standard,
    #                                                                            g_queries_list)
    #     plots.plot_heatmap(df_confusion, g_queries_list)

    summary_per_query, bin_metrics_per_query = evaluate_all(gold_standard,
                                                            queries_list,
                                                            args.filter,
                                                            args.remove_genomes,
                                                            args.keyword,
                                                            min_completeness, max_contamination,
                                                            args.output_dir)
    exit()
    df = pd.DataFrame.from_dict(summary_per_query)
    print(df.to_csv(sep='\t', index=False, float_format='%.3f'))
    df.to_csv(path_or_buf=os.path.join(args.output_dir, "summary.tsv"), sep='\t', index=False, float_format='%.3f')

    create_legend(summary_per_query, args.output_dir)
    plots.plot_avg_precision_recall(summary_per_query, args.output_dir)
    plots.plot_weighed_precision_recall(summary_per_query, args.output_dir)
    plots.plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, args.output_dir)

    plots.plot_boxplot(bin_metrics_per_query, labels, 'purity', args.output_dir)
    plots.plot_boxplot(bin_metrics_per_query, labels, 'completeness', args.output_dir)

    plot_by_genome.plot_by_genome2(bin_metrics_per_query, labels, args.output_dir)
    compute_rankings(summary_per_query, args.output_dir)

    precision_recall_files = []
    for query_file, label in zip(args.bin_files, labels):
        precision_recall_files.append(os.path.join(args.output_dir, label, "purity_completeness.tsv"))
    df = pd.DataFrame.from_dict(summary_per_query)
    df.rename(columns={utils_labels.TOOL: 'Tool'}, inplace=True)
    df.set_index('Tool', inplace=True)
    html_plots.build_html(precision_recall_files,
                          labels,
                          df,
                          os.path.join(args.output_dir, "summary.html"))


if __name__ == "__main__":
    main()
