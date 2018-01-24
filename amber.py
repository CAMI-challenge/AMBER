#!/usr/bin/env python

import argparse
import collections
import os
import matplotlib
from version import __version__
from src import accuracy
from src import genome_recovery
from src import html_plots
from src import plot_by_genome
from src import plots
from src import precision_recall_by_bpcount
from src import precision_recall_per_bin
from src import rand_index
from src import precision_recall_average
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from src.utils import exclude_genomes
from src.utils import load_data
from src.utils import argparse_parents
from src.utils import labels


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


def evaluate_all(gold_standard_file,
                 fasta_file,
                 query_files,
                 binning_labels,
                 filter_tail_percentage,
                 filter_genomes_file,
                 keyword,
                 map_by_recall,
                 min_completeness, max_contamination,
                 output_dir):
    gold_standard = load_data.get_genome_mapping(gold_standard_file, fasta_file)
    labels_iterator = iter(binning_labels)
    summary_per_query = []
    bin_metrics_per_query = []
    count = 0
    for query_file in query_files:
        tool_id = query_file.split('/')[-1]
        binning_label = next(labels_iterator)
        path = os.path.join(output_dir, tool_id)
        load_data.make_sure_path_exists(path)

        f = open(os.path.join(path, "label.txt"), 'w')
        f.write("#({}){}".format(count, binning_label))
        f.close()

        query = load_data.open_query(query_file)

        if map_by_recall:
            bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes = precision_recall_per_bin.map_genomes_by_recall(gold_standard, query.bin_id_to_list_of_sequence_id)
        else:
            bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes = precision_recall_per_bin.map_genomes(gold_standard, query.bin_id_to_list_of_sequence_id)

        df_confusion = precision_recall_per_bin.compute_confusion_matrix(
            bin_id_to_mapped_genome,
            bin_id_to_genome_id_to_total_length,
            gold_standard,
            query)
        plots.plot_heatmap(df_confusion, path)

        # PRECISION RECALL PER BIN
        bin_metrics = precision_recall_per_bin.compute_metrics(query, gold_standard, bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped_genomes)
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
        precision_recall_average.print_precision_recall(binning_label,
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        sem_precision,
                                                        sem_recall,
                                                        f)
        f.close()

        # PRECISION RECALL BY BP COUNTS
        precision, recall = precision_recall_by_bpcount.compute_metrics(query, gold_standard)
        f = open(os.path.join(path, "purity_completeness_by_bpcount.tsv"), 'w')
        precision_recall_by_bpcount.print_precision_recall_by_bpcount(precision, recall, f)
        f.close()

        # (ADJUSTED) RAND INDEX
        ri_by_seq, ri_by_bp, a_rand_index_by_bp, a_rand_index_by_seq, percent_assigned_bps = rand_index.compute_metrics(query, gold_standard)
        f = open(os.path.join(path, "rand_index.tsv"), 'w')
        rand_index.print_rand_indices(ri_by_seq, ri_by_bp, a_rand_index_by_bp, a_rand_index_by_seq, percent_assigned_bps, f)
        f.close()

        # GENOME RECOVERY
        genome_recovery_val = genome_recovery.calc_dict(bin_metrics, min_completeness, max_contamination)

        # ACCURACY
        acc = accuracy.compute_metrics(query, gold_standard)

        summary_per_query.append(collections.OrderedDict([(labels.TOOL, binning_label),
                                                          (labels.AVG_PRECISION, avg_precision),
                                                          (labels.STD_DEV_PRECISION, std_deviation_precision),
                                                          (labels.SEM_PRECISION, sem_precision),
                                                          (labels.AVG_RECALL, avg_recall),
                                                          (labels.STD_DEV_RECALL, std_deviation_recall),
                                                          (labels.SEM_RECALL, sem_recall),
                                                          (labels.PRECISION, precision),
                                                          (labels.RECALL, recall),
                                                          (labels.RI_BY_BP, ri_by_bp),
                                                          (labels.RI_BY_SEQ, ri_by_seq),
                                                          (labels.ARI_BY_BP, a_rand_index_by_bp),
                                                          (labels.ARI_BY_SEQ, a_rand_index_by_seq),
                                                          (labels.PERCENTAGE_ASSIGNED_BPS, percent_assigned_bps),
                                                          (labels.ACCURACY, acc)] +
                                                         [(k, v) for k, v in genome_recovery_val.items()]))
        bin_metrics_per_query.append(bin_metrics)
        count += 1
    return summary_per_query, bin_metrics_per_query


def create_legend(summary_per_query, output_dir):
    colors_iter = iter(plots.create_colors_list())
    binning_labels = []
    circles = []
    for summary in summary_per_query:
        binning_labels.append(summary[labels.TOOL])
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=10, markerfacecolor=next(colors_iter)))

    fig = plt.figure(figsize=(0.5, 0.5))
    fig.legend(circles, binning_labels, loc='center', frameon=False, ncol=5, handletextpad=0.1)
    fig.savefig(os.path.normpath(output_dir + '/legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


def compute_rankings(summary_per_query, output_dir):
    f = open(os.path.normpath(output_dir + '/rankings.txt'), 'w')
    f.write("Tool\tAverage purity\n")
    sorted_by = sorted(summary_per_query, key=lambda x: x[labels.AVG_PRECISION], reverse=True)
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary[labels.TOOL], summary[labels.AVG_PRECISION]))

    sorted_by = sorted(summary_per_query, key=lambda x: x[labels.AVG_RECALL], reverse=True)
    f.write("\nTool\tAverage completeness\n")
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary[labels.TOOL], summary[labels.AVG_RECALL]))

    sorted_by = sorted(summary_per_query, key=lambda x: x[labels.AVG_PRECISION] + x[labels.AVG_RECALL], reverse=True)
    f.write("\nTool\tAverage purity + Average completeness\tAverage purity\tAverage completeness\n")
    for summary in sorted_by:
        f.write("%s\t%1.3f\t%1.3f\t%1.3f\n" % (summary[labels.TOOL],
                                               summary[labels.AVG_PRECISION] + summary[labels.AVG_RECALL],
                                               summary[labels.AVG_PRECISION],
                                               summary[labels.AVG_RECALL]))
    f.close()


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics and figures for one or more binning files; output summary to screen and results per binning file to chosen directory",
                                     parents=[argparse_parents.PARSER_MULTI2], prog='AMBER')
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    parser.add_argument('-m', '--map_by_completeness', help=argparse_parents.HELP_MAP_BY_RECALL, action='store_true')
    parser.add_argument('-x', '--min_completeness', help=argparse_parents.HELP_THRESHOLDS_COMPLETENESS, required=False)
    parser.add_argument('-y', '--max_contamination', help=argparse_parents.HELP_THRESHOLDS_CONTAMINATION, required=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()

    min_completeness = None
    max_contamination = None
    if args.min_completeness:
        min_completeness = [int(x.strip())/100.0 for x in args.min_completeness.split(',')]
    if args.max_contamination:
        max_contamination = [int(x.strip())/100.0 for x in args.max_contamination.split(',')]

    binning_labels = get_labels(args.labels, args.bin_files)
    summary_per_query, bin_metrics_per_query = evaluate_all(args.gold_standard_file,
                                                            args.fasta_file,
                                                            args.bin_files,
                                                            binning_labels,
                                                            args.filter,
                                                            args.remove_genomes,
                                                            args.keyword,
                                                            args.map_by_completeness,
                                                            min_completeness, max_contamination,
                                                            args.output_dir)
    df = pd.DataFrame.from_dict(summary_per_query)
    print(df.to_csv(sep='\t', index=False, float_format='%.3f'))
    df.to_csv(path_or_buf=os.path.join(args.output_dir, "summary.tsv"), sep='\t', index=False, float_format='%.3f')

    create_legend(summary_per_query, args.output_dir)
    plots.plot_avg_precision_recall(summary_per_query, args.output_dir)
    plots.plot_weighed_precision_recall(summary_per_query, args.output_dir)
    plots.plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, args.output_dir)

    plots.plot_boxplot(bin_metrics_per_query, binning_labels, 'purity', args.output_dir)
    plots.plot_boxplot(bin_metrics_per_query, binning_labels, 'completeness', args.output_dir)

    plot_by_genome.plot_by_genome2(bin_metrics_per_query, binning_labels, args.output_dir)
    compute_rankings(summary_per_query, args.output_dir)

    precision_recall_files = []
    for query_file in args.bin_files:
        tool_id = query_file.split('/')[-1]
        precision_recall_files.append(os.path.join(args.output_dir, tool_id, "purity_completeness.tsv"))
    df = pd.DataFrame.from_dict(summary_per_query)
    df.rename(columns={labels.TOOL: 'Tool'}, inplace=True)
    df.set_index('Tool', inplace=True)
    html_plots.build_html(precision_recall_files,
                          binning_labels,
                          df,
                          os.path.join(args.output_dir, "summary.html"))


if __name__ == "__main__":
    main()
