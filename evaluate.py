#!/usr/bin/env python

import argparse
import os
import sys
import collections
import precision_recall_per_bin
import precision_recall_average
import precision_recall_by_bpcount
import rand_index
import genome_recovery
import accuracy
import plot_by_genome
import html_plots
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from utils import exclude_genomes
from utils import load_data
from utils import argparse_parents
from utils import labels


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
                 labels,
                 filter_tail_percentage,
                 genomes_file,
                 keyword,
                 map_by_recall,
                 output_dir):
    gold_standard = load_data.get_genome_mapping(gold_standard_file, fasta_file)
    labels_iterator = iter(labels)
    summary_per_query = []
    for query_file in query_files:
        tool_id = query_file.split('/')[-1]
        binning_label = next(labels_iterator)
        path = os.path.join(output_dir, tool_id)
        load_data.make_sure_path_exists(path)

        query = load_data.open_query(query_file)

        # PRECISION RECALL PER BIN
        bin_metrics = precision_recall_per_bin.compute_metrics(query, gold_standard, map_by_recall)
        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)
        f = open(path + "/precision_recall.tsv", 'w')
        precision_recall_per_bin.print_metrics(bin_metrics, f)
        # slow code disabled
        # plot_by_genome.plot_by_genome(bin_metrics, path + '/genomes_sorted_by_recall', 'recall')
        # plot_by_genome.plot_by_genome(bin_metrics, path + '/genomes_sorted_by_precision', 'precision')
        f.close()

        # AVG PRECISION RECALL
        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, sem_precision, sem_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_percentage)
        f = open(path + "/precision_recall_avg.tsv", 'w')
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
        f = open(path + "/precision_recall_by_bpcount.tsv", 'w')
        precision_recall_by_bpcount.print_precision_recall_by_bpcount(precision, recall, f)
        f.close()

        # (ADJUSTED) RAND INDEX
        ri_by_seq, ri_by_bp, a_rand_index_by_bp, a_rand_index_by_seq, percent_assigned_bps = rand_index.compute_metrics(query, gold_standard)
        f = open(path + "/rand_index.tsv", 'w')
        rand_index.print_rand_indices(ri_by_seq, ri_by_bp, a_rand_index_by_bp, a_rand_index_by_seq, percent_assigned_bps, f)
        f.close()

        # GENOME RECOVERY
        genome_recovery_val = genome_recovery.calc_table(bin_metrics)

        # ACCURACY
        acc = accuracy.compute_metrics(query, gold_standard)
        #acc = 0.0

        summary_per_query.append((collections.OrderedDict([('binning_label', binning_label),
                                   ('avg_precision', avg_precision),
                                   ('std_deviation_precision', std_deviation_precision),
                                   ('sem_precision', sem_precision),
                                   ('avg_recall', avg_recall),
                                   ('std_deviation_recall', std_deviation_recall),
                                   ('sem_recall', sem_recall),
                                   ('precision', precision),
                                   ('recall', recall),
                                   ('ri_by_bp', ri_by_bp),
                                   ('ri_by_seq', ri_by_seq),
                                   ('a_rand_index_by_bp', a_rand_index_by_bp),
                                   ('a_rand_index_by_seq', a_rand_index_by_seq),
                                   ('percent_assigned_bps', percent_assigned_bps),
                                   ('accuracy', acc),
                                   ('_05compl_01cont', genome_recovery_val[5]),
                                   ('_07compl_01cont', genome_recovery_val[3]),
                                   ('_09compl_01cont', genome_recovery_val[1]),
                                   ('_05compl_005cont', genome_recovery_val[4]),
                                   ('_07compl_005cont', genome_recovery_val[2]),
                                   ('_09compl_005cont', genome_recovery_val[0])]),
                                  bin_metrics))
    return summary_per_query


def convert_summary_to_tuples_of_strings(summary_per_query):
    tuples = []
    for summary in summary_per_query:
        tuples.append(((summary['binning_label']),
                      format(summary['avg_precision'], '.3f'),
                      format(summary['std_deviation_precision'], '.3f'),
                      format(summary['sem_precision'], '.3f'),
                      format(summary['avg_recall'], '.3f'),
                      format(summary['std_deviation_recall'], '.3f'),
                      format(summary['sem_recall'], '.3f'),
                      format(summary['precision'], '.3f'),
                      format(summary['recall'], '.3f'),
                      format(summary['ri_by_bp'], '.3f'),
                      format(summary['ri_by_seq'], '.3f'),
                      format(summary['a_rand_index_by_bp'], '.3f'),
                      format(summary['a_rand_index_by_seq'], '.3f'),
                      format(summary['percent_assigned_bps'], '.3f'),
                      format(summary['accuracy'], '.3f'),
                      str(summary['_05compl_01cont']),
                      str(summary['_07compl_01cont']),
                      str(summary['_09compl_01cont']),
                      str(summary['_05compl_005cont']),
                      str(summary['_07compl_005cont']),
                      str(summary['_09compl_005cont'])))
    return tuples


def create_legend(summary_per_query, output_dir):
    colors_iter = iter(plot_by_genome.create_colors_list())
    labels = []
    circles = []
    for summary in summary_per_query:
        labels.append(summary['binning_label'])
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=10, markerfacecolor=next(colors_iter)))

    fig = plt.figure(figsize=(0.5, 0.5))
    fig.legend(circles, labels, loc='center', frameon=False, ncol=3, handletextpad=0.1)
    fig.savefig(os.path.normpath(output_dir + '/legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


def plot_summary(summary_per_query, output_dir, plot_type, file_name, xlabel, ylabel):
    colors_list = plot_by_genome.create_colors_list()
    if len(summary_per_query) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axis to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    i = 0
    plot_labels = []
    if plot_type == 'e':
        for summary in summary_per_query:
            axs.errorbar(summary['avg_precision'], summary['avg_recall'], xerr=summary['sem_precision'], yerr=summary['sem_recall'],
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
            plot_labels.append(summary['binning_label'])
            i += 1
    if plot_type == 'w':
        for summary in summary_per_query:
            axs.plot(summary['precision'], summary['recall'], marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(summary['binning_label'])
            i += 1
    elif plot_type == 'p':
        for summary in summary_per_query:
            axs.plot(summary['a_rand_index_by_bp'], summary['percent_assigned_bps'], marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(summary['binning_label'])
            i += 1

    # turn on grid
    axs.minorticks_on()
    axs.grid(which='major', linestyle='-', linewidth='0.5')
    axs.grid(which='minor', linestyle=':', linewidth='0.5')

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}%'.format(x * 100) for x in vals])
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')
    lgd = plt.legend(plot_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)

    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def plot_avg_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'e',
                 'avg_precision_recall',
                 'Average precision per bin',
                 'Average recall per genome')


def plot_weighed_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'w',
                 'weighed_precision_recall',
                 'Precision per base pair',
                 'Recall per base pair')


def plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'p',
                 'ari_vs_assigned_bps',
                 'Adjusted Rand Index',
                 'Percentage of assigned base pairs')


def print_summary(summary_per_query, output_dir=None):
    if output_dir is None:
        stream=sys.stdout
    else:
        stream = open(output_dir + "/summary.tsv", 'w')
    stream.write("%s\n" % "\t".join((labels.TOOL,
                                     labels.AVG_PRECISION,
                                     labels.STD_DEV_PRECISION,
                                     labels.SEM_PRECISION,
                                     labels.AVG_RECALL,
                                     labels.STD_DEV_RECALL,
                                     labels.SEM_RECALL,
                                     labels.PRECISION,
                                     labels.RECALL,
                                     labels.RI_BY_BP,
                                     labels.RI_BY_SEQ,
                                     labels.ARI_BY_BP,
                                     labels.ARI_BY_SEQ,
                                     labels.PERCENTAGE_ASSIGNED_BPS,
                                     labels.ACCURACY,
                                     ">0.5compl<0.1cont",
                                     ">0.7compl<0.1cont",
                                     ">0.9compl<0.1cont",
                                     ">0.5compl<0.05cont",
                                     ">0.7compl<0.05cont",
                                     ">0.9compl<0.05cont")))
    for summary in summary_per_query:
        stream.write("%s\n" % "\t".join(summary))
    if output_dir is not None:
        stream.close()


def compute_rankings(summary_per_query, output_dir):
    f = open(os.path.normpath(output_dir + '/rankings.txt'), 'w')
    f.write("Tool\tAverage precision\n")
    sorted_by = sorted(summary_per_query, key=lambda x: x['avg_precision'], reverse=True)
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary['binning_label'], summary['avg_precision']))

    sorted_by = sorted(summary_per_query, key=lambda x: x['avg_recall'], reverse=True)
    f.write("\nTool\tAverage recall\n")
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary['binning_label'], summary['avg_recall']))

    sorted_by = sorted(summary_per_query, key=lambda x: x['avg_precision'] + x['avg_recall'], reverse=True)
    f.write("\nTool\tAverage precision + Average recall\tAverage precision\tAverage recall\n")
    for summary in sorted_by:
        f.write("%s\t%1.3f\t%1.3f\t%1.3f\n" % (summary['binning_label'],
                                               summary['avg_precision'] + summary['avg_recall'],
                                               summary['avg_precision'],
                                               summary['avg_recall']))
    f.close()


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics and figures for one or more binning files; output summary to screen and results per binning file to chosen directory",
                                     parents=[argparse_parents.PARSER_MULTI2])
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    parser.add_argument('-m', '--map_by_recall', help=argparse_parents.HELP_MAP_BY_RECALL, action='store_true')
    args = parser.parse_args()
    binning_labels = get_labels(args.labels, args.bin_files)
    summary_per_query = evaluate_all(args.gold_standard_file,
                                     args.fasta_file,
                                     args.bin_files,
                                     binning_labels,
                                     args.filter,
                                     args.remove_genomes,
                                     args.keyword,
                                     args.map_by_recall,
                                     args.output_dir)
    summary_dict = [x[0] for x in summary_per_query]
    summary_as_string = convert_summary_to_tuples_of_strings(summary_dict)
    print_summary(summary_as_string)
    print_summary(summary_as_string, args.output_dir)
    create_legend(summary_dict, args.output_dir)
    plot_avg_precision_recall(summary_dict, args.output_dir)
    plot_weighed_precision_recall(summary_dict, args.output_dir)
    plot_adjusted_rand_index_vs_assigned_bps(summary_dict, args.output_dir)
    plot_by_genome.plot_by_genome2(summary_per_query, args.output_dir)
    compute_rankings(summary_dict, args.output_dir)

    precision_recall_files = []
    for query_file in args.bin_files:
        tool_id = query_file.split('/')[-1]
        precision_recall_files.append(os.path.join(args.output_dir, tool_id) + "/precision_recall.tsv")
    df = pd.DataFrame.from_dict(summary_dict)
    df.set_index('binning_label', inplace=True)
    df.rename(columns={'binning_label': 'Tool'}, inplace=True)
    html_plots.build_html(precision_recall_files,
                          binning_labels,
                          df,
                          os.path.join(args.output_dir, "summary.html"))


if __name__ == "__main__":
    main()
