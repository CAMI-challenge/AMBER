#!/usr/bin/env python

import argparse
import os
import sys
import errno
import precision_recall_per_genome
import precision_recall_average
import precision_recall_by_bpcount
import rand_index
import genome_recovery
import plot_by_genome
import matplotlib.pyplot as plt
import numpy as np
from utils import exclude_genomes
from utils import load_data
from utils import argparse_parents
from utils import labels


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def evaluate_all(gold_standard_file, fasta_file, query_files, labels, filter_tail_percentage, genomes_file, keyword, output_dir):
    gold_standard = load_data.get_genome_mapping(gold_standard_file, fasta_file)
    labels_iterator = iter(labels)
    summary_per_query = []
    for query_file in query_files:
        tool_id = query_file.split('/')[-1]
        binning_label = next(labels_iterator) if len(labels) > 0 else tool_id
        path = os.path.join(output_dir, tool_id)
        make_sure_path_exists(path)

        query = load_data.open_query(query_file)

        # PRECISION RECALL PER GENOME
        bin_metrics = precision_recall_per_genome.compute_metrics(query, gold_standard)
        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)
        f = open(path + "/precision_recall.tsv", 'w')
        precision_recall_per_genome.print_metrics(bin_metrics, f)
        plot_by_genome.plot_by_genome(bin_metrics, path + '/genomes_sorted_by_recall', 'recall')
        plot_by_genome.plot_by_genome(bin_metrics, path + '/genomes_sorted_by_precision', 'precision')
        f.close()

        # AVG PRECISION RECALL
        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_percentage)
        f = open(path + "/precision_recall_avg.tsv", 'w')
        precision_recall_average.print_precision_recall_table_header(f)
        precision_recall_average.print_precision_recall(binning_label,
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        std_error_precision,
                                                        std_error_recall,
                                                        f)
        f.close()

        # PRECISION RECALL BY BP COUNTS
        precision, recall = precision_recall_by_bpcount.compute_metrics(query, gold_standard)
        f = open(path + "/precision_recall_by_bpcount.tsv", 'w')
        precision_recall_by_bpcount.print_precision_recall_by_bpcount(precision, recall, f)
        f.close()

        # (ADJUSTED) RAND INDEX
        ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq, percentage_of_assigned_bps = rand_index.compute_metrics(query, gold_standard)
        f = open(path + "/rand_index.tsv", 'w')
        rand_index.print_rand_indices(ri_by_seq, ri_by_bp, ari_by_bp, ari_by_seq, percentage_of_assigned_bps, f)
        f.close()

        # GENOME RECOVERY
        genome_recovery_val = genome_recovery.calc_table(bin_metrics)

        summary_per_query.append({'binning_label': binning_label,
                                  'avg_precision': avg_precision,
                                  'std_deviation_precision': std_deviation_precision,
                                  'std_error_precision': std_error_precision,
                                  'avg_recall': avg_recall,
                                  'std_deviation_recall': std_deviation_recall,
                                  'std_error_recall': std_error_recall,
                                  'precision': precision,
                                  'recall': recall,
                                  'ri_by_bp': ri_by_bp,
                                  'ri_by_seq': ri_by_seq,
                                  'ari_by_bp': ari_by_bp,
                                  'ari_by_seq': ari_by_seq,
                                  'percentage_of_assigned_bps': percentage_of_assigned_bps,
                                  '_05compl_01cont': genome_recovery_val[5],
                                  '_07compl_01cont': genome_recovery_val[3],
                                  '_09compl_01cont': genome_recovery_val[1],
                                  '_05compl_005cont': genome_recovery_val[4],
                                  '_07compl_005cont': genome_recovery_val[2],
                                  '_09compl_005cont': genome_recovery_val[0]})
    return summary_per_query


def convert_summary_to_tuples_of_strings(summary_per_query):
    tuples = []
    for summary in summary_per_query:
        tuples.append(((summary['binning_label']),
                      format(summary['avg_precision'], '.3f'),
                      format(summary['std_deviation_precision'], '.3f'),
                      format(summary['std_error_precision'], '.3f'),
                      format(summary['avg_recall'], '.3f'),
                      format(summary['std_deviation_recall'], '.3f'),
                      format(summary['std_error_recall'], '.3f'),
                      format(summary['precision'], '.3f'),
                      format(summary['recall'], '.3f'),
                      format(summary['ri_by_bp'], '.3f'),
                      format(summary['ri_by_seq'], '.3f'),
                      format(summary['ari_by_bp'], '.3f'),
                      format(summary['ari_by_seq'], '.3f'),
                      format(summary['percentage_of_assigned_bps'], '.3f'),
                      str(summary['_05compl_01cont']),
                      str(summary['_07compl_01cont']),
                      str(summary['_09compl_01cont']),
                      str(summary['_05compl_005cont']),
                      str(summary['_07compl_005cont']),
                      str(summary['_09compl_005cont'])))
    return tuples


def create_colors_list():
    colors_list = []
    for color in plt.cm.Set1(np.linspace(0, 1, 9)):
        colors_list.append(tuple(color))
    for color in plt.cm.Set2(np.linspace(0, 1, 8)):
        colors_list.append(tuple(color))
    for color in plt.cm.Set3(np.linspace(0, 1, 12)):
        colors_list.append(tuple(color))
    return colors_list


def plot_summary(summary_per_query, output_dir, plot_type, file_name, xlabel, ylabel):
    colors_list = create_colors_list()
    if len(summary_per_query) > len(colors_list):
        raise RuntimeError("Plot of precision and recall only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axis to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    i = 0
    plot_labels = []
    if plot_type == 'e':
        for summary in summary_per_query:
            axs.errorbar(summary['avg_precision'], summary['avg_recall'], xerr=summary['std_error_precision'], yerr=summary['std_error_recall'],
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3)
            plot_labels.append(summary['binning_label'])
            i += 1
    elif plot_type == 'p':
        for summary in summary_per_query:
            axs.plot(summary['ari_by_bp'], summary['percentage_of_assigned_bps'], marker='o', color=colors_list[i])
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

    lgd = plt.legend(plot_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')


def plot_avg_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'e',
                 'avg_precision_recall',
                 'Precision',
                 'Recall')


def plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'p',
                 'ari_vs_assigned_bps',
                 'Adjusted rand index',
                 'Percentage of assigned base pairs')


def print_summary(summary_per_query, stream=sys.stdout):
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
                                     ">0.5compl<0.1cont",
                                     ">0.7compl<0.1cont",
                                     ">0.9compl<0.1cont",
                                     ">0.5compl<0.05cont",
                                     ">0.7compl<0.05cont",
                                     ">0.9compl<0.05cont")))
    for summary in summary_per_query:
        stream.write("%s\n" % "\t".join(summary))


def compute_rankings(summary_per_query, output_dir):
    f = open(os.path.normpath(output_dir + '/rankings.txt'), 'w')
    f.write("Average precision\n")
    sorted_by = sorted(summary_per_query, key=lambda x: x['avg_precision'], reverse=True)
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary['binning_label'], summary['avg_precision']))

    sorted_by = sorted(summary_per_query, key=lambda x: x['avg_recall'], reverse=True)
    f.write("\nAverage recall\n")
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary['binning_label'], summary['avg_recall']))

    sorted_by = sorted(summary_per_query, key=lambda x: x['avg_precision'] + x['avg_recall'], reverse=True)
    f.write("\nAverage precision + average recall\n")
    for summary in sorted_by:
        f.write("%s \t %1.3f\n" % (summary['binning_label'], summary['avg_precision'] + summary['avg_recall']))
    f.close()


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics and figures for one or more binning files; output summary to screen and results per binning file to chosen directory",
                                     parents=[argparse_parents.PARSER_MULTI2])
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    args = parser.parse_args()
    binning_labels = []
    if args.labels:
        binning_labels = [x.strip() for x in args.labels.split(',')]
        if len(binning_labels) != len(args.bin_files):
            parser.error('number of labels does not match the number of binning files')
    summary_per_query = evaluate_all(args.gold_standard_file,
                                     args.fasta_file,
                                     args.bin_files,
                                     binning_labels,
                                     args.filter,
                                     args.genomes_file,
                                     args.keyword,
                                     args.output_dir)
    print_summary(convert_summary_to_tuples_of_strings(summary_per_query))
    plot_avg_precision_recall(summary_per_query, args.output_dir)
    plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, args.output_dir)
    compute_rankings(summary_per_query, args.output_dir)


if __name__ == "__main__":
    main()
