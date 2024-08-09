#!/usr/bin/env python3

# Copyright 2024 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from cami_amber import amber_html
from cami_amber import evaluate
from cami_amber import plots
from cami_amber import binning_classes
from cami_amber.utils import load_data
from cami_amber.utils import argparse_parents
from cami_amber.utils import labels as utils_labels
from version import __version__
from collections import defaultdict
import argparse
import errno
import logging
import os
import sys
import pandas as pd


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
    logging.getLogger('amber').info('Creating output directories')
    for sample_id in sample_id_to_queries_list:
        for query in sample_id_to_queries_list[sample_id]:
            make_sure_path_exists(os.path.join(output_dir, query.binning_type, query.label))


def get_labels(labels, bin_files):
    if labels:
        labels_list = [x.strip() for x in labels.split(',')]
        if len(set(labels_list)) != len(bin_files):
            logging.getLogger('amber').critical(
                'Number of different labels does not match the number of binning files. Please check parameter -l, --labels.')
            exit(1)
        return labels_list
    tool_id = []
    for bin_file in bin_files:
        tool_id.append(bin_file.split('/')[-1].split('.binning')[0])
    return tool_id


def get_num_threads(num_threads):
    try:
        return abs(int(num_threads)) if num_threads else None
    except ValueError:
        logging.getLogger('amber').error('Invalid value for number of threads')
        exit(1)


def save_metrics(sample_id_to_queries_list, df_summary, pd_bins, output_dir, stdout):
    logging.getLogger('amber').info('Saving computed metrics')
    df_summary.to_csv(os.path.join(output_dir, 'results.tsv'), sep='\t', index=False)
    pd_bins.to_csv(os.path.join(output_dir, 'bin_metrics.tsv'), index=False, sep='\t')
    if stdout:
        print(df_summary[[label for label in utils_labels.LABELS1 if label in df_summary.columns]].rename(columns=utils_labels.LABELS).to_string(index=False))
    for tool, pd_group in pd_bins[pd_bins['rank'] == 'NA'].groupby(utils_labels.TOOL):
        bins_columns = utils_labels.get_genome_bins_columns()
        table = pd_group[['sample_id'] + list(bins_columns.keys())].rename(columns=dict(bins_columns))
        table.to_csv(os.path.join(output_dir, 'genome', tool, 'metrics_per_bin.tsv'), sep='\t', index=False)
    for tool, pd_group in pd_bins[pd_bins['rank'] != 'NA'].groupby(utils_labels.TOOL):
        bins_columns = utils_labels.get_tax_bins_columns()
        if 'name' not in pd_bins.columns or pd_group['name'].isnull().any():
            del bins_columns['name']
        table = pd_group[['sample_id'] + list(bins_columns.keys())].rename(columns=dict(bins_columns))
        table.to_csv(os.path.join(output_dir, 'taxonomic', tool, 'metrics_per_bin.tsv'), sep='\t', index=False)

    pd_genomes_all = pd.DataFrame()
    for sample_id in sample_id_to_queries_list:
        pd_genomes_sample = pd.DataFrame()
        for query in sample_id_to_queries_list[sample_id]:
            if isinstance(query, binning_classes.GenomeQuery):
                query.recall_df_cami1[utils_labels.TOOL] = query.label
                pd_genomes_sample = pd.concat([pd_genomes_sample, query.recall_df_cami1], ignore_index=True, sort=False)
        pd_genomes_sample['sample_id'] = sample_id
        pd_genomes_all = pd.concat([pd_genomes_all, pd_genomes_sample], ignore_index=True, sort=False)
    if not pd_genomes_all.empty:
        pd_genomes_all.to_csv(os.path.join(output_dir, 'genome_metrics_cami1.tsv'), index=False, sep='\t')


def main(args=None):
    parser = argparse.ArgumentParser(description="AMBER: Assessment of Metagenome BinnERs",
                                     parents=[argparse_parents.PARSER_MULTI2], prog='AMBER')
    parser.add_argument('-p', '--filter', help=argparse_parents.HELP_FILTER)
    parser.add_argument('-n', '--min_length', help="Minimum length of sequences", type=int, required=False)
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    parser.add_argument('--stdout', help="Print summary to stdout", action='store_true')
    parser.add_argument('-d', '--desc', help="Description for HTML page", required=False)
    parser.add_argument('--colors', help="Color indices", required=False)
    parser.add_argument('--silent', help='Silent mode', action='store_true')
    parser.add_argument('-t', '--threads', help='Number of threads (default: number of CPUs)', required=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    group_g = parser.add_argument_group('genome binning-specific arguments')
    group_g.add_argument('-x', '--min_completeness', help=argparse_parents.HELP_THRESHOLDS_COMPLETENESS, required=False)
    group_g.add_argument('-y', '--max_contamination', help=argparse_parents.HELP_THRESHOLDS_CONTAMINATION,
                         required=False)
    group_g.add_argument('-r', '--remove_genomes', help=argparse_parents.HELP_GENOMES_FILE, required=False)
    group_g.add_argument('-k', '--keyword', help=argparse_parents.HELP_KEYWORD, required=False)
    group_g.add_argument('--genome_coverage', help='genome coverages', required=False)

    group_t = parser.add_argument_group('taxonomic binning-specific arguments')
    group_t.add_argument('--ncbi_dir', help="Directory containing the NCBI taxonomy database dump files nodes.dmp, merged.dmp, and names.dmp", required=False)
    # group_t.add_argument('--rank_as_genome_binning',
    #                      help="Assess taxonomic binning at a rank also as genome binning. Valid ranks: superkingdom, phylum, class, order, family, genus, species, strain",
    #                      required=False)

    args = parser.parse_args(args)
    output_dir = os.path.abspath(args.output_dir)
    logger = get_logger(output_dir, args.silent)

    labels = get_labels(args.labels, args.bin_files)

    genome_to_unique_common = load_data.load_unique_common(args.remove_genomes, args.keyword)

    options = binning_classes.Options(filter_tail_percentage=args.filter,
                                      genome_to_unique_common=genome_to_unique_common,
                                      filter_keyword=args.keyword,
                                      min_length=args.min_length,
                                      rank_as_genome_binning=None, #args.rank_as_genome_binning,
                                      output_dir=output_dir,
                                      min_completeness=args.min_completeness,
                                      max_contamination=args.max_contamination,
                                      ncbi_dir=args.ncbi_dir,
                                      skip_gs=False)
    options_gs = binning_classes.Options(filter_tail_percentage=.0,
                                         genome_to_unique_common=genome_to_unique_common,
                                         filter_keyword=args.keyword,
                                         min_length=args.min_length,
                                         rank_as_genome_binning=None, #args.rank_as_genome_binning,
                                         output_dir=output_dir,
                                         ncbi_dir=args.ncbi_dir,
                                         skip_gs=False)

    num_threads = get_num_threads(args.threads)
    sample_id_to_g_queries_list, sample_id_to_t_queries_list, sample_ids_list = load_data.load_queries_mthreaded(
        args.gold_standard_file, args.bin_files, labels, options, options_gs, num_threads)

    coverages_pd = load_data.open_coverages(args.genome_coverage)

    sample_id_to_queries_list = defaultdict(list)
    for sample_id in sample_id_to_g_queries_list | sample_id_to_t_queries_list:
        sample_id_to_queries_list[sample_id] += sample_id_to_g_queries_list[sample_id]
        sample_id_to_queries_list[sample_id] += sample_id_to_t_queries_list[sample_id]

    create_output_directories(output_dir, sample_id_to_queries_list)

    df_summary, pd_bins = evaluate.evaluate_samples_queries(sample_id_to_g_queries_list, sample_id_to_t_queries_list, num_threads)

    save_metrics(sample_id_to_queries_list, df_summary, pd_bins, output_dir, args.stdout)

    plots.plot_genome_binning(args.colors,
                              sample_id_to_g_queries_list,
                              df_summary,
                              pd_bins[pd_bins['rank'] == 'NA'],
                              labels,
                              coverages_pd,
                              output_dir)
    plots.plot_taxonomic_binning(args.colors, df_summary, pd_bins, labels, output_dir)

    amber_html.create_html(df_summary,
                           pd_bins,
                           [utils_labels.GS] + labels,
                           sample_ids_list,
                           options,
                           args.desc)
    logger.info('AMBER finished successfully. All results have been saved to {}'.format(output_dir))


if __name__ == "__main__":
    main()
