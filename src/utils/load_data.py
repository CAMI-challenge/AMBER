# Copyright 2020 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
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

import pandas as pd
import sys
import logging
import traceback
import os
from collections import defaultdict
from collections import OrderedDict
from src import binning_classes
from src.utils import labels as utils_labels

try:
    import load_ncbi_taxinfo
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import load_ncbi_taxinfo
    finally:
        sys.path.remove(os.path.dirname(__file__))


def open_coverages(file_path):
    logging.getLogger('amber').info('Loading coverage file')
    sample_id_to_coverages = {}
    with open(file_path) as read_handler:
        try:
            sample_id_prev = 0
            # for sample_id, (genome_id, coverage) in read_metadata(read_handler, file_path, get_column_indices_genome_coverage, read_row_genome_coverage, False):
            for sample_id, (genome_id, coverage) in read_metadata(file_path):
                if sample_id != sample_id_prev:
                    genome_id_to_coverage = {}
                    sample_id_to_coverages[sample_id] = genome_id_to_coverage
                sample_id_prev = sample_id
                genome_id_to_coverage[genome_id] = coverage
        except BaseException as e:
            traceback.print_exc()
            logging.getLogger('amber').critical("File {} is malformed. {}".format(file_path, e))
            exit(1)
    return sample_id_to_coverages


def load_unique_common(unique_common_file_path, args_keyword):
    if not unique_common_file_path:
        return None
    logging.getLogger('amber').info('Loading list of genomes to be removed')
    with open(unique_common_file_path) as f:
        line = f.readline()
        line_split = line.split('\t')
        f.seek(0)
        if len(line_split) == 1:
            return [line.strip('\n') for line in f]
        elif len(line_split) > 1:
            genome_to_unique_common = {}
            for line in f:
                line_split = line.split('\t')
                genome_to_unique_common[line_split[0]] = line_split[1].strip('\n')
            if args_keyword:
                return [genome_id for genome_id, keyword in genome_to_unique_common.items() if keyword == args_keyword]
            else:
                return list(genome_to_unique_common.keys())
    return None


def load_ncbi_info(ncbi_nodes_file, ncbi_names_file, ncbi_merged_file):
    if ncbi_nodes_file:
        logging.getLogger('amber').info('Loading NCBI files')
        binning_classes.TaxonomicQuery.tax_id_to_parent, binning_classes.TaxonomicQuery.tax_id_to_rank = load_ncbi_taxinfo.load_tax_info(ncbi_nodes_file)
        if ncbi_names_file:
            binning_classes.TaxonomicQuery.tax_id_to_name = load_ncbi_taxinfo.load_names(binning_classes.TaxonomicQuery.tax_id_to_rank, ncbi_names_file)
        if ncbi_merged_file:
            binning_classes.TaxonomicQuery.tax_id_to_tax_id = load_ncbi_taxinfo.load_merged(ncbi_merged_file)


def read_metadata(file_path_query):
    header = {}
    index_to_column_name = {}
    got_column_indices = False

    samples_metadata = []

    with open(file_path_query, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()

            if len(line) == 0 or line.startswith("#"):
                if got_column_indices and not reading_data:
                    data_start = i + 1
                continue

            # parse header with column indices
            if line.startswith('@@'):
                for index, column_name in enumerate(line[2:].split('\t')):
                    index_to_column_name[index] = column_name
                got_column_indices = True
                reading_data = False
                data_start = i + 1

            elif line.startswith('@'):
                # parse header with metadata
                if got_column_indices:
                    samples_metadata.append((data_start, data_end, header, index_to_column_name))
                    header = {}
                    index_to_column_name = {}
                key, value = line[1:].split(':', 1)
                header[key.upper()] = value.strip()
                got_column_indices = False
                reading_data = False

            else:
                reading_data = True
                data_end = i
    samples_metadata.append((data_start, data_end, header, index_to_column_name))
    return samples_metadata


def load_binnings(samples_metadata, file_path_query):
    sample_id_to_query_df = OrderedDict()
    for metadata in samples_metadata:
        nrows = metadata[1] - metadata[0] + 1

        df = pd.read_csv(file_path_query, sep='\t', comment='#', skiprows=metadata[0], nrows=nrows, header=None)
        df.rename(columns=metadata[3], inplace=True)
        df.rename(columns={'_LENGTH': 'LENGTH'}, inplace=True)
        sample_id_to_query_df[metadata[2]['SAMPLEID']] = df
    return sample_id_to_query_df


def open_query(file_path_query):
    try:
        samples_metadata = read_metadata(file_path_query)
        sample_id_to_query_df = load_binnings(samples_metadata, file_path_query)
    except BaseException as e:
        traceback.print_exc()
        logging.getLogger('amber').critical("File {} is malformed. {}".format(file_path_query, e))
        exit(1)
    return sample_id_to_query_df


def get_sample_id_to_num_genomes(sample_id_to_queries_list):
    sample_id_to_num_genomes = {}
    for sample_id in sample_id_to_queries_list:
        if isinstance(sample_id_to_queries_list[sample_id][0], binning_classes.GenomeQuery):
            gs_df = sample_id_to_queries_list[sample_id][0].df
            sample_id_to_num_genomes[sample_id] = gs_df['BINID'].nunique()
        else:
            sample_id_to_num_genomes[sample_id] = 0
    return sample_id_to_num_genomes


def get_rank_to_df(query_df, is_gs=False):
    query_df.sort_values(by=['TAXID'], inplace=True)
    rank_to_sequence_id_to_bin_id = defaultdict(dict)
    tax_id = None
    for row in zip(query_df.SEQUENCEID, query_df.TAXID):
        if row[1] != tax_id:
            tax_id = row[1]
            tax_id_path = load_ncbi_taxinfo.get_id_path(tax_id, binning_classes.TaxonomicQuery.tax_id_to_parent,
                                                        binning_classes.TaxonomicQuery.tax_id_to_rank,
                                                        binning_classes.TaxonomicQuery.tax_id_to_tax_id)
            if not tax_id_path:
                continue
        for tax_id2 in tax_id_path:
            if not tax_id2:  # tax_id may be empty
                continue
            rank = binning_classes.TaxonomicQuery.tax_id_to_rank[tax_id2]
            rank_to_sequence_id_to_bin_id[rank][row[0]] = tax_id2
    rank_to_df = dict()
    if is_gs:
        if 'LENGTH' not in query_df.columns:
            logging.getLogger('amber').critical("Sequences length could not be determined. Please add column LENGTH to gold standard.")
            exit(1)
        query_df = query_df[['SEQUENCEID', 'LENGTH']].set_index('SEQUENCEID')
    for rank in rank_to_sequence_id_to_bin_id:
        rank_to_df[rank] = pd.DataFrame.from_dict(rank_to_sequence_id_to_bin_id[rank], orient='index',
                                                  columns=['TAXID']).reset_index().rename(columns={'index': 'SEQUENCEID'}).set_index('SEQUENCEID')
        if is_gs:
            rank_to_df[rank] = rank_to_df[rank].join(query_df, on='SEQUENCEID', how='left', sort=False)
    return rank_to_df


def filter_by_length(sample_id_to_df, min_length):
    sample_id_to_df_filtered = OrderedDict()
    for sample_id in sample_id_to_df:
        df = sample_id_to_df[sample_id]
        sample_id_to_df_filtered[sample_id] = df[df['LENGTH'] >= min_length]
    return sample_id_to_df_filtered


def load_queries(gold_standard_file, query_files, labels, options, options_gs):
    logging.getLogger('amber').info('Loading {}'.format(utils_labels.GS))
    sample_id_to_gs_df = open_query(gold_standard_file)
    if options.min_length:
        sample_id_to_gs_df_filtered = filter_by_length(sample_id_to_gs_df, options.min_length)
    sample_id_to_gs_rank_to_df = defaultdict()
    sample_id_to_queries_list = defaultdict(list)
    sample_id_to_g_gs = defaultdict(list)
    sample_id_to_t_gs = defaultdict(list)

    for sample_id in sample_id_to_gs_df:
        gs_df = sample_id_to_gs_df_filtered[sample_id] if options.min_length else sample_id_to_gs_df[sample_id]
        if gs_df.empty:
            logging.getLogger('amber').critical("Gold standard for sample {} is empty. Exiting.".format(sample_id))
            exit(1)
        if 'BINID' in gs_df.columns:
            g_query_gs = binning_classes.GenomeQuery(gs_df, utils_labels.GS, sample_id, options_gs)
            g_query_gs.gold_standard = g_query_gs
            g_query_gs.gold_standard_df = gs_df
            sample_id_to_queries_list[sample_id].append(g_query_gs)
            sample_id_to_g_gs[sample_id] = g_query_gs
        if 'TAXID' in gs_df.columns and binning_classes.TaxonomicQuery.tax_id_to_rank:
            gs_rank_to_df = get_rank_to_df(gs_df, is_gs=True)
            sample_id_to_gs_rank_to_df[sample_id] = gs_rank_to_df
            t_query_gs = binning_classes.TaxonomicQuery(gs_rank_to_df, utils_labels.GS, sample_id, options_gs)
            t_query_gs.gold_standard = t_query_gs
            t_query_gs.gold_standard_df = gs_rank_to_df
            sample_id_to_queries_list[sample_id].append(t_query_gs)
            sample_id_to_t_gs[sample_id] = t_query_gs

    for query_file, label in zip(query_files, labels):
        logging.getLogger('amber').info('Loading {}'.format(label))
        sample_id_to_query_df = open_query(query_file)
        for sample_id in sample_id_to_query_df:
            if sample_id not in sample_id_to_gs_df:
                logging.getLogger('amber').critical("Sample ID {} in {} not found in the gold standard.".format(sample_id, label))
                exit(1)

            gs_df = sample_id_to_gs_df[sample_id]
            query_df = sample_id_to_query_df[sample_id]
            condition = query_df['SEQUENCEID'].isin(gs_df['SEQUENCEID'])
            if ~condition.all():
                logging.getLogger('amber').warning("{} sequences in {} not found in the gold standard.".format(query_df[~condition]['SEQUENCEID'].nunique(), label))
                query_df = query_df[condition]

            if options.min_length:
                gs_df = sample_id_to_gs_df_filtered[sample_id]
                condition = query_df['SEQUENCEID'].isin(gs_df['SEQUENCEID'])
                if ~condition.all():
                    query_df = query_df[condition]

            if 'BINID' in query_df.columns:
                g_query = binning_classes.GenomeQuery(query_df, label, sample_id, options)
                g_query.gold_standard = sample_id_to_g_gs[sample_id]
                g_query.gold_standard_df = gs_df
                sample_id_to_queries_list[sample_id].append(g_query)
                options.only_taxonomic_queries = options_gs.only_taxonomic_queries = False
            if 'TAXID' in query_df.columns and binning_classes.TaxonomicQuery.tax_id_to_rank:
                t_query = binning_classes.TaxonomicQuery(get_rank_to_df(query_df), label, sample_id, options)
                t_query.gold_standard = sample_id_to_t_gs[sample_id]
                t_query.gold_standard_df = sample_id_to_gs_rank_to_df[sample_id]
                sample_id_to_queries_list[sample_id].append(t_query)
                options.only_taxonomic_queries = options_gs.only_genome_queries = False

    return sample_id_to_queries_list, list(sample_id_to_gs_df.keys())

