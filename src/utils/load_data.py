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

import os
import sys
import traceback
import logging
from collections import defaultdict
from src import binning_classes
from src import new_classes
from src.utils import labels as utils_labels

try:
    import load_ncbi_taxinfo
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import load_ncbi_taxinfo
    finally:
        sys.path.remove(os.path.dirname(__file__))


def load_unique_common(unique_common_file_path):
    if not unique_common_file_path:
        return None
    logging.getLogger('amber').info('Loading list of genomes to be removed...')
    genome_to_unique_common = {}
    with open(unique_common_file_path) as read_handler:
        for line in read_handler:
            genome_to_unique_common[line.split("\t")[0]] = line.split("\t")[1].strip('\n')
    logging.getLogger('amber').info('done')
    return genome_to_unique_common


def read_cami_file(read_handler, file_path, function_get_index, function_read_row, is_gs=False):
    header = {}
    column_name_to_index = {}
    got_column_indices = False
    reading_data = False

    for line in read_handler:
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            continue

        # parse header with column indices
        if line.startswith("@@"):
            for index, column_name in enumerate(line[2:].split('\t')):
                column_name_to_index[column_name] = index
            indices = function_get_index(column_name_to_index, is_gs)
            got_column_indices = True
            reading_data = False
            continue

        # parse header with metadata
        if line.startswith("@"):
            if reading_data:
                header = {}
            key, value = line[1:].split(':', 1)
            header[key.upper()] = value.strip()
            got_column_indices = False
            reading_data = False
            continue

        if not got_column_indices:
            logging.getLogger('amber').critical("Header line starting with @@ in file {} is missing or at wrong position.\n".format(file_path))
            exit(1)

        if 'SAMPLEID' not in header:
            logging.getLogger('amber').critical("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID.\n".format(file_path))
            exit(1)

        reading_data = True
        row_values = function_read_row(line, indices, is_gs)
        yield header['SAMPLEID'], row_values


def read_row(line, index_vars, is_gs):
    index_seq_id, index_bin_id, index_tax_id, index_length = index_vars
    line = line.rstrip('\n')
    row_data = line.split('\t')
    try:
        seq_id = row_data[index_seq_id]
    except:
        logging.getLogger('amber').critical("Value in column SEQUENCEID could not be read.")
        exit(1)

    if index_bin_id:
        try:
            bin_id = row_data[index_bin_id]
        except:
            logging.getLogger('amber').critical("Value in column BINID could not be read.")
            exit(1)
    else:
        bin_id = None

    if index_tax_id:
        try:
            tax_id = row_data[index_tax_id]
        except:
            logging.getLogger('amber').critical("Value in column TAXID could not be read.")
            exit(1)
    else:
        tax_id = None

    if is_gs and index_length is not None:
        try:
            length = int(row_data[index_length])
        except:
            logging.getLogger('amber').critical("Value in column LENGTH could not be read.")
            exit(1)
        return seq_id, bin_id, tax_id, length
    else:
        return seq_id, bin_id, tax_id, int(0)


def read_row_genome_coverage(line, index_vars, is_gs=False):
    index_genome_id, index_coverage, = index_vars
    line = line.rstrip('\n')
    row_data = line.split('\t')
    try:
        genome_id = row_data[index_genome_id]
    except:
        logging.getLogger('amber').critical("Value in column GENOMEID could not be read.")
        exit(1)

    if index_coverage:
        try:
            coverage = row_data[index_coverage]
        except:
            logging.getLogger('amber').critical("Value in column COVERAGE could not be read.")
            exit(1)
    return genome_id, float(coverage)


def get_column_indices(column_name_to_index, is_gs):
    if "SEQUENCEID" not in column_name_to_index:
        logging.getLogger('amber').critical("Column not found: {}".format("SEQUENCEID"))
        exit(1)
    if "BINID" not in column_name_to_index and "TAXID" not in column_name_to_index:
        logging.getLogger('amber').critical("Column not found: {}".format("BINID/TAXID"))
        exit(1)
    if is_gs and "LENGTH" not in column_name_to_index and "_LENGTH" not in column_name_to_index:
        logging.getLogger('amber').critical("Column not found: {}".format("LENGTH"))
        exit(1)
    index_seq_id = column_name_to_index["SEQUENCEID"]

    index_bin_id = None
    index_tax_id = None
    if "BINID" in column_name_to_index:
        index_bin_id = column_name_to_index["BINID"]
    if "TAXID" in column_name_to_index:
        index_tax_id = column_name_to_index["TAXID"]

    index_length = None
    if "LENGTH" in column_name_to_index:
        index_length = column_name_to_index["LENGTH"]
    elif "_LENGTH" in column_name_to_index:
        index_length = column_name_to_index["_LENGTH"]
    return index_seq_id, index_bin_id, index_tax_id, index_length


def get_column_indices_genome_coverage(column_name_to_index, is_gs=False):
    if "GENOMEID" not in column_name_to_index:
        logging.getLogger('amber').critical("Column not found: {}".format("GENOMEID"))
        exit(1)
    if "COVERAGE" not in column_name_to_index:
        logging.getLogger('amber').critical("Column not found: {}".format("COVERAGE"))
        exit(1)
    return column_name_to_index["GENOMEID"], column_name_to_index["COVERAGE"]


def initialize_query(sample_id_to_g_gold_standard, sample_id_to_t_gold_standard, is_gs, sample_id, options, label):
    g_query = binning_classes.GenomeQuery()
    t_query = binning_classes.TaxonomicQuery()
    g_query.options = options
    t_query.options = options
    t_query.label = label
    g_query.label = label
    g_sequence_ids = {}
    t_sequence_ids = {}
    sequence_id_to_length = None
    if is_gs:
        g_query.sequence_id_to_length = t_query.sequence_id_to_length = sequence_id_to_length = {}
        g_query.gold_standard = g_query
        t_query.gold_standard = t_query
    else:
        if sample_id_to_g_gold_standard and sample_id in sample_id_to_g_gold_standard:
            g_query.gold_standard = sample_id_to_g_gold_standard[sample_id]
            sequence_id_to_length = sample_id_to_g_gold_standard[sample_id].sequence_id_to_length
            g_sequence_ids = sample_id_to_g_gold_standard[sample_id].get_sequence_ids()
        if sample_id_to_t_gold_standard and sample_id in sample_id_to_t_gold_standard:
            t_query.gold_standard = sample_id_to_t_gold_standard[sample_id]
            sequence_id_to_length = sample_id_to_t_gold_standard[sample_id].sequence_id_to_length
            t_sequence_ids = sample_id_to_t_gold_standard[sample_id].get_sequence_ids()
    if sequence_id_to_length is None:
        logging.getLogger('amber').critical("Sample ID {} not found in the gold standard.".format(sample_id))
        exit(1)

    return g_query, t_query, g_sequence_ids, t_sequence_ids, sequence_id_to_length


def open_query(file_path_query, is_gs, sample_id_to_g_gold_standard, sample_id_to_t_gold_standard, options, label):
    sample_id_to_g_query = {}
    sample_id_to_t_query = {}
    sample_ids_list = []

    with open(file_path_query) as read_handler:
        try:
            sample_id_prev = 0
            for sample_id, (sequence_id, bin_id, tax_id, length) in read_cami_file(read_handler, file_path_query, get_column_indices, read_row, is_gs):

                if sample_id != sample_id_prev:
                    sample_ids_list.append(sample_id)
                    g_query, t_query, g_sequence_ids, t_sequence_ids, sequence_id_to_length = \
                        initialize_query(sample_id_to_g_gold_standard, sample_id_to_t_gold_standard, is_gs, sample_id, options, label)
                    sample_id_to_g_query[sample_id] = g_query
                    sample_id_to_t_query[sample_id] = t_query

                    if is_gs and not length:
                        logging.getLogger('amber').critical("Sequences length could not be determined. Please add column _LENGTH to gold standard.")
                        exit(1)

                sample_id_prev = sample_id

                if is_gs:
                    sequence_id_to_length[sequence_id] = length
                elif sequence_id not in sequence_id_to_length:
                    logging.getLogger('amber').warning("Ignoring sequence {} - length unknown (file {})".format(sequence_id, file_path_query))
                    continue

                if sequence_id_to_length[sequence_id] < options.min_length:
                    logging.getLogger('amber').warning("Ignoring sequence {} - shorter than {} bps: (file {})".format(sequence_id, options.min_length, file_path_query))
                    continue

                if bin_id:
                    if not is_gs and sequence_id not in g_sequence_ids:
                        logging.getLogger('amber').warning("Ignoring sequence {} - not found in the genome binning gold standard, sample {}: (file {})".format(sequence_id, sample_id, file_path_query))
                    else:
                        if bin_id not in g_query.get_bin_ids():
                            bin = binning_classes.GenomeBin(bin_id)
                            g_query.add_bin(bin)
                        else:
                            bin = g_query.get_bin_by_id(bin_id)
                        g_query.sequence_id_to_bin_id = (sequence_id, bin_id)
                        if is_gs:
                            bin.mapping_id = bin_id

                if tax_id:
                    if not binning_classes.TaxonomicQuery.tax_id_to_parent:
                        if is_gs:
                            continue
                        else:
                            logging.getLogger('amber').critical("Taxonomic binning cannot be assessed. Please provide an NCBI nodes file using option --ncbi_nodes_file.")
                            exit(1)

                    if tax_id not in binning_classes.TaxonomicQuery.tax_id_to_rank:
                        if binning_classes.TaxonomicQuery.tax_id_to_tax_id and tax_id in binning_classes.TaxonomicQuery.tax_id_to_tax_id:
                            tax_id = binning_classes.TaxonomicQuery.tax_id_to_tax_id[tax_id]
                        else:
                            logging.getLogger('amber').warning("Ignoring sequence {} - not a valid NCBI taxonomic ID: {} (file {})".format(sequence_id, tax_id, file_path_query))
                            continue

                    if not is_gs and sequence_id not in t_sequence_ids:
                        logging.getLogger('amber').warning("Ignoring sequence {} - not found in the taxonomic binning gold standard, sample {}: (file {})".format(sequence_id, sample_id, file_path_query))
                        continue

                    tax_id_path = load_ncbi_taxinfo.get_id_path(tax_id, binning_classes.TaxonomicQuery.tax_id_to_parent, binning_classes.TaxonomicQuery.tax_id_to_rank)

                    if not tax_id_path:
                        continue
                    for tax_id in tax_id_path:
                        if not tax_id: # tax_id may be empty
                            continue
                        if tax_id not in t_query.get_bin_ids():
                            bin = binning_classes.TaxonomicBin(tax_id)
                            bin.rank = binning_classes.TaxonomicQuery.tax_id_to_rank[tax_id]
                            t_query.add_bin(bin)
                        t_query.rank_to_sequence_id_to_bin_id = (binning_classes.TaxonomicQuery.tax_id_to_rank[tax_id], sequence_id, tax_id)
        except BaseException as e:
            traceback.print_exc()
            logging.getLogger('amber').critical("File {} is malformed. {}".format(file_path_query, e))
            exit(1)

    empty = []
    for sample_id, g_query in sample_id_to_g_query.items():
        if not g_query.bins:
            empty.append(sample_id)
    [sample_id_to_g_query.pop(sample_id) for sample_id in empty]

    empty = []
    for sample_id, t_query in sample_id_to_t_query.items():
        if not t_query.bins:
            empty.append(sample_id)
    [sample_id_to_t_query.pop(sample_id) for sample_id in empty]

    if not sample_id_to_g_query:
        sample_id_to_g_query = None
    if not sample_id_to_t_query:
        sample_id_to_t_query = None

    return sample_ids_list, sample_id_to_g_query, sample_id_to_t_query


def open_coverages(file_path):
    logging.getLogger('amber').info('Loading coverage file...')
    sample_id_to_coverages = {}
    with open(file_path) as read_handler:
        try:
            sample_id_prev = 0
            for sample_id, (genome_id, coverage) in read_cami_file(read_handler, file_path, get_column_indices_genome_coverage, read_row_genome_coverage, False):
                if sample_id != sample_id_prev:
                    genome_id_to_coverage = {}
                    sample_id_to_coverages[sample_id] = genome_id_to_coverage
                sample_id_prev = sample_id
                genome_id_to_coverage[genome_id] = coverage
        except BaseException as e:
            traceback.print_exc()
            logging.getLogger('amber').critical("File {} is malformed. {}".format(file_path, e))
            exit(1)
    logging.getLogger('amber').info('done')
    return sample_id_to_coverages


# def load_ncbi_info(ncbi_nodes_file, ncbi_names_file, ncbi_merged_file):
#     if ncbi_nodes_file:
#         logging.getLogger('amber').info('Loading NCBI files...')
#         binning_classes.TaxonomicQuery.tax_id_to_parent, binning_classes.TaxonomicQuery.tax_id_to_rank = \
#             load_ncbi_taxinfo.load_tax_info(ncbi_nodes_file)
#         if ncbi_names_file:
#             binning_classes.TaxonomicQuery.tax_id_to_name = load_ncbi_taxinfo.load_names(binning_classes.TaxonomicQuery.tax_id_to_rank, ncbi_names_file)
#         if ncbi_merged_file:
#             binning_classes.TaxonomicQuery.tax_id_to_tax_id = load_ncbi_taxinfo.load_merged(ncbi_merged_file)
#         logging.getLogger('amber').info('done')


def load_ncbi_info(ncbi_nodes_file, ncbi_names_file, ncbi_merged_file):
    if ncbi_nodes_file:
        logging.getLogger('amber').info('Loading NCBI files...')
        new_classes.TaxonomicQueryNEW.tax_id_to_parent, new_classes.TaxonomicQueryNEW.tax_id_to_rank = \
            load_ncbi_taxinfo.load_tax_info(ncbi_nodes_file)
        # if ncbi_names_file:
        #     binning_classes.TaxonomicQuery.tax_id_to_name = load_ncbi_taxinfo.load_names(binning_classes.TaxonomicQuery.tax_id_to_rank, ncbi_names_file)
        # if ncbi_merged_file:
        #     binning_classes.TaxonomicQuery.tax_id_to_tax_id = load_ncbi_taxinfo.load_merged(ncbi_merged_file)
        logging.getLogger('amber').info('done')


def create_genome_queries_from_taxonomic_queries(rank, sample_id_to_g_query, sample_id_to_queries_list):
    sample_id_to_genome_queries = defaultdict(list)
    for sample_id in sample_id_to_queries_list:
        for query in sample_id_to_queries_list[sample_id]:
            if isinstance(query, binning_classes.GenomeQuery):
                continue
            if sample_id not in sample_id_to_g_query:
                logging.getLogger('amber').warning("No genome binning gold standard available for sample " + sample_id)
                continue
            g_query = binning_classes.GenomeQuery()
            sample_id_to_genome_queries[sample_id].append(g_query)
            g_query.options = query.options
            g_query.label = query.label
            g_query.gold_standard = sample_id_to_g_query[sample_id]

            for sequence_id in query.rank_to_sequence_id_to_bin_id[rank]:
                bin_id = query.rank_to_sequence_id_to_bin_id[rank][sequence_id]
                if bin_id not in g_query.get_bin_ids():
                    bin = binning_classes.GenomeBin(bin_id)
                    g_query.add_bin(bin)
                g_query.sequence_id_to_bin_id = (sequence_id, bin_id)

    for sample_id in sample_id_to_genome_queries:
        sample_id_to_queries_list[sample_id].extend(sample_id_to_genome_queries[sample_id])


def get_gs_sample_id_to_num_genomes(sample_id_to_g_query, options):
    if not sample_id_to_g_query:
        return None
    sample_id_to_num_genomes = {}
    if options.genome_to_unique_common:
        for sample_id, g_query in sample_id_to_g_query.items():
            count = 0
            for bin in g_query.bins:
                if bin.id not in options.genome_to_unique_common or (options.filter_keyword and options.genome_to_unique_common[bin.id] != options.filter_keyword):
                    count += 1
            sample_id_to_num_genomes[sample_id] = count
    else:
        for sample_id, g_query in sample_id_to_g_query.items():
            sample_id_to_num_genomes[sample_id] = len(g_query.bins)
    return sample_id_to_num_genomes


def load_queries(gold_standard_file, query_files, options, labels):
    logging.getLogger('amber').info('Loading binnings...')
    sample_ids_list, sample_id_to_g_query, sample_id_to_t_query = \
        open_query(gold_standard_file,
                   True,
                   None, None,
                   options,
                   utils_labels.GS)
    sample_id_to_num_genomes = get_gs_sample_id_to_num_genomes(sample_id_to_g_query, options)

    sample_id_to_queries_list = defaultdict(list)
    for query_file, label in zip(query_files, labels):
        query = open_query(query_file,
                           False,
                           sample_id_to_g_query, sample_id_to_t_query,
                           options,
                           label)
        for sample_id_to_query in [query[1], query[2]]:
            if not sample_id_to_query:
                continue
            for sample_id in sample_id_to_query:
                sample_id_to_queries_list[sample_id].append(sample_id_to_query[sample_id])

    if options.rank_as_genome_binning:
        create_genome_queries_from_taxonomic_queries(options.rank_as_genome_binning, sample_id_to_g_query, sample_id_to_queries_list)

    sample_id_to_gs_list = defaultdict(list)
    for sample_id_to_query in [sample_id_to_g_query, sample_id_to_t_query]:
        if not sample_id_to_query:
            continue
        for sample_id in sample_id_to_query:
            sample_id_to_gs_list[sample_id].append(sample_id_to_query[sample_id])

    logging.getLogger('amber').info('done')
    return sample_id_to_gs_list, sample_ids_list, sample_id_to_num_genomes, sample_id_to_queries_list
