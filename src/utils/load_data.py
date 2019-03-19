#!/usr/bin/env python

import os
import sys
import traceback
import logging
from collections import defaultdict
from src import binning_classes

try:
    import load_ncbi_taxinfo
    import add_length_column
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import load_ncbi_taxinfo
        import add_length_column
    finally:
        sys.path.remove(os.path.dirname(__file__))


def load_unique_common(unique_common_file_path):
    genome_to_unique_common = {}
    with open(unique_common_file_path) as read_handler:
        for line in read_handler:
            genome_to_unique_common[line.split("\t")[0]] = line.split("\t")[1].strip('\n')
    return genome_to_unique_common


def read_binning_file(read_handler, file_path, is_gs):
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
            index_seq_id, index_bin_id, index_tax_id, index_length = get_column_indices(column_name_to_index)
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
            raise RuntimeError

        if 'SAMPLEID' not in header:
            logging.getLogger('amber').critical("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID.\n".format(file_path))
            raise RuntimeError

        reading_data = True
        sequence_id, bin_id, tax_id, length = read_row(line, index_seq_id, index_bin_id, index_tax_id, index_length, is_gs)
        yield header['SAMPLEID'], sequence_id, bin_id, tax_id, length


def get_genome_mapping_without_lenghts(mapping_file, remove_genomes_file=None, keyword=None):
    gold_standard = binning_classes.Query()
    gold_standard.bin_id_to_list_of_sequence_ids = {}
    gold_standard.sequence_id_to_genome_id = {}

    filtering_genomes_to_keyword = {}
    if remove_genomes_file:
        filtering_genomes_to_keyword = load_unique_common(remove_genomes_file)

    with open(mapping_file, 'r') as read_handler:
        try:
            for sample_id, sequence_id, genome_id, length in read_binning_file(read_handler, mapping_file):
                if genome_id in filtering_genomes_to_keyword and (not keyword or filtering_genomes_to_keyword[genome_id] == keyword):
                    continue
                gold_standard.sequence_id_to_genome_id[sequence_id] = genome_id
                if genome_id not in gold_standard.bin_id_to_list_of_sequence_ids:
                    gold_standard.bin_id_to_list_of_sequence_ids[genome_id] = []
                gold_standard.bin_id_to_list_of_sequence_ids[genome_id].append(sequence_id)
        except:
            exit("Error. File {} is malformed.".format(mapping_file))

    if len(gold_standard.bin_id_to_list_of_sequence_ids) == 0:
        exit('All bins of the gold standard have been filtered out due to option --remove_genomes.')

    return gold_standard


def read_header(input_stream):
    """
    @Version:0.9.1
    @SampleID:RH_S1

    @@SEQUENCEID    BINID   TAXID
    """
    header = {}
    column_names = {}
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip('\n')
        if line.startswith("@@"):
            for index, column_name in enumerate(line[2:].split('\t')):
                column_names[column_name] = index
            return header, column_names
        if line.startswith("@"):
            key, value = line[1:].split(':', 1)
            header[key] = value.strip()


def read_row(line, index_seq_id, index_bin_id, index_tax_id, index_length, is_gs):
    line = line.rstrip('\n')
    row_data = line.split('\t')
    try:
        seq_id = row_data[index_seq_id]
    except:
        print("Value in column SEQUENCEID could not be read.")
        raise

    if index_bin_id:
        try:
            bin_id = row_data[index_bin_id]
        except:
            print("Value in column BINID could not be read.")
            raise
    else:
        bin_id = None

    if index_tax_id:
        try:
            tax_id = row_data[index_tax_id]
        except:
            print("Value in column TAXID could not be read.")
            raise
    else:
        tax_id = None

    if is_gs and index_length is not None:
        try:
            length = int(row_data[index_length])
        except:
            print("Value in column _LENGTH could not be read. Please provide a value or remove column altogether (and provide a FASTA or FASTQ file instead - see README).")
            raise
        return seq_id, bin_id, tax_id, length
    else:
        return seq_id, bin_id, tax_id, int(0)


def is_length_column_available(input_stream):
    header, column_names = read_header(input_stream)
    input_stream.seek(0)
    return "_LENGTH" in column_names


def get_column_indices(column_name_to_index):
    if "SEQUENCEID" not in column_name_to_index:
        raise RuntimeError("Column not found: {}".format("SEQUENCEID"))
    if "BINID" not in column_name_to_index and "TAXID" not in column_name_to_index:
        raise RuntimeError("Column not found: {}".format("BINID/TAXID"))
    index_seq_id = column_name_to_index["SEQUENCEID"]

    index_bin_id = None
    index_tax_id = None
    if "BINID" in column_name_to_index:
        index_bin_id = column_name_to_index["BINID"]
    if "TAXID" in column_name_to_index:
        index_tax_id = column_name_to_index["TAXID"]

    index_length = None
    if "_LENGTH" in column_name_to_index:
        index_length = column_name_to_index["_LENGTH"]
    return index_seq_id, index_bin_id, index_tax_id, index_length


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
        exit("Error. Sample ID {} not found in the gold standard.".format(sample_id))

    return g_query, t_query, g_sequence_ids, t_sequence_ids, sequence_id_to_length


def open_query(file_path_query, is_gs, fastx_file, sample_id_to_g_gold_standard, sample_id_to_t_gold_standard, options, label):
    sample_id_to_g_query = {}
    sample_id_to_t_query = {}
    sample_ids_list = []

    with open(file_path_query) as read_handler:
        try:
            sample_id_prev = 0
            for sample_id, sequence_id, bin_id, tax_id, length in read_binning_file(read_handler, file_path_query, is_gs):

                if sample_id != sample_id_prev:
                    sample_ids_list.append(sample_id)
                    g_query, t_query, g_sequence_ids, t_sequence_ids, sequence_id_to_length = \
                        initialize_query(sample_id_to_g_gold_standard, sample_id_to_t_gold_standard, is_gs, sample_id, options, label)
                    sample_id_to_g_query[sample_id] = g_query
                    sample_id_to_t_query[sample_id] = t_query

                    # TODO: re-enable check
                    # if is_gs and not binning_classes.Bin.sequence_id_to_length and is_length_column_available(read_handler):
                    #     if not fastx_file:
                    #         exit("Sequences length could not be determined. Please provide a FASTA or FASTQ file using option -f or add column _LENGTH to gold standard.")
                    #     binning_classes.Bin.sequence_id_to_length = add_length_column.read_lengths_from_fastx_file(fastx_file)

                sample_id_prev = sample_id

                if is_gs:
                    sequence_id_to_length[sequence_id] = length
                elif sequence_id not in sequence_id_to_length:
                    print("Ignoring sequence {} - length unknown (file {})".format(sequence_id, file_path_query), file=sys.stderr)
                    continue

                if bin_id:
                    if sequence_id_to_length[sequence_id] >= options.min_length:
                        if not is_gs and sequence_id not in g_sequence_ids:
                            print("Ignoring sequence {} - not found in the genome binning gold standard, sample {}: (file {})".format(sequence_id, sample_id, file_path_query), file=sys.stderr)
                        else:
                            if bin_id not in g_query.get_bin_ids():
                                bin = binning_classes.GenomeBin(bin_id)
                                g_query.add_bin(bin)
                            else:
                                bin = g_query.get_bin_by_id(bin_id)
                            g_query.sequence_id_to_bin_id = (sequence_id, bin_id)
                            bin.add_sequence_id(sequence_id, sequence_id_to_length[sequence_id])
                            if is_gs:
                                bin.mapping_id = bin_id
                    else:
                        print("Ignoring sequence {} - shorter than {} bps: (file {})".format(sequence_id, options.min_length, file_path_query), file=sys.stderr)

                if tax_id:
                    if not binning_classes.TaxonomicQuery.tax_id_to_parent:
                        if is_gs:
                            continue
                        else:
                            exit("Taxonomic binning cannot be assessed. Please provide an NCBI nodes file using option --ncbi_nodes_file.")

                    if tax_id not in binning_classes.TaxonomicQuery.tax_id_to_rank:
                        print("Ignoring sequence {} - not a valid NCBI taxonomic ID: {} (file {})".format(sequence_id, tax_id, file_path_query), file=sys.stderr)
                        continue

                    if not is_gs and sequence_id not in t_sequence_ids:
                        print("Ignoring sequence {} - not found in the taxonomic binning gold standard, sample {}: (file {})".format(sequence_id, sample_id, file_path_query), file=sys.stderr)
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
                        else:
                            bin = t_query.get_bin_by_id(tax_id)
                        t_query.rank_to_sequence_id_to_bin_id = (binning_classes.TaxonomicQuery.tax_id_to_rank[tax_id], sequence_id, tax_id)
                        bin.add_sequence_id(sequence_id, sequence_id_to_length[sequence_id])
        except BaseException as e:
            traceback.print_exc()
            exit("Error. File {} is malformed. {}".format(file_path_query, e))

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


def load_ncbi_info(ncbi_nodes_file, ncbi_names_file):
    if ncbi_nodes_file:
        binning_classes.TaxonomicQuery.tax_id_to_parent, binning_classes.TaxonomicQuery.tax_id_to_rank = \
            load_ncbi_taxinfo.load_tax_info(ncbi_nodes_file)
        if ncbi_names_file:
            binning_classes.TaxonomicQuery.tax_id_to_name = load_ncbi_taxinfo.load_names(binning_classes.TaxonomicQuery.tax_id_to_rank, ncbi_names_file)


def load_queries(gold_standard_file, fastx_file, query_files, options, labels):
    gold_standard = open_query(gold_standard_file,
                               True,
                               fastx_file,
                               None, None,
                               options,
                               None)
    sample_id_to_num_genomes = None
    if gold_standard[1]:
        sample_id_to_num_genomes = {}
        for sample_id, g_query in gold_standard[1].items():
            sample_id_to_num_genomes[sample_id] = len(g_query.bins)

    sample_id_to_queries_list = defaultdict(list)
    for query_file, label in zip(query_files, labels):
        query = open_query(query_file,
                           False,
                           None,
                           gold_standard[1], gold_standard[2],
                           options,
                           label)
        for sample_id_to_query in [query[1], query[2]]:
            if not sample_id_to_query:
                continue
            for sample_id in sample_id_to_query:
                sample_id_to_queries_list[sample_id].append(sample_id_to_query[sample_id])

    # TODO if there is a g_query (t_query), there must be a g_gold_standard (t_gold_standard)

    return [sample_id for sample_id in gold_standard[0]], sample_id_to_num_genomes, sample_id_to_queries_list
