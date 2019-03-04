#!/usr/bin/env python

import os
import sys
import numpy as np
import traceback
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


def load_tsv_table(stream):
    data = []
    next(stream)
    for line in stream:
        line = line.strip()
        if len(line) == 0 or line.startswith("@"):
            continue
        row_data = line.split('\t')

        mapped_genome = row_data[0]
        real_size = int(float(row_data[5]))
        predicted_size = int(float(row_data[3]))
        correctly_predicted = int(float(row_data[4]))

        if row_data[1] != "NA" and predicted_size > 0:
            precision = float(row_data[1])
        else:
            precision = np.nan
        data.append({'mapped_genome': mapped_genome, 'precision': precision, 'recall': float(row_data[2]),
                     'predicted_size': predicted_size, 'correctly_predicted': correctly_predicted, 'true_size': real_size})
    return data


def get_genome_mapping_without_lenghts(mapping_file, remove_genomes_file=None, keyword=None):
    gold_standard = binning_classes.Query()
    gold_standard.bin_id_to_list_of_sequence_ids = {}
    gold_standard.sequence_id_to_genome_id = {}

    filtering_genomes_to_keyword = {}
    if remove_genomes_file:
        filtering_genomes_to_keyword = load_unique_common(remove_genomes_file)

    with open(mapping_file, 'r') as read_handler:
        try:
            for sequence_id, genome_id, length in read_binning_file(read_handler):
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


def read_rows(input_stream, index_seq_id, index_bin_id, index_tax_id, index_length, is_gs):
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
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
            yield seq_id, bin_id, tax_id, length
        else:
            yield seq_id, bin_id, tax_id, int(0)


def is_length_column_available(input_stream):
    header, column_names = read_header(input_stream)
    input_stream.seek(0)
    return "_LENGTH" in column_names


def get_column_indices(input_stream):
    """

    :param input_stream: 
    :return: 
    """
    header, column_names = read_header(input_stream)
    if "SEQUENCEID" not in column_names:
        raise RuntimeError("Column not found: {}".format("SEQUENCEID"))
    if "BINID" not in column_names and "TAXID" not in column_names:
        raise RuntimeError("Column not found: {}".format("BINID/TAXID"))
    index_seq_id = column_names["SEQUENCEID"]

    index_bin_id = None
    index_tax_id = None
    if "BINID" in column_names:
        index_bin_id = column_names["BINID"]
    if "TAXID" in column_names:
        index_tax_id = column_names["TAXID"]

    index_length = None
    if "_LENGTH" in column_names:
        index_length = column_names["_LENGTH"]
    return index_seq_id, index_bin_id, index_tax_id, index_length


def read_binning_file(input_stream, is_gs):
    index_seq_id, index_bin_id, index_tax_id, index_length = get_column_indices(input_stream)
    return read_rows(input_stream, index_seq_id, index_bin_id, index_tax_id, index_length, is_gs)


def open_query(file_path_query, is_gs, fastx_file, g_gold_standard, t_gold_standard, options):
    g_query = binning_classes.GenomeQuery()
    t_query = binning_classes.TaxonomicQuery()
    if is_gs:
        g_query.gold_standard = g_query
        t_query.gold_standard = t_query
        min_length = options.min_length
    else:
        g_query.gold_standard = g_gold_standard
        t_query.gold_standard = t_gold_standard
        min_length = 0

    with open(file_path_query) as read_handler:
        if is_gs and not binning_classes.Bin.sequence_id_to_length:
            if not is_length_column_available(read_handler):
                if not fastx_file:
                    exit("Sequences length could not be determined. Please provide a FASTA or FASTQ file using option -f or add column _LENGTH to gold standard.")
                binning_classes.Bin.sequence_id_to_length = add_length_column.read_lengths_from_fastx_file(fastx_file)

        g_sequence_ids = {}
        t_sequence_ids = {}
        if g_gold_standard:
            g_sequence_ids = g_gold_standard.get_sequence_ids()
        if t_gold_standard:
            t_sequence_ids = t_gold_standard.get_sequence_ids()

        try:
            for sequence_id, bin_id, tax_id, length in read_binning_file(read_handler, is_gs):
                if is_gs:
                    binning_classes.Bin.sequence_id_to_length[sequence_id] = length
                elif sequence_id not in binning_classes.Bin.sequence_id_to_length:
                    print("Ignoring sequence {} - length unknown (file {})".format(sequence_id, file_path_query), file=sys.stderr)
                    continue

                if bin_id:
                    if binning_classes.Bin.sequence_id_to_length[sequence_id] >= min_length:
                        if not is_gs and sequence_id not in g_sequence_ids:
                            print("Ignoring sequence {} - not found in the genome binning gold standard: (file {})".format(sequence_id, file_path_query), file=sys.stderr)
                        else:
                            if bin_id not in g_query.get_bin_ids():
                                bin = binning_classes.GenomeBin(bin_id)
                                g_query.add_bin(bin)
                            else:
                                bin = g_query.get_bin_by_id(bin_id)
                            g_query.sequence_id_to_bin_id = (sequence_id, bin_id)
                            bin.add_sequence_id(sequence_id)
                            if is_gs:
                                bin.mapping_id = bin_id

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
                        print("Ignoring sequence {} - not found in the taxonomic binning gold standard: (file {})".format(sequence_id, file_path_query), file=sys.stderr)
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
                        bin.add_sequence_id(sequence_id)
        except BaseException as e:
            traceback.print_exc()
            exit("Error. File {} is malformed. {}".format(file_path_query, e))

    if g_query.bins:
        g_query.options = options
    else:
        g_query = None

    if t_query.bins:
        t_query.options = options
    else:
        t_query = None

    return g_query, t_query


def load_ncbi_info(ncbi_nodes_file, ncbi_names_file):
    if ncbi_nodes_file:
        binning_classes.TaxonomicQuery.tax_id_to_parent, binning_classes.TaxonomicQuery.tax_id_to_rank = \
            load_ncbi_taxinfo.load_tax_info(ncbi_nodes_file)
        if ncbi_names_file:
            binning_classes.TaxonomicQuery.tax_id_to_name = load_ncbi_taxinfo.load_names(binning_classes.TaxonomicQuery.tax_id_to_rank, ncbi_names_file)


def load_queries(gold_standard_file, fastx_file, query_files, options, labels):
    g_gold_standard, t_gold_standard = open_query(gold_standard_file,
                                                  True,
                                                  fastx_file,
                                                  None, None,
                                                  options)
    queries_list = []
    for query_file, label in zip(query_files, labels):
        g_query, t_query = open_query(query_file,
                                      False,
                                      None,
                                      g_gold_standard, t_gold_standard,
                                      options)
        for query in [g_query, t_query]:
            if not query:
                continue
            query.label = label
            queries_list.append(query)

    # TODO if there is a g_query (t_query), there must be a g_gold_standard (t_gold_standard)

    return queries_list
