#!/usr/bin/env python

import os
import gzip
import mimetypes
import sys
from Bio import SeqIO
import numpy as np
import traceback
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src import binning_classes

try:
    import load_ncbi_taxinfo
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import load_ncbi_taxinfo
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
                     'predicted_size': predicted_size, 'correctly_predicted': correctly_predicted, 'real_size': real_size})
    return data


def read_lengths_from_fastx_file(fastx_file):
    """

    @param fastx_file: file path
    @type fastx_file: str
    @rtype: dict[str, int]
    """
    file_type = mimetypes.guess_type(fastx_file)[1]
    if file_type == 'gzip':
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError("Unknown type of file: '{}".format(fastx_file))

    length = {}
    if os.path.getsize(fastx_file) == 0:
        return length

    file_format = None
    line = f.readline()
    if line.startswith('@'):
        file_format = "fastq"
    elif line.startswith(">"):
        file_format = "fasta"
    f.seek(0)

    if not file_format:
        raise RuntimeError("Invalid sequence file: '{}".format(fastx_file))

    for seq_record in SeqIO.parse(f, file_format):
        length[seq_record.id] = len(seq_record.seq)

    f.close()

    return length


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


def open_query(file_path_query, is_gs, fastx_file, tax_id_to_parent, tax_id_to_rank, gold_standard, min_length=0):
    g_query = binning_classes.GenomeQuery()
    t_query = binning_classes.TaxonomicQuery()
    t_query.tax_id_to_rank = tax_id_to_rank
    is_length_column_av = False

    with open(file_path_query) as read_handler:
        if is_gs and not binning_classes.Bin.sequence_id_to_length:
            is_length_column_av = is_length_column_available(read_handler)
            if not is_length_column_av:
                if not fastx_file:
                    exit("Sequences length could not be determined. Please provide a FASTA or FASTQ file using option -f or add column _LENGTH to gold standard.")
                binning_classes.Bin.sequence_id_to_length = read_lengths_from_fastx_file(fastx_file)

        g_sequence_ids = {}
        t_sequence_ids = {}
        if gold_standard:
            if gold_standard.genome_query:
                g_sequence_ids = gold_standard.genome_query.get_sequence_ids()
            if gold_standard.taxonomic_query:
                t_sequence_ids = gold_standard.taxonomic_query.get_sequence_ids()

        try:
            for sequence_id, bin_id, tax_id, length in read_binning_file(read_handler, is_gs):
                if is_gs and is_length_column_av:
                    binning_classes.Bin.sequence_id_to_length[sequence_id] = length

                if bin_id:
                    if not is_gs and sequence_id in g_sequence_ids or binning_classes.Bin.sequence_id_to_length[sequence_id] >= min_length:
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
                    if not is_gs and sequence_id not in t_sequence_ids:
                        continue
                    if not tax_id_to_parent:
                        if is_gs:
                            continue
                        else:
                            exit("Taxonomic binning cannot be assessed. Please provide an NCBI nodes file using option --ncbi_nodes_file.")
                    tax_id_path = load_ncbi_taxinfo.get_id_path(tax_id, tax_id_to_parent, tax_id_to_rank)
                    # TODO: tax id not found in ncbi or has "no rank". handle properly.
                    if not tax_id_path:
                        continue
                    for tax_id in tax_id_path:
                        if not tax_id: # tax_id may be empty
                            continue
                        if tax_id not in t_query.get_bin_ids():
                            bin = binning_classes.TaxonomicBin(tax_id)
                            t_query.add_bin(bin)
                        else:
                            bin = t_query.get_bin_by_id(tax_id)
                        t_query.rank_to_sequence_id_to_bin_id = (tax_id_to_rank[tax_id], sequence_id, tax_id)
                        bin.rank = tax_id_to_rank[tax_id]
                        bin.add_sequence_id(sequence_id)
        except BaseException as e:
            traceback.print_exc()
            exit("Error. File {} is malformed. {}".format(file_path_query, e))

    if not g_query.bins:
        g_query = None
    if not t_query.bins:
        t_query = None
    return g_query, t_query


def load_queries(gold_standard_file, fastx_file, query_files, map_by_completeness, filter_tail_percentage,
                 filter_genomes_file, filter_keyword, ncbi_nodes_file, min_length, labels):
    if not min_length:
        min_length = 0

    if ncbi_nodes_file:
        tax_id_to_parent, tax_id_to_rank = load_ncbi_taxinfo.load_tax_info(ncbi_nodes_file)
    else:
        tax_id_to_parent = tax_id_to_rank = None

    g_gold_standard, t_gold_standard = open_query(gold_standard_file,
                                                  True,
                                                  fastx_file,
                                                  tax_id_to_parent, tax_id_to_rank,
                                                  None,
                                                  min_length)
    g_gold_standard.filter_genomes_file = filter_genomes_file
    g_gold_standard.filter_keyword = filter_keyword
    gold_standard = binning_classes.GoldStandard(g_gold_standard, t_gold_standard)

    queries_list = []
    for query_file, label in zip(query_files, labels):
        g_query, t_query = open_query(query_file,
                                      False,
                                      None,
                                      tax_id_to_parent, tax_id_to_rank,
                                      gold_standard,
                                      0)
        if g_query:
            g_query.label = label
            g_query.map_by_completeness = map_by_completeness
            g_query.filter_tail_percentage = filter_tail_percentage
            g_query.filter_genomes_file = filter_genomes_file
            g_query.filter_keyword = filter_keyword
            queries_list.append(g_query)
        if t_query:
            t_query.label = label
            queries_list.append(t_query)

    # TODO if there is a g_query (t_query), there must be a g_gold_standard (t_gold_standard)

    return gold_standard, queries_list
