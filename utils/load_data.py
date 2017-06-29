#!/usr/bin/env python

import io
import os
import gzip
import mimetypes
from Bio import SeqIO
import numpy as np


class GoldStandard:
    def __init__(self):
        pass


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


def get_genome_mapping_without_lenghts(mapping_file):
    """
    @param mapping_file:
    @type mapping_file: str | unicode

    @return:
    """
    genome_id_to_list_of_contigs = {}
    anonymous_contig_id_to_genome_id = {}

    with open(mapping_file, 'r') as read_handler:
        for anonymous_contig_id, genome_id, length in read_binning_file(read_handler):
            anonymous_contig_id_to_genome_id[anonymous_contig_id] = genome_id
            if genome_id not in genome_id_to_list_of_contigs:
                genome_id_to_list_of_contigs[genome_id] = []
            genome_id_to_list_of_contigs[genome_id].append(anonymous_contig_id)

    return genome_id_to_list_of_contigs, anonymous_contig_id_to_genome_id


def get_genome_mapping(mapping_file, fastx_file):
    """
    @param mapping_file:
    @type mapping_file: str | unicode

    @return:
    """
    gold_standard = GoldStandard()
    gold_standard.genome_id_to_total_length = {}
    gold_standard.genome_id_to_list_of_contigs = {}
    gold_standard.contig_id_to_genome_id = {}
    gold_standard.contig_id_to_lengths = {}

    with open(mapping_file, 'r') as read_handler:
        is_length_column_av = is_length_column_available(read_handler)
        sequence_length = {}
        if not is_length_column_av:
            if not fastx_file:
                raise RuntimeError("Sequences' length could not be determined")
            sequence_length = read_lengths_from_fastx_file(fastx_file)
            # print(sorted(sequence_length.keys()))
        for anonymous_contig_id, genome_id, length in read_binning_file(read_handler):
            total_length = length if is_length_column_av else sequence_length[anonymous_contig_id]
            gold_standard.contig_id_to_lengths[anonymous_contig_id] = total_length
            gold_standard.contig_id_to_genome_id[anonymous_contig_id] = genome_id
            if genome_id not in gold_standard.genome_id_to_total_length:
                gold_standard.genome_id_to_total_length[genome_id] = 0
                gold_standard.genome_id_to_list_of_contigs[genome_id] = []
            gold_standard.genome_id_to_total_length[genome_id] += total_length
            gold_standard.genome_id_to_list_of_contigs[genome_id].append(anonymous_contig_id)

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


def read_rows(input_stream, index_key, index_value, index_length):
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip('\n')
        row_data = line.split('\t')
        key = row_data[index_key]
        value = row_data[index_value]
        if index_length is not None:
            length = row_data[index_length]
            yield key, value, length
        else:
            yield key, value, 0


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
    index_key = column_names["SEQUENCEID"]
    if "BINID" in column_names:
        index_value = column_names["BINID"]
    elif "TAXID" in column_names:
        index_value = column_names["TAXID"]
    else:
        raise RuntimeError("Column not found: {}".format("BINID/TAXID"))
    index_length = None
    if "_LENGTH" in column_names:
        index_length = column_names["_LENGTH"]
    return index_key, index_value, index_length


def read_binning_file(input_stream):
    index_key, index_value, index_length = get_column_indices(input_stream)
    return read_rows(input_stream, index_key, index_value, index_length)


def open_query(file_path_query):
    bin_id_to_list_of_sequence_id = {}
    sequence_id_to_bin_id = {}
    with open(file_path_query) as read_handler:
        for sequence_id, predicted_bin, length in read_binning_file(read_handler):
            if predicted_bin not in bin_id_to_list_of_sequence_id:
                bin_id_to_list_of_sequence_id[predicted_bin] = []
            bin_id_to_list_of_sequence_id[predicted_bin].append(sequence_id)
            sequence_id_to_bin_id[sequence_id] = predicted_bin
    return bin_id_to_list_of_sequence_id, sequence_id_to_bin_id
