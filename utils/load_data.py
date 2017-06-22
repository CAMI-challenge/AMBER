#!/usr/bin/python

import io
from utils import compression_handler


class GoldStandard:
    def __init__(self):
        pass


def read_lengths_from_fastx_file(fastx_file):
    length = {}
    f = compression_handler.get_compressed_file(fastx_file)
    br = io.BufferedReader(f.accessor)

    is_fastq = False
    for line in br:
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            is_fastq = True
        elif line.startswith(">"):
            br.seek(0)
            break

    if is_fastq:
        for line in br:
            line = line.strip()
            if not line:
                continue
            if line.startswith('@'):
                sequence_id = line[1:].strip()
                sequence = f.accessor.readline().strip()
                length[sequence_id] = len(sequence)
                next(f.accessor)
    else:
        for line in br:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                sequence_id = line[1:]
                if sequence_id not in length:
                    length[sequence_id] = 0
                continue
            sequence = line
            length[sequence_id] += len(sequence)
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
    gold_standard.sequence_id_to_genome_id = {}
    gold_standard.contig_id_to_lengths = {}

    with open(mapping_file, 'r') as read_handler:
        is_length_column_av = is_length_column_available(read_handler)
        if not is_length_column_av:
            if not fastx_file:
                raise RuntimeError("Sequences' length could not be determined")
            sequence_length = read_lengths_from_fastx_file(fastx_file)

        for anonymous_contig_id, genome_id, length in read_binning_file(read_handler):
            total_length = length if is_length_column_av else sequence_length[anonymous_contig_id]
            gold_standard.contig_id_to_lengths[anonymous_contig_id] = total_length
            gold_standard.sequence_id_to_genome_id[anonymous_contig_id] = genome_id
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
