#!/usr/bin/python

# Script for computing precision and recall. It takes as input:
# - a gold standard file in bioboxes format
# (https://github.com/bioboxes/rfc/blob/4bb19a633a6a969c2332f1f298852114c5f89b1b/data-format/binning.mkd)
# with optional column _LENGTH
# - a (compressed) fasta or fastq file, required if _LENGTH is not present in the gold standard file
# - the bins to be evaluated in the same format as above
# It outputs a tsv file containing precision and recall for each bin.

import argparse
import io
from collections import Counter
from collections import defaultdict
from exceptions import RuntimeError
from utils import compression_handler


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


def get_genome_mapping(mapping_file, fastx_file):
    """
    @param mapping_file:
    @type mapping_file: str | unicode

    @return:
    """
    genome_id_to_total_length = {}
    genome_id_to_list_of_contigs = {}
    anonymous_contig_id_to_genome_id = {}
    anonymous_contig_id_to_lengths = {}

    with open(mapping_file, 'r') as read_handler:
        is_length_column_av = is_length_column_available(read_handler)
        if not is_length_column_av:
            if not fastx_file:
                raise RuntimeError("Sequences' length could not be determined")
            sequence_length = read_lengths_from_fastx_file(fastx_file)

        for anonymous_contig_id, genome_id, length in read_binning_file(read_handler):
            total_length = length if is_length_column_av else sequence_length[anonymous_contig_id]
            anonymous_contig_id_to_lengths[anonymous_contig_id] = total_length
            anonymous_contig_id_to_genome_id[anonymous_contig_id] = genome_id
            if genome_id not in genome_id_to_total_length:
                genome_id_to_total_length[genome_id] = 0
                genome_id_to_list_of_contigs[genome_id] = []
            genome_id_to_total_length[genome_id] += total_length
            genome_id_to_list_of_contigs[genome_id].append(anonymous_contig_id)

    return genome_id_to_total_length, genome_id_to_list_of_contigs, anonymous_contig_id_to_genome_id, anonymous_contig_id_to_lengths


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


def calcPrecisionRecall(bin_id_to_genome_id_to_total_length, bin_id_to_total_lengths, genome_id_to_total_length):
    # now calculate precision looping bins
    true_positives_p = 0
    denominator_p = sum(bin_id_to_total_lengths.values())
    false_positives_p = 0.0

    for predicted_bin, genome_assigns in bin_id_to_genome_id_to_total_length.iteritems():
        # get maximal genome assignment
        if len(genome_assigns) > 0:
            maxAssign = max(genome_assigns.values())
        else:
            maxAssign = 0.0

        true_positives_p += maxAssign

    # now calculate precision as TP/FP + TP
    precision = true_positives_p / denominator_p

    # now calculate recall looping genomes
    true_positives_r = 0.0
    false_negatives_r = 0.0
    for genome_id in genome_id_to_total_length:
        # now loop bins
        bin_assigns = []
        for bin_id in bin_id_to_total_lengths:
            if genome_id in bin_id_to_genome_id_to_total_length[bin_id]:
                bin_assigns.append(bin_id_to_genome_id_to_total_length[bin_id][genome_id])
        if len(bin_assigns) > 0:
            maxAssign = max(bin_assigns)
        else:
            maxAssign = 0.0
        true_positives_r += maxAssign
        false_negatives_r += genome_id_to_total_length[genome_id] - maxAssign

    denominator_r = sum(genome_id_to_total_length.values())
    recall = true_positives_r / denominator_r

    return (precision, recall)


def validate_genomes(file_path_query, file_path_mapping, file_path_output):
    """
    This script calculates precision, recall, accuracy and ARI for binned contigs
    
    @param file_path_query:
    @param file_path_mapping:
    @param file_path_output:
    @return:
    """
    # import ipdb; ipdb.set_trace()
    genome_id_to_total_length, genome_id_to_list_of_contigs, sequence_id_to_genome_id, anonymouse_contig_id_to_length = get_genome_mapping(
        file_path_mapping)
    bin_id_to_list_of_sequence_id = {}
    bin_id_to_total_length = {}
    with open(file_path_query) as read_handler:
        for sequence_id, predicted_bin, length in read_binning_file(read_handler):
            if predicted_bin not in bin_id_to_total_length:
                bin_id_to_list_of_sequence_id[predicted_bin] = []
                bin_id_to_total_length[predicted_bin] = 0
            bin_id_to_list_of_sequence_id[predicted_bin].append(sequence_id)
            bin_id_to_total_length[predicted_bin] += anonymouse_contig_id_to_length[sequence_id]

        # now compute confusion matrix as dictionary of Counters
    bin_id_to_genome_id_to_total_length = defaultdict(Counter)
    for predicted_bin in bin_id_to_list_of_sequence_id:
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            genome_id = sequence_id_to_genome_id[sequence_id]
            bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += anonymouse_contig_id_to_length[sequence_id]

    (precision, recall) = calcPrecisionRecall(bin_id_to_genome_id_to_total_length, bin_id_to_total_length,
                                              genome_id_to_total_length)

    with open(file_path_output + "/by_genome_weighted.tsv", 'w') as write_handler:
        write_handler.write("Precision\tRecall\n")
        write_handler.write("%s \t %s\n" % (precision, recall))


def map_genomes(file_path_mapping, file_path_query, file_fasta):
    """
    This script mapps a predicted bin to the genome with the highest recall

    @attention: In case of reads, read ids might not be paired read id and cause error: ReadID/1 ReadID/2

    @param file_path_query:
    @param file_path_mapping:
    @return:
    """
    genome_id_to_total_length, genome_id_to_list_of_contigs, sequence_id_to_genome_id, anonymous_contig_id_to_lengths = get_genome_mapping(file_path_mapping, file_fasta)
    bin_id_to_list_of_sequence_id = {}
    bin_id_to_total_lengths = {}

    with open(file_path_query) as read_handler:
        for sequence_id, predicted_bin, length in read_binning_file(read_handler):
            if predicted_bin not in bin_id_to_total_lengths:
                bin_id_to_list_of_sequence_id[predicted_bin] = []
                bin_id_to_total_lengths[predicted_bin] = 0
            bin_id_to_list_of_sequence_id[predicted_bin].append(sequence_id)
            bin_id_to_total_lengths[predicted_bin] += anonymous_contig_id_to_lengths[sequence_id]

    bin_metrics = []
    bin_id_to_genome_id_to_total_length = {}
    mapped = set()

    for predicted_bin in bin_id_to_list_of_sequence_id:
        bin_id_to_genome_id_to_total_length[predicted_bin] = {}
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            genome_id = sequence_id_to_genome_id[sequence_id]
            if genome_id not in bin_id_to_genome_id_to_total_length[predicted_bin]:
                bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] = 0
            bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += anonymous_contig_id_to_lengths[sequence_id]
        max_length = 0
        best_genome_id = ""
        for genome_id in bin_id_to_genome_id_to_total_length[predicted_bin]:
            if max_length < bin_id_to_genome_id_to_total_length[predicted_bin][genome_id]:
                max_length = bin_id_to_genome_id_to_total_length[predicted_bin][genome_id]
                best_genome_id = genome_id
        mapped.add(best_genome_id)
        # length of genome in bin divided by bin size
        precision = float(bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id]) / float(bin_id_to_total_lengths[predicted_bin])
        recall = float(bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id]) / float(genome_id_to_total_length[best_genome_id])
        bin_metrics.append([predicted_bin, best_genome_id, precision, recall,
                                      bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id],
                                      genome_id_to_total_length[best_genome_id]])

    bin_metrics = sorted(bin_metrics, key=lambda t: t[3], reverse=True)

    print "@@genome\tprecision\trecall\tpredicted_size\tcorrectly_predicted\treal_size"
    for metrics in bin_metrics:
        print "%s\t%s\t%s\t%s\t%s\t%s" % (
            metrics[1],
            metrics[2],
            metrics[3],
            bin_id_to_total_lengths[metrics[0]],
            metrics[4],
            metrics[5])
    for genome_id in genome_id_to_list_of_contigs:
        if genome_id not in mapped:
            print "%s\t%s\t%s\t%s\t%s\t%s" % (
            genome_id, 'NA', .0, 0, 0, genome_id_to_total_length[genome_id])  # precision is NA for unpredicted bins


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gold_standard_file", help="gold standard - ground truth - file", required=True)
    parser.add_argument("-q", "--query_file", help="query file", required=True)
    parser.add_argument("-f", "--fast_file",
                        help="FASTA or FASTQ file w/ sequences of gold standard - required if gold standard file misses column _LENGTH")
    args = parser.parse_args()
    if not args.gold_standard_file or not args.query_file:
        parser.print_help()
        parser.exit(1)
    map_genomes(file_path_mapping=args.gold_standard_file,
                file_path_query=args.query_file,
                file_fasta=args.fast_file)


if __name__ == "__main__":
    main()


# validate_genomes(
#     file_path_query=sys.argv[1],
#     file_path_mapping=sys.argv[2],
#     file_path_output=sys.argv[3]
# )
