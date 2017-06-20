#!/usr/bin/python

# Script for computing precision and recall. It takes as input:
# - a gold standard file in bioboxes format
# (https://github.com/bioboxes/rfc/blob/4bb19a633a6a969c2332f1f298852114c5f89b1b/data-format/binning.mkd)
# with optional column _LENGTH
# - a (compressed) fasta or fastq file, required if _LENGTH is not present in the gold standard file
# - the bins to be evaluated in the same format as above
# It writes to standard output a table containing precision and recall for each bin.

import argparse
from utils import load_data


def map_genomes(sequence_id_to_genome_id, bin_id_to_list_of_sequence_id, anonymous_contig_id_to_lengths):
    """
        This script maps a predicted bin to the genome with the highest recall

        @attention: In case of reads, read ids might not be paired read id and cause error: ReadID/1 ReadID/2

        @param sequence_id_to_genome_id:
        @param anonymous_contig_id_to_lengths:
        @param bin_id_to_list_of_sequence_id
        @return:
        """
    bin_id_to_genome_id_to_total_length = {}
    mapped = set()
    bin_id_to_mapped_genome = {}
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
        bin_id_to_mapped_genome[predicted_bin] = best_genome_id
    return bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped


def compute_precision_recall(genome_id_to_total_length,
                             genome_id_to_list_of_contigs,
                             bin_id_to_list_of_sequence_id,
                             bin_id_to_mapped_genome,
                             bin_id_to_genome_id_to_total_length,
                             anonymous_contig_id_to_lengths,
                             mapped):
    bin_id_to_total_lengths = {}
    for predicted_bin in bin_id_to_list_of_sequence_id:
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            if predicted_bin not in bin_id_to_total_lengths:
                bin_id_to_total_lengths[predicted_bin] = 0
            bin_id_to_total_lengths[predicted_bin] += anonymous_contig_id_to_lengths[sequence_id]

    bin_metrics = []
    for predicted_bin in bin_id_to_list_of_sequence_id:
        best_genome_id = bin_id_to_mapped_genome[predicted_bin]
        # length of genome in bin divided by bin size
        precision = float(bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id]) / float(bin_id_to_total_lengths[predicted_bin])
        recall = float(bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id]) / float(genome_id_to_total_length[best_genome_id])
        bin_metrics.append({'mapped_genome': best_genome_id,
                            'precision': precision,
                            'recall': recall,
                            'predicted_size': bin_id_to_total_lengths[predicted_bin],
                            'correctly_predicted': bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id],
                            'real_size': genome_id_to_total_length[best_genome_id]})
    # sort bins by recall
    bin_metrics = sorted(bin_metrics, key=lambda t: t['recall'], reverse=True)

    print "@@genome\tprecision\trecall\tpredicted_size\tcorrectly_predicted\treal_size"
    for bin in bin_metrics:
        print "%s\t%s\t%s\t%s\t%s\t%s" % (
            bin['mapped_genome'],
            bin['precision'],
            bin['recall'],
            bin['predicted_size'],
            bin['correctly_predicted'],
            bin['real_size'])
    for genome_id in genome_id_to_list_of_contigs:
        if genome_id not in mapped:
            print "%s\t%s\t%s\t%s\t%s\t%s" % (
            genome_id, 'NA', .0, 0, 0, genome_id_to_total_length[genome_id])  # precision is NA for unpredicted bins


def compute_metrics(file_path_mapping, file_path_query, file_fasta):
    genome_id_to_total_length, genome_id_to_list_of_contigs, sequence_id_to_genome_id, anonymous_contig_id_to_lengths = \
        load_data.get_genome_mapping(file_path_mapping, file_fasta)
    bin_id_to_list_of_sequence_id, sequence_id_to_bin_id = load_data.open_query(file_path_query)
    bin_id_to_mapped_genome, bin_id_to_genome_id_to_total_length, mapped = map_genomes(sequence_id_to_genome_id,
                                                                                       bin_id_to_list_of_sequence_id,
                                                                                       anonymous_contig_id_to_lengths)
    compute_precision_recall(genome_id_to_total_length,
                             genome_id_to_list_of_contigs,
                             bin_id_to_list_of_sequence_id,
                             bin_id_to_mapped_genome,
                             bin_id_to_genome_id_to_total_length,
                             anonymous_contig_id_to_lengths,
                             mapped)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gold_standard_file", help="gold standard - ground truth - file", required=True)
    parser.add_argument("-q", "--query_file", help="query file", required=True)
    parser.add_argument("-f", "--fasta_file",
                        help="FASTA or FASTQ file w/ sequences of gold standard - required if gold standard file misses column _LENGTH")
    args = parser.parse_args()
    if not args.gold_standard_file or not args.query_file:
        parser.print_help()
        parser.exit(1)
    compute_metrics(file_path_mapping=args.gold_standard_file,
                    file_path_query=args.query_file,
                    file_fasta=args.fasta_file)


if __name__ == "__main__":
    main()
