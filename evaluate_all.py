#!/usr/bin/python

import argparse
import precision_recall_per_genome
import precision_recall_average
from utils import exclude_genomes


def evaluate_all(gold_standard_file, query_files, fasta_file, filter_tail_bool, genomes_file, keyword):
    for query_file in query_files:
        print "Evaluating: %s" % query_file
        bin_metrics = precision_recall_per_genome.compute_metrics(gold_standard_file,
                                                                  query_file,
                                                                  fasta_file)

        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)

        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_bool)

        print "PRECISION"
        print "Precision:\t\t%1.3f" % avg_precision
        print "Standard deviation:\t%1.3f" % std_deviation_precision
        print "Standard error of mean:\t%1.3f" % std_error_precision
        print "RECALL"
        print "Recall:\t\t\t%1.3f" % avg_recall
        print "Standard deviation:\t%1.3f" % std_deviation_recall
        print "Standard error of mean:\t%1.3f" % std_error_recall

    return 0


def main():
    parser = argparse.ArgumentParser(description="Compute metrics for multiple binning files at once")
    parser.add_argument("query_files", nargs='+', help="Query files")
    parser.add_argument("-g", "--gold_standard_file", help="gold standard - ground truth - file", required=True)
    parser.add_argument("-f", "--fasta_file",
                        help="FASTA or FASTQ file w/ sequences of gold standard - required if gold standard file misses column _LENGTH")
    parser.add_argument('-s', '--filter', action="store_true", help="Filter out 1%% smallest bins - default is false")
    parser.add_argument('-r', '--genomes_file', help="File with list of genomes to be removed", required=False)
    parser.add_argument('-k', '--keyword', help="Keyword in second column of input for bins to be removed (no keyword=remove all in list)", required=False)
    args = parser.parse_args()
    evaluate_all(args.gold_standard_file,
                 args.query_files,
                 args.fasta_file,
                 args.filter,
                 args.genomes_file,
                 args.keyword)


if __name__ == "__main__":
    main()
