#!/usr/bin/env python

import argparse

HELP_GOLD_STANDARD_FILE = "Gold standard - ground truth - file"
HELP_FILE = "File containing precision and recall for each genome"
HELP_QUERY_FILE = "Binning file"
HELP_QUERY_FILES = "Binning files"
HELP_FASTA_FILE = "FASTA or FASTQ file with sequences of gold standard (required if gold standard file misses column _LENGTH)"
HELP_FILTER = "Filter out [FILTER]%% smallest bins (default: 0)"
HELP_GENOMES_FILE = "File with list of genomes to be removed"
HELP_LABEL = "Binning name"
HELP_LABELS = "Comma-separated binning names"
HELP_KEYWORD = "Keyword in second column of input for bins to be removed (no keyword=remove all in list)"
HELP_MAP_BY_RECALL = "Map genomes to bins by maximizing recall"

PARSER_GS = argparse.ArgumentParser(add_help=False)
PARSER_GS.add_argument("bin_file", help=HELP_QUERY_FILE)
PARSER_GS.add_argument("-g", "--gold_standard_file", help=HELP_GOLD_STANDARD_FILE, required=True)
PARSER_GS.add_argument("-f", "--fasta_file", help=HELP_FASTA_FILE)

PARSER_MULTI = argparse.ArgumentParser(add_help=False)
PARSER_MULTI.add_argument('file', nargs='?', type=argparse.FileType('r'), help=HELP_FILE)
PARSER_MULTI.add_argument('-p', '--filter', help=HELP_FILTER)
PARSER_MULTI.add_argument('-l', '--label', help=HELP_LABEL, required=False)

PARSER_MULTI2 = argparse.ArgumentParser(add_help=False)
PARSER_MULTI2.add_argument("bin_files", nargs='+', help=HELP_QUERY_FILES)
PARSER_MULTI2.add_argument("-g", "--gold_standard_file", help=HELP_GOLD_STANDARD_FILE, required=True)
PARSER_MULTI2.add_argument("-f", "--fasta_file", help=HELP_FASTA_FILE)
PARSER_MULTI2.add_argument('-l', '--labels', help=HELP_LABELS, required=False)
PARSER_MULTI2.add_argument('-p', '--filter', help=HELP_FILTER)
PARSER_MULTI2.add_argument('-r', '--genomes_file', help=HELP_GENOMES_FILE, required=False)
PARSER_MULTI2.add_argument('-k', '--keyword', help=HELP_KEYWORD, required=False)
