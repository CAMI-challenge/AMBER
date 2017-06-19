#!/usr/bin/python

import sys
import argparse


def load_unique_common(unique_common_file_path):
    genome_to_unique_common = {}
    with open(unique_common_file_path) as read_handler:
        for line in read_handler:
            genome_to_unique_common[line.split("\t")[0]] = line.split("\t")[1].strip('\n')
    return genome_to_unique_common


def load_data(stream, unique_common_file_path):
    genome_to_unique_common = load_unique_common(unique_common_file_path)
    for line in stream:
        line = line.strip()
        if len(line) == 0 or line.startswith("@"):
            continue
        bin = line.split('\t')[0]
        if bin in genome_to_unique_common and genome_to_unique_common[bin] == "circular element":
            continue
        print line


def main():
    parser = argparse.ArgumentParser(description="Exclude bins from file or standard input")
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), help="File containing precision and recall for each bin")
    parser.add_argument('-b', '--bins', help="File with list of bins", required=True)
    parser.add_argument('-k', '--keyword', help="Keyword in second column of input for bins to be removed (none=remove all in list)", required=False)
    args = parser.parse_args()
    print args.keyword
    if not args.bins and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    load_data(sys.stdin if not sys.stdin.isatty() else args.file,
              args.bins,
              args.keyword)


if __name__ == "__main__":
    main()