#!/usr/bin/env python

import argparse
import sys
import os

try:
    import argparse_parents
    import load_data
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import argparse_parents
        import load_data
    finally:
        sys.path.remove(os.path.dirname(__file__))


def add_column(mapping_file, fasta_file):
    lengths = load_data.read_lengths_from_fastx_file(fasta_file)
    sequence_id_column_index = 0
    with open(mapping_file, 'r') as read_handler:
        for line in read_handler:
            line = line.rstrip('\n')
            if line.startswith('@@'):
                for index, column_name in enumerate(line[2:].split('\t')):
                    if column_name == 'SEQUENCEID':
                        sequence_id_column_index = index
                print('%s\t_LENGTH' % line)
                continue

            if len(line.strip()) == 0 or line.startswith('#') or line.startswith('@'):
                print(line)
                continue

            row_data = line.split('\t')
            sequence_id = row_data[sequence_id_column_index]
            print('%s\t%d' % (line, lengths[sequence_id]))


def main():
    parser = argparse.ArgumentParser(description="Add length column _LENGTH to gold standard mapping and print mapping on the standard output")
    parser.add_argument("-g", "--gold_standard_file", help=argparse_parents.HELP_GOLD_STANDARD_FILE, required=True)
    parser.add_argument("-f", "--fasta_file", help="FASTA or FASTQ file with sequences of gold standard", required=True)
    args = parser.parse_args()
    add_column(args.gold_standard_file, args.fasta_file)


if __name__ == "__main__":
    main()