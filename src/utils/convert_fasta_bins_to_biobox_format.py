#!/usr/bin/env python

import argparse


def read_fasta_file(fasta_file):
    with open(fasta_file, 'r') as read_handler:
        sequence_id = ""
        for line in read_handler:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(sequence_id) > 0:
                    yield sequence_id, fasta_file
                sequence_id = line[1:]
                continue
        if len(sequence_id) > 0:
            yield sequence_id, fasta_file


def convert(paths, output_file):
    with open(output_file, 'w') as write_handler:
        write_handler.write("#CAMI Format for Binning\n@Version:0.9.0\n@SampleID:_SAMPLEID_\n@@SEQUENCEID\tBINID\n")
        for file in paths:
            for sequence_id, bin_id in read_fasta_file(file):
                write_handler.write("%s\t%s\n" % (sequence_id, bin_id))


def main():
    parser = argparse.ArgumentParser(description="Convert bins in FASTA files to CAMI tsv format")
    parser.add_argument("paths", nargs='+', help="FASTA files including full paths")
    parser.add_argument("-o", "--output_file", required=False, help="Output file")
    args = parser.parse_args()
    convert(args.paths, args.output_file)


if __name__ == "__main__":
    main()