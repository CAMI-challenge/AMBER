[![CircleCI](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master)

# Scripts for computing precision and recall. 

## Requirements

* python2

## User Guide

It takes as input:

* a gold standard file in bioboxes format
(https://github.com/bioboxes/rfc/blob/4bb19a633a6a969c2332f1f298852114c5f89b1b/data-format/binning.mkd)
with optional column _LENGTH

* a (compressed) fasta or fastq file, required if _LENGTH is not present in the gold standard file

* the bins to be evaluated in the same format as above
It writes to standard output a table containing precision and recall for each bin.

~~~BASH
usage: precision_recall_per_genome.py [-h] -g GOLD_STANDARD_FILE -q QUERY_FILE
                                   [-f FAST_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        gold standard - ground truth - file
  -q QUERY_FILE, --query_file QUERY_FILE
                        query file
  -f FAST_FILE, --fast_file FAST_FILE
                        FASTA or FASTQ file w/ sequences of gold standard -
                        required if gold standard file misses column _LENGTH
~~~

### Example:

Download gold standard assembly from https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz

Then run:

~~~BASH
./precision_recall_per_genome.py -g test/gsa_mapping.bin -q test/admiring_curie_3 -f CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz
~~~


## Developer Guide

We are using [tox]((https://tox.readthedocs.io/en/latest/)) for project automation.

### Tests

If you want to run tests just type tox in project dir:

~~~BASH
tox
~~~