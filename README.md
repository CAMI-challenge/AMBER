[![CircleCI](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master)

# Requirements

* The examples below require the gold standard assembly from
[https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz](https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz). 
Please download it to the _test_ directory.
* python2

# User Guide

## precision_recall.py

~~~BASH
usage: precision_recall.py [-h] [-l LABELS] -g GOLD_STANDARD_FILE
                           [-f FASTA_FILE] [-p FILTER] [-r GENOMES_FILE]
                           [-k KEYWORD]
                           query_files [query_files ...]

Compute precision and recall, including standard deviation and standard error
of the mean, for binning files

positional arguments:
  query_files           Query files

optional arguments:
  -h, --help            show this help message and exit
  -l LABELS, --labels LABELS
                        Comma-separated binning names
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file w/ sequences of gold standard -
                        required if gold standard file misses column _LENGTH
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest bins - default is 0
  -r GENOMES_FILE, --genomes_file GENOMES_FILE
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in second column of input for bins to be
                        removed (no keyword=remove all in list)
~~~
**Example:**
~~~BASH
./precision_recall.py -g /home/fmeyer/cami/data/gs_low/gsa_mapping.bin \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
-r test/unique_common.tsv -k "circular element" \
-p 1 \
-l "MaxBin 2.0, CONCOCT, MetaBAT" \
test/naughty_carson_2 test/goofy_hypatia_2 test/elated_franklin_0
~~~
**Output:**
~~~BASH
tool       precision std_dev_precision sem_precision recall std_dev_recall sem_recall
MaxBin 2.0 0.948     0.095             0.016         0.799  0.364          0.058
CONCOCT    0.837     0.266             0.052         0.517  0.476          0.069
MetaBAT    0.822     0.256             0.047         0.570  0.428          0.065
~~~

## precision_recall_per_genome.py

~~~BASH
usage: precision_recall_per_genome.py [-h] -g GOLD_STANDARD_FILE
                                      [-f FASTA_FILE]
                                      query_file

Compute table of precision and recall per genome bin

positional arguments:
  query_file            Query file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file w/ sequences of gold standard -
                        required if gold standard file misses column _LENGTH
~~~
**Example:**
~~~BASH
./precision_recall_per_genome.py -g /home/fmeyer/cami/data/gs_low/gsa_mapping.bin \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2
~~~
**Output:**
~~~BASH
genome          precision      recall  predicted_size correctly_predicted real_size
1049005         0.998458014811 1.0     4114177        4107833             4107833
evo_1035930.032 1.0            1.0     2294831        2294831             2294831
Sample18_8      0.409809009246 1.0     121786         49909               49909
evo_1035930.029 0.995223021915 1.0     2423708        2412130             2412130
...
~~~

## utils/exclude_genomes.py
~~~BASH
usage: exclude_genomes.py [-h] -r GENOMES_FILE [-k KEYWORD] [file]

Exclude genome bins from table file of precision and recall or standard input

positional arguments:
  file                  File containing precision and recall for each genome

optional arguments:
  -h, --help            show this help message and exit
  -r GENOMES_FILE, --genomes_file GENOMES_FILE
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in second column of input for bins to be
                        removed (no keyword=remove all in list)
~~~
**Example:**

The example computes the table of precision and recall and pipes it to utils/exclude_genomes.py.
~~~BASH
./precision_recall_per_genome.py -g /home/fmeyer/cami/data/gs_low/gsa_mapping.bin \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2 | \
./utils/exclude_genomes.py -r test/unique_common.tsv -k "circular element"
~~~

**Output:**

The output the is the table from precision_recall_per_genome.py without the excluded genomes.

## precision_recall_average.py
~~~BASH
usage: precision_recall_average.py [-h] [-p FILTER] [-l LABEL] [file]

Compute precision and recall, including standard deviation and standard error
of the mean, from table of precision and recall per genome provided as a file
or via the standard input

positional arguments:
  file                  File containing precision and recall for each genome

optional arguments:
  -h, --help            show this help message and exit
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest bins - default is 0
  -l LABEL, --label LABEL
                        Binning name                        
~~~
**Example:**
~~~BASH
./precision_recall_per_genome.py -g /home/fmeyer/cami/data/gs_low/gsa_mapping.bin \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2 | \
./utils/exclude_genomes.py -r test/unique_common.tsv -k "circular element" | \
./precision_recall_average.py -p 1 -l "MaxBin 2.0"
~~~
**Output:**
~~~BASH
tool       precision std_dev_precision sem_precision recall std_dev_recall sem_recall
MaxBin 2.0 0.948     0.095             0.016         0.799  0.364          0.058
~~~

## precision_recall_by_bpcount.py
~~~BASH
usage: precision_recall_by_bpcount.py [-h] -g GOLD_STANDARD_FILE
                                      [-f FASTA_FILE]
                                      query_file

Compute precision and recall weighed by base pair counts - not averaged over
genome bins - from binning file

positional arguments:
  query_file            Query file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file w/ sequences of gold standard -
                        required if gold standard file misses column _LENGTH
~~~
**Example:**
~~~BASH
./precision_recall_by_bpcount.py -g test/gsa_mapping.bin \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
/home/fmeyer/cami/data/binning/naughty_carson_2
~~~
**Output:**
~~~BASH
precision recall
0.934     0.838
~~~

## ari.py
~~~BASH
usage: ari.py [-h] -g GOLD_STANDARD_FILE query_file

Compute adjusted rand index from binning file

positional arguments:
  query_file            Query file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        gold standard - ground truth - file
~~~
**Example:**
~~~BASH
./ari.py -g test/gsa_mapping.bin test/naughty_carson_2
~~~
**Output:**
~~~BASH
0.782
~~~

## genome_recovery.py
~~~BASH
usage: genome_recovery.py [-h] [-p FILTER] [-l LABEL] [file]

Compute number of genomes in ranges of completeness and contamination

positional arguments:
  file                  File containing precision and recall for each genome

optional arguments:
  -h, --help            show this help message and exit
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest bins - default is 0
  -l LABEL, --label LABEL
                        Binning name
~~~
**Example:**
~~~BASH
./precision_recall_per_genome.py -g test/gsa_mapping.bin \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2 | \
./genome_recovery.py -l "MaxBin 2.0" -p 1
~~~
**Output:**
~~~BASH
MaxBin 2.0         >50% complete >70% complete >90% complete
<10% contamination 28            28            24
<5% contamination  23            23            21
~~~


# Developer Guide

We are using [tox]((https://tox.readthedocs.io/en/latest/)) for project automation.

### Tests

If you want to run tests just type tox in project dir:

~~~BASH
tox
~~~