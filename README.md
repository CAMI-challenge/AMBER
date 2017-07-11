[![CircleCI](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master)

# Requirements

* The examples below require the gold standard assembly from
[https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz](https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz). 
Please download it to the _test_ directory.
* python3.5
* tox (only for automatic tests)

# User Guide

## List of metrics and abbreviations

* **avg_precision**: precision averaged over genome bins
* **std_dev_precision**: standard deviation of precision averaged over genome bins
* **sem_precision**: standard error of the mean of precision averaged over genome bins
* **avg_recall**: recall averaged over genome bins
* **std_dev_recall**: standard deviation of recall averaged over genome bins
* **sem_recall**: standard error of the mean of recall averaged over genome bins
* **precision**: precision weighed by base pairs
* **recall**: recall weighed by base pairs
* **rand_index_by_bp**: rand index weighed by base pairs
* **rand_index_by_seq**: rand index weighed by sequence counts
* **a_rand_index_by_bp**: adjusted rand index weighed by base pairs
* **a_rand_index_by_seq**: adjusted rand index weighed by sequence counts
* **percent_assigned_bps**: percentage of base pairs that were assigned to bins
* **\>0.5compl<0.1cont**: number of genomes with more than 50% completeness and less than 10% contamination
* **\>0.7compl<0.1cont**: number of genomes with more than 70% completeness and less than 10% contamination
* **\>0.9compl<0.1cont**: number of genomes with more than 90% completeness and less than 10% contamination
* **\>0.5compl<0.05cont**: number of genomes with more than 50% completeness and less than 5% contamination
* **\>0.7compl<0.05cont**: number of genomes with more than 70% completeness and less than 5% contamination
* **\>0.9compl<0.05cont**: number of genomes with more than 90% completeness and less than 5% contamination

## evaluate.py
~~~BASH
usage: evaluate.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE] [-l LABELS]
                   [-p FILTER] [-r GENOMES_FILE] [-k KEYWORD] -o OUTPUT_DIR
                   query_files [query_files ...]

Compute all metrics for binning files; output summary to screen and results
per query binning file to chosen directory

positional arguments:
  query_files           Query binning files

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard -
                        required if gold standard file misses column _LENGTH
  -l LABELS, --labels LABELS
                        Comma-separated binning names
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest bins - default is 0
  -r GENOMES_FILE, --genomes_file GENOMES_FILE
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in second column of input for bins to be
                        removed (no keyword=remove all in list)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write the results per query to
~~~
**Example:**
~~~BASH
./evaluate.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
-l "MaxBin 2.0, CONCOCT, MetaBAT" \
-p 1 \
-r test/unique_common.tsv \
-k "circular element" \
test/naughty_carson_2 \
test/goofy_hypatia_2 \
test/elated_franklin_0 \
-o output_dir/
~~~
**Output:**
~~~BASH
tool       avg_precision std_dev_precision sem_precision avg_recall std_dev_recall sem_recall precision recall rand_index_by_bp rand_index_by_seq a_rand_index_by_bp a_rand_index_by_seq percent_assigned_bps >0.5compl<0.1cont >0.7compl<0.1cont >0.9compl<0.1cont >0.5compl<0.05cont >0.7compl<0.05cont >0.9compl<0.05cont
MaxBin 2.0 0.948         0.095             0.016         0.799      0.364          0.058      0.934     0.838  0.995            0.951             0.917              0.782               0.864                28                28                24                23                 23                 21
CONCOCT    0.837         0.266             0.052         0.517      0.476          0.069      0.684     0.936  0.972            0.946             0.644              0.751               0.967                18                17                15                16                 16                 14
MetaBAT    0.822         0.256             0.047         0.57       0.428          0.065      0.724     0.825  0.976            0.965             0.674              0.860               0.917                17                16                12                17                 16                 12
~~~
Additionally, in directory _output_dir_, directories _naughty_carson_2_, _goofy_hypatia_2_, and _elated_franklin_0_ are created with the following files:
* _rand_index.tsv_: contains value of (adjusted) rand index and percentage of assigned/binned bases. Rand index is weighed and unweighed by base pairs
* _precision_recall.tsv_: contains precision and recall per genome bin
* _precision_recall_avg.tsv_: contains precision and recall averaged over genome bins. Includes standard deviation and standard error of the mean
* _precision_recall_by_bpcount.tsv_: contains precision and recall weighed by base pairs

## precision_recall.py
~~~BASH
usage: precision_recall.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE]
                           [-l LABELS] [-p FILTER] [-r GENOMES_FILE]
                           [-k KEYWORD]
                           query_files [query_files ...]

Compute precision and recall, including standard deviation and standard error
of the mean, for binning files

positional arguments:
  query_files           Query binning files

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file w/ sequences of gold standard -
                        required if gold standard file misses column _LENGTH
  -l LABELS, --labels LABELS
                        Comma-separated binning names
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
./precision_recall.py -g test/gsa_mapping.binning \
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
  query_file            Query binning file

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
./precision_recall_per_genome.py -g test/gsa_mapping.binning \
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

Exclude genome bins from table of precision and recall per genome. The table
can be provided as file or via the standard input

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
./precision_recall_per_genome.py -g test/gsa_mapping.binning \
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
of the mean, from table of precision and recall per genome. The table can be
provided as file or via the standard input

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
./precision_recall_per_genome.py -g test/gsa_mapping.binning \
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

Compute precision and recall weighed by base pair counts (not averaged over
genome bins) from binning file

positional arguments:
  query_file            Query binning file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard -
                        required if gold standard file misses column _LENGTH
~~~
**Example:**
~~~BASH
./precision_recall_by_bpcount.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2
~~~
**Output:**
~~~BASH
precision recall
0.934     0.838
~~~

## rand_index.py
~~~BASH
usage: rand_index.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE] query_file

Compute (adjusted) rand index from binning file, unweighed and weighed by base
pairs, and percentage of binned base pairs

positional arguments:
  query_file            Query binning file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard -
                        required if gold standard file misses column _LENGTH
~~~
**Example:**
~~~BASH
./rand_index.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2
~~~
**Output:**
~~~BASH
rand_index_by_bp rand_index_by_seq a_rand_index_by_bp a_rand_index_by_seq percent_assigned_bps
0.995            0.951             0.917              0.782               0.864
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
./precision_recall_per_genome.py -g test/gsa_mapping.binning \
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


## Run the tool as a Biobox

1. Build the container running the following command:

~~~BASH
docker build  -t genome_binning_evaluation .
~~~

2. Run the container by specifying the executing the following command

~~~BASH
docker run -v $(pwd)/input/gold_standard.fasta:/bbx/input/gold_standard.fasta -v $(pwd)/input/gsa_mapping.binning:/bbx/input/gsa_mapping.binning  -v  $(pwd)/input/test_query.binning:/bbx/input/test_query.binning  -v  $(pwd)/output:/bbx/output -v $(pwd)/input/biobox.yaml:/bbx/input/biobox.yaml genome_binning_evaluation default
~~~

where biobox.yaml contains the following values:

~~~YAML
version: 0.11.0
arguments:
  - fasta:
      value: /bbx/input/gold_standard.fasta
      type: contig
  - labels:
      value: /bbx/input/gsa_mapping.binning
      type: binning
  - predictions:
      value: /bbx/input/test_query.binning
      type: binning
~~~ 

# Developer Guide

We are using [tox]((https://tox.readthedocs.io/en/latest/)) for project automation.

### Tests

If you want to run tests, just type _tox_ in the project's root directory:

~~~BASH
tox
~~~
