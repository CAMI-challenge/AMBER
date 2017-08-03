[![CircleCI](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/genome_binning_evaluation/tree/master)

# Introduction
AMBER (Assessment of Metagenome BinnERs) is an evaluation package for the comparative assessment of genome
reconstructions from metagenome benchmark datasets. It provides performance metrics, results rankings, and
comparative visualizations for assessing multiple programs or parameter effects. The provided metrics were
used in the first community benchmarking challenge of the initiative for the [Critical Assessment of Metagenomic
Interpretation](http://www.cami-challenge.org/).

## Inputs
As input, the main tool _evaluate.py_ of AMBER uses three files:
1. a gold standard mapping of contigs or read IDs to genomes in the
[CAMI binning Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format);
see [here](https://github.com/CAMI-challenge/genome_binning_evaluation/blob/master/test/gsa_mapping.binning) example (note: only columns SEQUENCEID and BINID are required)
2. one or more files with bin assignments for the sequences also in the
[CAMI binning Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format), with each file
containing all the bin assignments from the run of a binning program. A tool for converting FASTA files, such that each file represents a bin,
is available (see _utils/convert_fasta_bins_to_biobox_format.py_ below)
3. a FASTA or FASTQ file with the sequences for obtaining their lengths. Optionally, the lenghts may be added to the
gold standard mapping file using tool _utils/add_length_column.py_ (see below). In this way,
_evaluate.py_ no longer requires a FASTA or FASTQ file

Additional parameters may be specified - see below.

# Requirements

* python &ge; 3.5
* python-tk
* numpy &ge; v1.13.0
* biopython &ge; v1.69.0
* matplotlib &ge; v2.0.2
* tox (only for automatic tests)
* The examples below require a gold standard assembly. Please [download it](https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz) to the _test_ directory.

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
* **rand_index_by_bp**: Rand index weighed by base pairs
* **rand_index_by_seq**: Rand index weighed by sequence counts
* **a_rand_index_by_bp**: adjusted Rand index weighed by base pairs
* **a_rand_index_by_seq**: adjusted Rand index weighed by sequence counts
* **percent_assigned_bps**: percentage of base pairs that were assigned to bins
* **\>0.5compl<0.1cont**: number of bins with more than 50% completeness and less than 10% contamination
* **\>0.7compl<0.1cont**: number of bins with more than 70% completeness and less than 10% contamination
* **\>0.9compl<0.1cont**: number of bins with more than 90% completeness and less than 10% contamination
* **\>0.5compl<0.05cont**: number of bins with more than 50% completeness and less than 5% contamination
* **\>0.7compl<0.05cont**: number of bins with more than 70% completeness and less than 5% contamination
* **\>0.9compl<0.05cont**: number of bins with more than 90% completeness and less than 5% contamination

## evaluate.py
~~~BASH
usage: evaluate.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE] [-l LABELS]
                   [-p FILTER] [-r GENOMES_FILE] [-k KEYWORD] -o OUTPUT_DIR
                   bin_files [bin_files ...]

Compute all metrics and figures for one or more binning files; output summary
to screen and results per binning file to chosen directory

positional arguments:
  bin_files             Binning files

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard
                        (required if gold standard file misses column _LENGTH)
  -l LABELS, --labels LABELS
                        Comma-separated binning names
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest bins (default: 0)
  -r GENOMES_FILE, --genomes_file GENOMES_FILE
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in second column of input for bins to be
                        removed (no keyword=remove all in list)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write the results to
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
Additionally, directory _output_dir_ will contain figures **avg_precision_recall.png + .pdf** (average precision vs. average recall)
and **ari_vs_assigned_bps.png + .pdf** (adjusted Rand index vs. percentage of assigned base pairs), and **rankings.txt** (binnings sorted by average precision, average recall, and average precision + recall).
In the same directory, subdirectories _naughty_carson_2_, _goofy_hypatia_2_, and _elated_franklin_0_ will be created with the following files:
* **rand_index.tsv**: contains value of (adjusted) Rand index and percentage of assigned/binned bases. Rand index is both weighed and unweighed by base pairs
* **precision_recall.tsv**: contains precision and recall per genome bin
* **precision_recall_avg.tsv**: contains precision and recall averaged over genome bins. Includes standard deviation and standard error of the mean
* **precision_recall_by_bpcount.tsv**: contains precision and recall weighed by base pairs
* **genomes_sorted_by_precision.png + .pdf**: figure of precision and recall per genome with genomes sorted by precision
* **genomes_sorted_by_recall.png + .pdf**: figure of precision and recall per genome with genomes sorted by recall

## precision_recall.py
~~~BASH
usage: precision_recall.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE]
                           [-l LABELS] [-p FILTER] [-r GENOMES_FILE]
                           [-k KEYWORD]
                           bin_files [bin_files ...]

Compute precision and recall, including standard deviation and standard error
of the mean, for binning files

positional arguments:
  bin_files             Binning files

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard
                        (required if gold standard file misses column _LENGTH)
  -l LABELS, --labels LABELS
                        Comma-separated binning names
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest bins (default: 0)
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
                                      bin_file

Compute table of precision and recall per genome bin

positional arguments:
  bin_file              Binning file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard
                        (required if gold standard file misses column _LENGTH)
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
                        Filter out [FILTER]% smallest bins (default: 0)
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
                                      bin_file

Compute precision and recall weighed by base pair counts (not averaged over
genome bins) from binning file

positional arguments:
  bin_file              Binning file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard
                        (required if gold standard file misses column _LENGTH)
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
usage: rand_index.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE] bin_file

Compute (adjusted) Rand index from binning file, unweighed and weighed by base
pairs, and percentage of binned base pairs

positional arguments:
  bin_file              Binning file

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard
                        (required if gold standard file misses column _LENGTH)
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
                        Filter out [FILTER]% smallest bins (default: 0)
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

## plot_by_genome.py
~~~BASH
usage: plot_by_genome.py [-h] [-s {recall,precision}] [-o OUT_FILE] [file]

Plot precision and recall per genome. Genomes can be sorted by recall
(default) or precision

positional arguments:
  file                  File containing precision and recall for each genome

optional arguments:
  -h, --help            show this help message and exit
  -s {recall,precision}, --sort_by {recall,precision}
                        Sort by either precision or recall (default: recall)
  -o OUT_FILE, --out_file OUT_FILE
                        Path to store image (default: only show image)
~~~
**Example:**
~~~BASH
./precision_recall_per_genome.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2 | \
./plot_by_genome.py
~~~
**Output:**
Figure is shown on screen.

## utils/convert_fasta_bins_to_biobox_format.py
~~~BASH
usage: convert_fasta_bins_to_biobox_format.py [-h] [-o OUTPUT_FILE]
                                              paths [paths ...]

Convert bins in FASTA files to CAMI tsv format

positional arguments:
  paths                 FASTA files including full paths

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file
~~~
**Example:**
~~~BASH
./utils/convert_fasta_bins_to_cami.py \
/path/to/file/maxbin.out.001.fasta \
/path/to/file/maxbin.out.002.fasta \
/path/to/file/maxbin.out.003.fasta \
/path/to/file/maxbin.out.004.fasta \
/path/to/file/maxbin.out.005.fasta \
-o bins.tsv
~~~
Alternatively:
~~~BASH
./utils/convert_fasta_bins_to_cami.py /path/to/file/maxbin.out.0* -o bins.tsv
~~~
**Output:**
File bins.tsv is created in the working directory.

## utils/add_length_column.py
~~~BASH
usage: add_length_column.py [-h] -g GOLD_STANDARD_FILE -f FASTA_FILE

Add length column _LENGTH to gold standard mapping and print mapping on the
standard output

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        FASTA or FASTQ file with sequences of gold standard
~~~
**Example:**
~~~BASH
./utils/add_length_column.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz
~~~
**Output:**
~~~BASH
@Version:0.9.1
@SampleID:gsa

@@SEQUENCEID    BINID           TAXID  _contig_id                                    _number_reads _LENGTH
RL|S1|C10817    Sample18_57     45202  Sample18_57_from_2_to_20519_total_20518       44394         20518
RL|S1|C11497    Sample22_57     10239  Sample22_57_from_4_to_37675_total_37672       18432         37672
RL|S1|C6571     evo_1286_AP.033 1385   contig_1_4_from_3_to_69916_total_69914        30978         69914
RL|S1|C10560    evo_1286_AP.033 1385   contig_1_4_from_69981_to_1065637_total_995657 443334        995657
...
~~~
Note that only columns SEQUENCEID and BINID are required in a gold standard mapping file. The added
optional column _LENGTH, however, eliminates the need for a FASTA or FASTQ file when
evaluating binnings.

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
