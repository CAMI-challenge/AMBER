[![CircleCI](https://circleci.com/gh/CAMI-challenge/AMBER/tree/master.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/AMBER/tree/master)

# Introduction
AMBER (Assessment of Metagenome BinnERs) is an evaluation package for the comparative assessment of genome reconstructions from metagenome benchmark datasets. It provides performance metrics, results rankings, and comparative visualizations for assessing multiple programs or parameter effects. The provided metrics were used in the first community benchmarking challenge of the initiative for the [Critical Assessment of Metagenomic Interpretation](http://www.cami-challenge.org/).

AMBER produces among other files interactive plots. 
Example: **https://cami-challenge.github.io/AMBER/**

# Requirements

* python &ge; 3.5
* numpy &ge; 1.13.0
* biopython &ge; 1.69.0
* matplotlib &ge; 2.0.2
* bokeh &ge; 0.12.6
* pandas &ge; 0.20.3

Optional:

* tox, for automatic tests
* LaTeX, for combining plots into a PDF file with tool [_create_summary_pdf.py_](README_TOOLS.md#create_summary_pdfpy)

## Installation

Install AMBER by installing pip first 

For Example on Ubuntu 16.04:

~~~BASH
sudo apt install python3-pip
~~~

and running then 

~~~BASH
pip3 install  https://github.com/CAMI-challenge/AMBER/archive/tag.tar.gz 
~~~

where **tag** is the release number you can find on the [releases page](https://github.com/CAMI-challenge/AMBER/releases)

You can also run [AMBER as a Biobox](#run-amber-as-a-biobox). 

# User Guide

## Input
As input, AMBER's main tool [_evaluate.py_](#evaluatepy) uses three files:
1. A gold standard mapping of contigs or read IDs to genomes in the [CAMI binning Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format). Columns are tab separated. Example:
~~~BASH
@Version:0.9.1
@SampleID:gsa

@@SEQUENCEID BINID      _LENGTH
RH|P|C37126  Sample6_89 25096
RH|P|C3274   Sample9_91 10009
RH|P|C26099  1053046    689201
RH|P|C35075  1053046    173282
RH|P|C20873  1053046    339258
~~~
See [here](./test/gsa_mapping.binning) another example. Note: column _LENGTH is optional, but eliminates the need for a FASTA or FASTQ file (input 3 below).

2. One or more files with bin assignments for the sequences also in the [CAMI binning Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format), with each file containing all the bin assignments from a binning program. A tool for converting FASTA files, such that each file represents a bin, is available (see [_utils/convert_fasta_bins_to_biobox_format.py_](README_TOOLS.md#utilsconvert_fasta_bins_to_biobox_formatpy)).
3. A FASTA or FASTQ file with the sequences for obtaining their lengths. Optionally, the lenghts may be added to the gold standard mapping file at column _LENGTH using tool [_utils/add_length_column.py_](README_TOOLS.md#utilsadd_length_columnpy). In this way, _evaluate.py_ no longer requires a FASTA or FASTQ file.

Additional parameters may be specified - see below.

## List of metrics and abbreviations

* **avg_precision**: precision averaged over genome bins
* **std_dev_precision**: standard deviation of precision averaged over genome bins
* **sem_precision**: standard error of the mean of precision averaged over genome bins
* **avg_recall**: recall averaged over genome bins
* **std_dev_recall**: standard deviation of recall averaged over genome bins
* **sem_recall**: standard error of the mean of recall averaged over genome bins
* **precision**: average precision per base pair
* **recall**: average recall per base pair
* **rand_index_by_bp**: Rand index weighed by base pairs
* **rand_index_by_seq**: Rand index weighed by sequence counts
* **a_rand_index_by_bp**: adjusted Rand index weighed by base pairs
* **a_rand_index_by_seq**: adjusted Rand index weighed by sequence counts
* **percent_assigned_bps**: percentage of base pairs that were assigned to bins
* **accuracy**: accuracy
* **\>0.5compl<0.1cont**: number of bins with more than 50% completeness and less than 10% contamination
* **\>0.7compl<0.1cont**: number of bins with more than 70% completeness and less than 10% contamination
* **\>0.9compl<0.1cont**: number of bins with more than 90% completeness and less than 10% contamination
* **\>0.5compl<0.05cont**: number of bins with more than 50% completeness and less than 5% contamination
* **\>0.7compl<0.05cont**: number of bins with more than 70% completeness and less than 5% contamination
* **\>0.9compl<0.05cont**: number of bins with more than 90% completeness and less than 5% contamination

## Tools

### evaluate.py
~~~BASH
usage: amber.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE] [-l LABELS]
                   [-p FILTER] [-r REMOVE_GENOMES] [-k KEYWORD] -o OUTPUT_DIR
                   [-m]
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
  -r REMOVE_GENOMES, --remove_genomes REMOVE_GENOMES
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in the second column of file with list of
                        genomes to be removed (no keyword=remove all genomes
                        in list)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write the results to
  -m, --map_by_recall   Map genomes to bins by maximizing recall
~~~
**Example:**
~~~BASH
./amber.py -g test/gsa_mapping.binning \
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
tool       avg_precision std_dev_precision sem_precision avg_recall std_dev_recall sem_recall precision recall rand_index_by_bp rand_index_by_seq a_rand_index_by_bp a_rand_index_by_seq percent_assigned_bps accuracy >0.5compl<0.1cont >0.7compl<0.1cont >0.9compl<0.1cont >0.5compl<0.05cont >0.7compl<0.05cont >0.9compl<0.05cont
MaxBin 2.0 0.948         0.095             0.016         0.799      0.364          0.058      0.934     0.838  0.995            0.951             0.917              0.782               0.864                0.807    28                28                24                23                 23                 21
CONCOCT    0.837         0.266             0.052         0.517      0.476          0.069      0.684     0.936  0.972            0.946             0.644              0.751               0.967                0.661    18                17                15                16                 16                 14
MetaBAT    0.822         0.256             0.047         0.57       0.428          0.065      0.724     0.825  0.976            0.965             0.674              0.860               0.917                0.664    17                16                12                17                 16                 12
~~~
Directory _output_dir_ will contain:
* **summary.html**: HTML page with results summary and interactive graphs
* **avg_precision_recall.png + .pdf**: figure of average precision vs. average recall
* **ari_vs_assigned_bps.png + .pdf**: figure of adjusted Rand index weighed by number of base pairs vs. percentage of assigned base pairs
* **rankings.txt**: tools sorted by average precision, average recall, and sum of average precision and average recall

In the same directory, subdirectories _naughty_carson_2_, _goofy_hypatia_2_, and _elated_franklin_0_ will be created with the following files:
* **rand_index.tsv**: contains value of (adjusted) Rand index and percentage of assigned/binned bases. Rand index is both weighed and unweighed by base pairs
* **precision_recall.tsv**: contains precision and recall per genome bin
* **precision_recall_avg.tsv**: contains precision and recall averaged over genome bins. Includes standard deviation and standard error of the mean
* **precision_recall_by_bpcount.tsv**: contains precision and recall weighed by base pairs
<!---* **genomes_sorted_by_precision.png + .pdf**: figure of precision and recall per genome with genomes sorted by precision-->
<!---* **genomes_sorted_by_recall.png + .pdf**: figure of precision and recall per genome with genomes sorted by recall-->

### For a complete list of tools, see [README_TOOLS.md](./README_TOOLS.md).

## Run AMBER as a Biobox


Run the amber docker image by specifying the executing the following command

~~~BASH
docker run -v $(pwd)/input/gold_standard.fasta:/bbx/input/gold_standard.fasta -v $(pwd)/input/gsa_mapping.binning:/bbx/input/gsa_mapping.binning  -v  $(pwd)/input/test_query.binning:/bbx/input/test_query.binning  -v  $(pwd)/output:/bbx/output -v $(pwd)/input/biobox.yaml:/bbx/input/biobox.yaml cami/amber:latest default
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

We are using [tox](https://tox.readthedocs.io/en/latest/) for project automation.

### Tests

If you want to run tests, just type _tox_ in the project's root directory:

~~~BASH
tox
~~~

By running tox you can use all libraries that amber depends on by running 

~~~BASH
source  <project_directory>/.tox/py35/bin/activate
~~~

### Update GitHub page

In order to update **https://cami-challenge.github.io/AMBER**
modify the file index.html.

### Make a Release

If the dev branch is merged into the master branch

1. Update the [version.py](version.py) according to (semantic versioning)[semver.org] on the dev branch.

2. Merge the dev branch into the master branch.

3. Make a release on GitHub with the same version number provided in [version.py](version.py) 

The tool can now be installed by using the master branch 

~~~BASH
pip3 install https://github.com/CAMI-challenge/AMBER/archive/master.zip
~~~

or a specific tag

~~~BASH
pip3 install https://github.com/CAMI-challenge/AMBER/archive/tag.tar.gz 
~~~
