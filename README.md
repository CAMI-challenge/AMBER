[![CircleCI](https://circleci.com/gh/CAMI-challenge/AMBER/tree/master.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/AMBER/tree/master)

# Introduction
AMBER (Assessment of Metagenome BinnERs) is an evaluation package for the comparative assessment of genome reconstructions from metagenome benchmark datasets. It provides performance metrics, results rankings, and comparative visualizations for assessing multiple programs or parameter effects. The provided metrics were used in the first community benchmarking challenge of the initiative for the [Critical Assessment of Metagenomic Interpretation](http://www.cami-challenge.org/).

AMBER produces print-ready and interactive plots. See a produced example page at *https://cami-challenge.github.io/AMBER/*

# Requirements

See [default.txt](requirements/default.txt) for all dependencies.

## Installation

Install pip first (tested on Linux Ubuntu 18.04):

~~~BASH
sudo apt install python3-pip
~~~

Then run:

~~~BASH
pip3 install cami-amber 
~~~

Make sure to add AMBER to your PATH:

~~~BASH
echo 'PATH=$PATH:${HOME}/.local/bin' >> ~/.bashrc
source ~/.bashrc
~~~

You can also run [AMBER as a Biobox](#run-amber-as-a-biobox). 

# User Guide

## Input
As input, AMBER uses three files:
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

2. One or more files with bin assignments for the sequences also in the [CAMI binning Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format), with each file containing all the bin assignments from a binning program. A tool for converting FASTA files, such that each file represents a bin, is available (see [_src/utils/convert_fasta_bins_to_biobox_format.py_](README_TOOLS.md#srcutilsconvert_fasta_bins_to_biobox_formatpy)).
3. A FASTA or FASTQ file with the sequences for obtaining their lengths. Optionally, the lenghts may be added to the gold standard mapping file at column _LENGTH using tool [_src/utils/add_length_column.py_](README_TOOLS.md#srcutilsadd_length_columnpy). In this way, _amber.py_ no longer requires a FASTA or FASTQ file.

Additional parameters may be specified - see below.

## Computed metrics

* Average purity (averaged over recovered genome bins)
* Average completeness (averaged over recovered genome bins)
* Average purity per base pair
* Average completeness per base pair
* (Adjusted) Rand index by sequence and base pair counts
* Percentage of assigned base pairs
* Accuracy
* Number of genomes recovered within levels of completeness and contamination

## Running _amber.py_

~~~BASH
usage: amber.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE] [-l LABELS]
                [-p FILTER] [-r REMOVE_GENOMES] [-k KEYWORD] -o OUTPUT_DIR
                [-m] [-x MIN_COMPLETENESS] [-y MAX_CONTAMINATION]
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
  -m, --map_by_completeness
                        Map genomes to bins by maximizing completeness
  -x MIN_COMPLETENESS, --min_completeness MIN_COMPLETENESS
                        Comma-separated list of min. completeness thresholds
                        (default %: 50,70,90)
  -y MAX_CONTAMINATION, --max_contamination MAX_CONTAMINATION
                        Comma-separated list of max. contamination thresholds
                        (default %: 10,5)
~~~
**Example:**
~~~BASH
python3 amber.py -g test/gsa_mapping.binning \
-l "MaxBin 2.0, CONCOCT, MetaBAT" \
-p 1 \
-r test/unique_common.tsv \
-k "circular element" \
test/naughty_carson_2 \
test/goofy_hypatia_2 \
test/elated_franklin_0 \
-o output_dir/
~~~

## Running AMBER as a Biobox

Build and run the AMBER docker image with the commands:

~~~BASH
docker build -t cami/amber:latest .
docker run -v $(pwd)/input/gold_standard.fasta:/bbx/input/gold_standard.fasta -v $(pwd)/input/gsa_mapping.binning:/bbx/input/gsa_mapping.binning  -v  $(pwd)/input/test_query.binning:/bbx/input/test_query.binning  -v  $(pwd)/output:/bbx/output -v $(pwd)/input/biobox.yaml:/bbx/input/biobox.yaml cami/amber:latest default
~~~

where biobox.yaml contains the following:

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

# Utilities

### src/utils/add_length_column.py
Adds column _LENGTH to the gold standard mapping file, eliminating the need to provide a FASTA or FASTQ file to amber.py.

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
File CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz used in the example can be downloaded [here](https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz).
~~~BASH
python3 src/utils/add_length_column.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz
~~~
**Output:**
~~~BASH
@Version:0.9.1
@SampleID:gsa

@@SEQUENCEID BINID           _LENGTH
RL|S1|C10817 Sample18_57     20518
RL|S1|C11497 Sample22_57     37672
RL|S1|C6571  evo_1286_AP.033 69914
RL|S1|C10560 evo_1286_AP.033 995657
...
~~~

### src/utils/convert_fasta_bins_to_biobox_format.py
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
python3 src/utils/convert_fasta_bins_to_cami.py \
/path/to/file/maxbin.out.001.fasta \
/path/to/file/maxbin.out.002.fasta \
/path/to/file/maxbin.out.003.fasta \
/path/to/file/maxbin.out.004.fasta \
/path/to/file/maxbin.out.005.fasta \
-o bins.tsv
~~~
Alternatively:
~~~BASH
python3 src/utils/convert_fasta_bins_to_cami.py /path/to/file/maxbin.out.0* -o bins.tsv
~~~
**Output:**
File bins.tsv is created in the working directory.

# Developer Guide

We are using [tox](https://tox.readthedocs.io/en/latest/) for project automation.

### Tests

If you want to run tests, just type _tox_ in the project's root directory:

~~~BASH
tox
~~~

You can use all libraries that AMBER depends on by activating tox's virtual environment with the command: 

~~~BASH
source  <project_directory>/.tox/py35/bin/activate
~~~

### Update GitHub page

In order to update *https://cami-challenge.github.io/AMBER*, modify file index.html.

### Make a Release

If the dev branch is merged into the master branch:

1. Update [version.py](version.py) according to [semantic versioning](semver.org) on the dev branch.

2. Merge the dev branch into the master branch.

3. Make a release on GitHub with the same version number provided in [version.py](version.py) .

4. Create package and upload it to PyPI:

~~~BASH
python3 setup.py sdist bdist_wheel
twine upload dist/*
~~~

# Citation

Please cite AMBER as:
* Fernando Meyer, Peter Hofmann, Peter Belmann, Ruben Garrido-Oter, Adrian Fritz, Alexander Sczyrba, and Alice C. McHardy. (2018). **AMBER: Assessment of Metagenome BinnERs.** *GigaScience*, giy069. doi:[10.1093/gigascience/giy069](https://doi.org/10.1093/gigascience/giy069)

The metrics implemented in AMBER were used and described in the CAMI manuscript, thus you may also cite:
* Sczyrba, Hofmann, Belmann, et al. (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)
