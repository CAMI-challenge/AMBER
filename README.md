[![CircleCI](https://circleci.com/gh/CAMI-challenge/AMBER/tree/master.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/AMBER/tree/master)

# AMBER: Assessment of Metagenome BinnERs
AMBER is an evaluation package for the comparative assessment of genome reconstructions and taxonomic assignments from metagenome benchmark datasets. It provides performance metrics, results rankings, and comparative visualizations for assessing multiple programs or parameter effects. The provided metrics were used in the first community benchmarking challenge of the initiative for the [Critical Assessment of Metagenomic Interpretation](http://www.cami-challenge.org/).

**Metrics computed per bin**

* Predicted bin size in bps and sequences
* True positives
* (Average) Purity
* (Average) Completeness

**Metrics computed per sample**

* Accuracy
* Misclassification rate (contamination)
* Purity
* Completeness
* (Adjusted) Rand index
* Percentage of binned base pairs and sequences
* Number of genomes recovered within levels of completeness and contamination
* UniFrac (for taxonomic binning)


**Example pages produced by AMBER**

* [CAMI I high complexity challenge dataset](https://cami-challenge.github.io/AMBER/cami_i_hc/)
* [CAMI II mouse gut toy dataset](https://cami-challenge.github.io/AMBER/cami_ii_mouse_gut/)
* [CAMI II mouse gut toy dataset (using option --filter 1)](https://cami-challenge.github.io/AMBER/cami_ii_mouse_gut_filter1/)

# Installation

## Requirements

AMBER 2.0.4 has been tested with Python 3.11.

See [requirements.txt](requirements.txt) for all dependencies.

## Installation options

There are several options to install AMBER:

* [Bioconda](#bioconda)
* [Python pip](#python-pip)
* [Docker](#docker)

## Bioconda

Install and configure [Bioconda](https://bioconda.github.io/index.html) if not already installed. Then use the following command to create a Conda environment and install AMBER:

~~~BASH
conda create --name amber cami-amber
~~~

Activate the Conda environment with:

~~~BASH
conda activate amber
~~~

## Python pip

Install pip if not already installed (tested on Linux Ubuntu 22.04):

~~~BASH
sudo apt install python3-pip
~~~
Should you receive the message `Unable to locate package python3-pip`, enter the following commands and repeat the previous step.

~~~BASH
sudo add-apt-repository universe
sudo apt update
~~~

Then run:

~~~BASH
pip install cami-amber 
~~~

Make sure to add AMBER to your PATH:

~~~BASH
echo 'PATH=$PATH:${HOME}/.local/bin' >> ~/.bashrc
source ~/.bashrc
~~~

Alternatively, download or git-clone AMBER from GitHub. In AMBER's directory, install all requirements with the command:

~~~BASH
pip install -r requirements.txt 
~~~

# Docker

You can pull a pre-built [AMBER Docker BioContainer](https://bioconda.github.io/recipes/cami-amber/README.html) as follows:

~~~BASH
docker pull quay.io/biocontainers/cami-amber:<tag>
~~~

See [valid values for &lt;tag&gt;](https://quay.io/repository/biocontainers/cami-amber?tab=tags).

Alternatively, download or git-clone AMBER from GitHub. In AMBER's directory, build the Docker image with the command:

~~~BASH
docker build -t amber .
~~~

See bellow an example of how to [run AMBER using Docker](#running-amberpy-using-docker).


# User guide

## Input
As input, AMBER uses three files and an additional file for assessing taxonomic binning:
1. A gold standard mapping of contigs or read IDs to genomes and/or taxon IDs in the [CAMI binning Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format). Columns are tab separated. Example:
~~~BASH
@Version:0.9.1
@SampleID:gsa

@@SEQUENCEID BINID      TAXID  LENGTH
RH|P|C37126  Sample6_89 45202  25096
RH|P|C3274   Sample9_91 32644  10009
RH|P|C26099  1053046    765201 689201
RH|P|C35075  1053046    765201 173282
RH|P|C20873  1053046    765201 339258
~~~
See [here](./test/gsa_mapping.binning) another example. Observations:
* The value of the SampleID header tag must uniquely identify a sample and be the same in the gold standard and the predictions (input 2 below).
* Column BINID (TAXID) is required to assess genome (taxonomic) binning.
* Column LENGTH can be added to a mapping file using tool [_src/utils/add_length_column.py_](#srcutilsadd_length_columnpy).

2. One or more files, each containing the bin assignments from a binning program, also in the [CAMI binning Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format). Column LENGTH is not required  (LENGTH is only required in the gold standard).

Note: a tool for converting FASTA files, such that each file represents a bin, is available (see [_src/utils/convert_fasta_bins_to_biobox_format.py_](#srcutilsconvert_fasta_bins_to_biobox_formatpy)).

3. For assessing **taxonomic binning**, AMBER also requires the file **nodes.dmp** from NCBI. Download taxdump.tar.gz from [ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz), extract nodes.tmp, and provide it to AMBER with option `--ncbi_dir`.

## Input format for multiple samples

Binnings of datasets with multiple samples are supported by AMBER. For each binning program, simply concatenate the binnings of the different samples into a single file to obtain one binning file per program. The gold standard must also consist in one file for all samples. Remember: binnings for the same sample must have the same SampleID.


## Running _amber.py_

~~~BASH
usage: AMBER [-h] -g GOLD_STANDARD_FILE [-l LABELS] [-p FILTER] [-n MIN_LENGTH] -o OUTPUT_DIR [--stdout] [-d DESC] [--colors COLORS] [--silent] [-v] [-x MIN_COMPLETENESS]
             [-y MAX_CONTAMINATION] [-r REMOVE_GENOMES] [-k KEYWORD] [--genome_coverage GENOME_COVERAGE] [--ncbi_dir NCBI_DIR]
             bin_files [bin_files ...]

AMBER: Assessment of Metagenome BinnERs

positional arguments:
  bin_files             Binning files

options:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -l LABELS, --labels LABELS
                        Comma-separated binning names
  -p FILTER, --filter FILTER
                        Filter out [FILTER]% smallest genome bins (default: 0)
  -n MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum length of sequences
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write the results to
  --stdout              Print summary to stdout
  -d DESC, --desc DESC  Description for HTML page
  --silent              Silent mode
  -v, --version         show program's version number and exit

genome binning-specific arguments:
  -x MIN_COMPLETENESS, --min_completeness MIN_COMPLETENESS
                        Comma-separated list of min. completeness thresholds (default %: 50,70,90)
  -y MAX_CONTAMINATION, --max_contamination MAX_CONTAMINATION
                        Comma-separated list of max. contamination thresholds (default %: 10,5)
  -r REMOVE_GENOMES, --remove_genomes REMOVE_GENOMES
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in the second column of file with list of genomes to be removed (no keyword=remove all genomes in list)
  --genome_coverage GENOME_COVERAGE
                        genome coverages

taxonomic binning-specific arguments:
  --ncbi_dir NCBI_DIR   Directory containing the NCBI taxonomy database dump files nodes.dmp, merged.dmp, and names.dmp
~~~
**Example:**
~~~BASH
amber.py -g test/gsa_mapping.binning \
-l "MaxBin 2.0, CONCOCT, MetaBAT" \
-p 1 \
-r test/unique_common.tsv \
-k "circular element" \
test/naughty_carson_2 \
test/goofy_hypatia_2 \
test/elated_franklin_0 \
-o output_dir/
~~~

## Running _amber.py_ using Docker

_amber.py_ can be run with the `docker run` command. Example:

~~~BASH
docker run -v $(pwd):/host amber \
amber.py \
-l "CONCOCT (CAMI), MaxBin 2.0.2 (CAMI)" \
-p 1 \
-r /host/test/unique_common.tsv \
-k "circular element" \
-g /host/test/gsa_mapping.binning \
/host/test/goofy_hypatia_2 \
/host/test/naughty_carson_2 \
-o /host/output_dir
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
File CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz used in the example can be downloaded [here](https://data.cami-challenge.org/participate).
~~~BASH
python3 src/utils/add_length_column.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz
~~~
**Output:**
~~~BASH
@Version:0.9.1
@SampleID:gsa

@@SEQUENCEID BINID           LENGTH
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

# Developer guide

We are using [tox](https://tox.readthedocs.io/en/latest/) for project automation.

### Tests

If you want to run tests, just type _tox_ in the project's root directory:

~~~BASH
tox
~~~

You can use all libraries that AMBER depends on by activating tox's virtual environment with the command: 

~~~BASH
source  <project_directory>/.tox/py311/bin/activate
~~~

### Update GitHub page

In order to update *https://cami-challenge.github.io/AMBER*, modify file index.html.

### Make a release

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
* Meyer, F., Hofmann, P., Belmann, P., Garrido-Oter, R., Fritz, A., Sczyrba, A., McHardy, A.C., **AMBER: Assessment of Metagenome BinnERs**, GigaScience 7, giy069 (2018). [https://doi.org/10.1093/gigascience/giy069](https://doi.org/10.1093/gigascience/giy069)

The metrics implemented in AMBER were used and described in the CAMI manuscript, thus you may also cite:
* Sczyrba, A., Hofmann, P., Belmann, P. et al. **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** Nat Methods 14, 1063–1071 (2017). [https://doi.org/10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)

or

* Meyer, F., Fritz, A., Deng, ZL. et al. **Critical Assessment of Metagenome Interpretation: the second round of challenges.** Nat Methods 19, 429–440 (2022). [https://doi.org/10.1038/s41592-022-01431-4](https://doi.org/10.1038/s41592-022-01431-4)

# License

AMBER 2 is licensed under GPL v3.