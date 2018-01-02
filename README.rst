|CircleCI|

Introduction
============

AMBER (Assessment of Metagenome BinnERs) is an evaluation package for
the comparative assessment of genome reconstructions from metagenome
benchmark datasets. It provides performance metrics, results rankings,
and comparative visualizations for assessing multiple programs or
parameter effects. The provided metrics were used in the first community
benchmarking challenge of the initiative for the `Critical Assessment of
Metagenomic Interpretation <http://www.cami-challenge.org/>`__.

AMBER produces print-ready and interactive plots. See a produced example
page at *https://cami-challenge.github.io/AMBER/*

Requirements
============

See `default.txt <requirements/default.txt>`__ for all dependencies.

Optional:

-  tox, for automatic tests
-  LaTeX, for combining plots into a PDF file with tool
   `*src/create\_summary\_pdf.py* <README_TOOLS.md#srccreate_summary_pdfpy>`__

Installation
------------

Install pip first (tested on Linux Ubuntu 16.04):

.. code:: bash

    sudo apt install python3-pip

Then run:

.. code:: bash

    pip3 install https://github.com/CAMI-challenge/AMBER/archive/v0.6.3.tar.gz 

Make sure to add AMBER to your PATH:

.. code:: bash

    echo 'PATH=$PATH:${HOME}/.local/bin' >> ~/.bashrc
    source ~/.bashrc

You can also run `AMBER as a Biobox <#run-amber-as-a-biobox>`__.

User Guide
==========

Input
-----

As input, AMBER's main tool `*amber.py* <#evaluatepy>`__ uses three
files: 1. A gold standard mapping of contigs or read IDs to genomes in
the `CAMI binning Bioboxes
format <https://github.com/bioboxes/rfc/tree/master/data-format>`__.
Columns are tab separated. Example: :sub:`~`\ BASH @Version:0.9.1
@SampleID:gsa

@@SEQUENCEID BINID \_LENGTH RH\|P\|C37126 Sample6\_89 25096 RH\|P\|C3274
Sample9\_91 10009 RH\|P\|C26099 1053046 689201 RH\|P\|C35075 1053046
173282 RH\|P\|C20873 1053046 339258 :sub:`~` See
`here <./test/gsa_mapping.binning>`__ another example. Note: column
\_LENGTH is optional, but eliminates the need for a FASTA or FASTQ file
(input 3 below).

2. One or more files with bin assignments for the sequences also in the
   `CAMI binning Bioboxes
   format <https://github.com/bioboxes/rfc/tree/master/data-format>`__,
   with each file containing all the bin assignments from a binning
   program. A tool for converting FASTA files, such that each file
   represents a bin, is available (see
   `*src/utils/convert\_fasta\_bins\_to\_biobox\_format.py* <README_TOOLS.md#srcutilsconvert_fasta_bins_to_biobox_formatpy>`__).
3. A FASTA or FASTQ file with the sequences for obtaining their lengths.
   Optionally, the lenghts may be added to the gold standard mapping
   file at column \_LENGTH using tool
   `*src/utils/add\_length\_column.py* <README_TOOLS.md#srcutilsadd_length_columnpy>`__.
   In this way, *amber.py* no longer requires a FASTA or FASTQ file.

Additional parameters may be specified - see below.

List of metrics and abbreviations
---------------------------------

-  **avg\_purity**: purity averaged over genome bins
-  **std\_dev\_purity**: standard deviation of purity averaged over
   genome bins
-  **sem\_purity**: standard error of the mean of purity averaged over
   genome bins
-  **avg\_completeness**: completeness averaged over genome bins
-  **std\_dev\_completeness**: standard deviation of completeness
   averaged over genome bins
-  **sem\_completeness**: standard error of the mean of completeness
   averaged over genome bins
-  **avg\_purity\_per\_bp**: average purity per base pair
-  **avg\_completeness\_per\_bp**: average completeness per base pair
-  **rand\_index\_by\_bp**: Rand index weighed by base pairs
-  **rand\_index\_by\_seq**: Rand index weighed by sequence counts
-  **a\_rand\_index\_by\_bp**: adjusted Rand index weighed by base pairs
-  **a\_rand\_index\_by\_seq**: adjusted Rand index weighed by sequence
   counts
-  **percent\_assigned\_bps**: percentage of base pairs that were
   assigned to bins
-  **accuracy**: accuracy
-  **>0.5compl<0.1cont**: number of bins with more than 50% completeness
   and less than 10% contamination
-  **>0.7compl<0.1cont**: number of bins with more than 70% completeness
   and less than 10% contamination
-  **>0.9compl<0.1cont**: number of bins with more than 90% completeness
   and less than 10% contamination
-  **>0.5compl<0.05cont**: number of bins with more than 50%
   completeness and less than 5% contamination
-  **>0.7compl<0.05cont**: number of bins with more than 70%
   completeness and less than 5% contamination
-  **>0.9compl<0.05cont**: number of bins with more than 90%
   completeness and less than 5% contamination

Tools
-----

amber.py
~~~~~~~~

.. code:: bash

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
      -m, --map_by_completeness   Map genomes to bins by maximizing completeness

**Example:** :sub:`[STRIKEOUT:BASH ./amber.py -g
test/gsa\_mapping.binning
-l "MaxBin 2.0, CONCOCT, MetaBAT"
-p 1
-r test/unique\_common.tsv
-k "circular element"
test/naughty\_carson\_2
test/goofy\_hypatia\_2
test/elated\_franklin\_0
-o output\_dir/]` **Output:** :sub:`[STRIKEOUT:BASH tool avg\_purity
std\_dev\_purity sem\_purity avg\_completeness std\_dev\_completeness
sem\_completeness avg\_purity\_per\_bp avg\_completeness\_per\_bp
rand\_index\_by\_bp rand\_index\_by\_seq a\_rand\_index\_by\_bp
a\_rand\_index\_by\_seq percent\_assigned\_bps accuracy
>0.5compl<0.1cont >0.7compl<0.1cont >0.9compl<0.1cont >0.5compl<0.05cont
>0.7compl<0.05cont >0.9compl<0.05cont MaxBin 2.0 0.948 0.095 0.016 0.799
0.364 0.058 0.934 0.838 0.995 0.951 0.917 0.782 0.864 0.807 28 28 24 23
23 21 CONCOCT 0.837 0.266 0.052 0.517 0.476 0.069 0.684 0.936 0.972
0.946 0.644 0.751 0.967 0.661 18 17 15 16 16 14 MetaBAT 0.822 0.256
0.047 0.57 0.428 0.065 0.724 0.825 0.976 0.965 0.674 0.860 0.917 0.664
17 16 12 17 16 12]` Directory *output\_dir* will contain: \*
**summary.tsv**: contains the same table as the output above with
tab-separated values \* **summary.html**: HTML page with results summary
and interactive graphs \* **avg\_purity\_completeness.png + .pdf**:
figure of average purity vs. average completeness \*
**avg\_purity\_completeness\_per\_bp.png + .pdf**: figure of purity vs.
completeness per base pair \* **ari\_vs\_assigned\_bps.png + .pdf**:
figure of adjusted Rand index weighed by number of base pairs vs.
percentage of assigned base pairs \* **rankings.txt**: tools sorted by
average purity, average completeness, and sum of average purity and
average completeness

In the same directory, subdirectories *naughty\_carson\_2*,
*goofy\_hypatia\_2*, and *elated\_franklin\_0* will be created with the
following files: \* **rand\_index.tsv**: contains value of (adjusted)
Rand index and percentage of assigned/binned bases. Rand index is both
weighed and unweighed by base pairs \* **purity\_completeness.tsv**:
contains purity and completeness per genome bin \*
**purity\_completeness\_avg.tsv**: contains purity and completeness
averaged over genome bins. Includes standard deviation and standard
error of the mean \* **purity\_completeness\_by\_bpcount.tsv**: contains
purity and completeness weighed by base pairs \* **heatmap.png + .pdf**:
heatmap representing base pair assignments to predicted bins vs. their
true origins from the underlying genomes

For a complete list of tools, see `README\_TOOLS.md <./README_TOOLS.md>`__.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run AMBER as a Biobox
---------------------

Run the AMBER docker image with the command:

.. code:: bash

    docker run -v $(pwd)/input/gold_standard.fasta:/bbx/input/gold_standard.fasta -v $(pwd)/input/gsa_mapping.binning:/bbx/input/gsa_mapping.binning  -v  $(pwd)/input/test_query.binning:/bbx/input/test_query.binning  -v  $(pwd)/output:/bbx/output -v $(pwd)/input/biobox.yaml:/bbx/input/biobox.yaml cami/amber:latest default

where biobox.yaml contains the following:

.. code:: yaml

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

Developer Guide
===============

We are using `tox <https://tox.readthedocs.io/en/latest/>`__ for project
automation.

Tests
~~~~~

If you want to run tests, just type *tox* in the project's root
directory:

.. code:: bash

    tox

You can use all libraries that AMBER depends on by activating tox's
virtual environment with the command:

.. code:: bash

    source  <project_directory>/.tox/py35/bin/activate

Update GitHub page
~~~~~~~~~~~~~~~~~~

In order to update *https://cami-challenge.github.io/AMBER*, modify file
index.html.

Make a Release
~~~~~~~~~~~~~~

If the dev branch is merged into the master branch:

1. Update `version.py <version.py>`__ according to `semantic
   versioning <semver.org>`__ on the dev branch.

2. Merge the dev branch into the master branch.

3. Make a release on GitHub with the same version number provided in
   `version.py <version.py>`__ .

The tool can now be installed by using the master branch

.. code:: bash

    pip3 install https://github.com/CAMI-challenge/AMBER/archive/master.zip

or a specific tag

.. code:: bash

    pip3 install https://github.com/CAMI-challenge/AMBER/archive/tag.tar.gz 

.. |CircleCI| image:: https://circleci.com/gh/CAMI-challenge/AMBER/tree/master.svg?style=svg
   :target: https://circleci.com/gh/CAMI-challenge/AMBER/tree/master
