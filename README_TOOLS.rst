src/precision\_recall.py
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    usage: precision_recall.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE]
                               [-l LABELS] [-p FILTER] [-r REMOVE_GENOMES]
                               [-k KEYWORD] [-m]
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
      -r REMOVE_GENOMES, --remove_genomes REMOVE_GENOMES
                            File with list of genomes to be removed
      -k KEYWORD, --keyword KEYWORD
                            Keyword in the second column of file with list of
                            genomes to be removed (no keyword=remove all genomes
                            in list)
      -m, --map_by_recall   Map genomes to bins by maximizing recall

| **Example:** :sub:`~`\ BASH python3 src/precision\_recall.py -g
  test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| -r test/unique\_common.tsv -k "circular element"
| -p 1
| -l "MaxBin 2.0, CONCOCT, MetaBAT"
| test/naughty\_carson\_2 test/goofy\_hypatia\_2
  test/elated\_franklin\_0 :sub:`:sub:`:sub:` **Output:** ```\ BASH tool
  precision std\_dev\_precision sem\_precision recall std\_dev\_recall
  sem\_recall MaxBin 2.0 0.948 0.095 0.016 0.799 0.364 0.058 CONCOCT
  0.837 0.266 0.052 0.517 0.476 0.069 MetaBAT 0.822 0.256 0.047 0.570
  0.428 0.065 :sub:`~`

src/precision\_recall\_per\_bin.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    usage: precision_recall_per_bin.py [-h] -g GOLD_STANDARD_FILE [-f FASTA_FILE]
                                       [-m]
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
      -m, --map_by_recall   Map genomes to bins by maximizing recall

| **Example:** :sub:`~`\ BASH python3 precision\_recall\_per\_bin.py -g
  test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| test/naughty\_carson\_2 :sub:`:sub:`:sub:` **Output:** ```\ BASH
  genome precision recall predicted\_size correctly\_predicted
  real\_size 1049005 0.998458014811 1.0 4114177 4107833 4107833
  evo\_1035930.032 1.0 1.0 2294831 2294831 2294831 Sample18\_8
  0.409809009246 1.0 121786 49909 49909 evo\_1035930.029 0.995223021915
  1.0 2423708 2412130 2412130 ... :sub:`~`

src/utils/exclude\_genomes.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    usage: exclude_genomes.py [-h] -r REMOVE_GENOMES [-k KEYWORD] [file]

    Exclude genome bins from table of precision and recall per genome. The table
    can be provided as file or via the standard input

    positional arguments:
      file                  File containing precision and recall for each genome

    optional arguments:
      -h, --help            show this help message and exit
      -r REMOVE_GENOMES, --remove_genomes REMOVE_GENOMES
                            File with list of genomes to be removed
      -k KEYWORD, --keyword KEYWORD
                            Keyword in the second column of file with list of
                            genomes to be removed (no keyword=remove all genomes
                            in list)

**Example:**

| The example computes the table of precision and recall and pipes it to
  utils/exclude\_genomes.py. :sub:`~`\ BASH python3
  src/precision\_recall\_per\_bin.py -g test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| test/naughty\_carson\_2 \|
| python3 src/utils/exclude\_genomes.py -r test/unique\_common.tsv -k
  "circular element" :sub:`~`

**Output:**

The output the is the table from precision\_recall\_per\_bin.py without
the excluded genomes.

src/precision\_recall\_average.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

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

| **Example:** :sub:`~`\ BASH python3 src/precision\_recall\_per\_bin.py
  -g test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| test/naughty\_carson\_2 \|
| python3 src/utils/exclude\_genomes.py -r test/unique\_common.tsv -k
  "circular element" \|
| python3 src/precision\_recall\_average.py -p 1 -l "MaxBin 2.0"
  :sub:`:sub:`:sub:` **Output:** ```\ BASH tool precision
  std\_dev\_precision sem\_precision recall std\_dev\_recall sem\_recall
  MaxBin 2.0 0.948 0.095 0.016 0.799 0.364 0.058 :sub:`~`

src/precision\_recall\_by\_bpcount.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

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

| **Example:** :sub:`~`\ BASH python3
  src/precision\_recall\_by\_bpcount.py -g test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| test/naughty\_carson\_2 :sub:`:sub:`:sub:` **Output:** ```\ BASH
  precision recall 0.934 0.838 :sub:`~`

src/rand\_index.py
~~~~~~~~~~~~~~~~~~

.. code:: bash

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

| **Example:** :sub:`~`\ BASH python3 src/rand\_index.py -g
  test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| test/naughty\_carson\_2 :sub:`:sub:`:sub:` **Output:** ```\ BASH
  rand\_index\_by\_bp rand\_index\_by\_seq a\_rand\_index\_by\_bp
  a\_rand\_index\_by\_seq percent\_assigned\_bps 0.995 0.951 0.917 0.782
  0.864 :sub:`~`

src/genome\_recovery.py
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    usage: genome_recovery.py [-h] [-p FILTER] [-l LABEL] [-x MIN_COMPLETENESS]
                              [-y MAX_CONTAMINATION]
                              [file]

    Calculate number of genome bins recovered with more than the specified
    thresholds of completeness and contamination. Default: >50%, >70%, >90%
    completeness vs. <10%, <5% contamination

    positional arguments:
      file                  File containing precision and recall for each genome

    optional arguments:
      -h, --help            show this help message and exit
      -p FILTER, --filter FILTER
                            Filter out [FILTER]% smallest bins (default: 0)
      -l LABEL, --label LABEL
                            Binning name
      -x MIN_COMPLETENESS, --min_completeness MIN_COMPLETENESS
                            Comma-separated list of min. completeness thresholds
                            (default: 50,70,90)
      -y MAX_CONTAMINATION, --max_contamination MAX_CONTAMINATION
                            Comma-separated list of max. contamination thresholds
                            (default: 10,5)

| **Example:** :sub:`~`\ BASH python3 src/precision\_recall\_per\_bin.py
  -g test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| test/naughty\_carson\_2 \|
| python3 src/genome\_recovery.py -l "MaxBin 2.0" -p 1
  :sub:`:sub:`:sub:` **Output:** ```\ BASH MaxBin 2.0 >50% complete >70%
  complete >90% complete <10% contamination 28 28 24 <5% contamination
  23 23 21 :sub:`~`

src/plot\_by\_genome.py
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

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

| **Example:** :sub:`~`\ BASH python3 src/precision\_recall\_per\_bin.py
  -g test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
| test/naughty\_carson\_2 \|
| python3 src/plot\_by\_genome.py :sub:`~` **Output:** Figure is shown
  on screen.

src/utils/convert\_fasta\_bins\_to\_biobox\_format.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    usage: convert_fasta_bins_to_biobox_format.py [-h] [-o OUTPUT_FILE]
                                                  paths [paths ...]

    Convert bins in FASTA files to CAMI tsv format

    positional arguments:
      paths                 FASTA files including full paths

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_FILE, --output_file OUTPUT_FILE
                            Output file

**Example:** :sub:`[STRIKEOUT:BASH python3
src/utils/convert\_fasta\_bins\_to\_cami.py
/path/to/file/maxbin.out.001.fasta
/path/to/file/maxbin.out.002.fasta
/path/to/file/maxbin.out.003.fasta
/path/to/file/maxbin.out.004.fasta
/path/to/file/maxbin.out.005.fasta
-o bins.tsv]` Alternatively: :sub:`[STRIKEOUT:BASH python3
src/utils/convert\_fasta\_bins\_to\_cami.py /path/to/file/maxbin.out.0\*
-o bins.tsv]` **Output:** File bins.tsv is created in the working
directory.

src/utils/add\_length\_column.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    usage: add_length_column.py [-h] -g GOLD_STANDARD_FILE -f FASTA_FILE

    Add length column _LENGTH to gold standard mapping and print mapping on the
    standard output

    optional arguments:
      -h, --help            show this help message and exit
      -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                            Gold standard - ground truth - file
      -f FASTA_FILE, --fasta_file FASTA_FILE
                            FASTA or FASTQ file with sequences of gold standard

| **Example:** File
  CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz used
  in the example can be downloaded
  `here <https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz>`__.
  :sub:`~`\ BASH python3 src/utils/add\_length\_column.py -g
  test/gsa\_mapping.binning
| -f
  test/CAMI\_low\_RL\_S001\_\_insert\_270\_GoldStandardAssembly.fasta.gz
  :sub:`:sub:`:sub:` **Output:** ```\ BASH @Version:0.9.1 @SampleID:gsa

@@SEQUENCEID BINID \_LENGTH RL\|S1\|C10817 Sample18\_57 20518
RL\|S1\|C11497 Sample22\_57 37672 RL\|S1\|C6571 evo\_1286\_AP.033 69914
RL\|S1\|C10560 evo\_1286\_AP.033 995657 ... :sub:`~` Column \_LENGTH in
the gold standard mapping eliminates the need for a FASTA or FASTQ file
in program evaluate.py.

src/create\_summary\_pdf.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

*create\_summary\_pdf.py* must be run after tool *evaluate.py*. The
input directory of *create\_summary\_pdf.py* must be the output
directory of *evaluate.py*. :sub:`~`\ BASH usage:
create\_summary\_pdf.py [-h] [-c FIGURE\_CODES] -i INPUT\_DIR [-o
OUTPUT\_DIR] [-g GOLD\_STANDARD\_FILE] [-r REMOVE\_GENOMES] [-k KEYWORD]

Combine figures and table of completeness and contamination into file
summary.pdf

optional arguments: -h, --help show this help message and exit -c
FIGURE\_CODES, --figure\_codes FIGURE\_CODES Comma-separated figure
codes without spaces (1=average precision vs. average recall, 2=ari vs.
%assigned bps, 3=weighed precision vs. weighed recall, 4=table of
contamination and completeness, 5=precision vs. recall per bin; default
1,2,3,4) -i INPUT\_DIR, --input\_dir INPUT\_DIR Directory to read file
summary.tsv and figures from -o OUTPUT\_DIR, --output\_dir OUTPUT\_DIR
Directory to write pdf and temporary files (if different from
input\_dir) -g GOLD\_STANDARD\_FILE, --gold\_standard\_file
GOLD\_STANDARD\_FILE Gold standard - ground truth - file -r
REMOVE\_GENOMES, --remove\_genomes REMOVE\_GENOMES File with list of
genomes to be removed -k KEYWORD, --keyword KEYWORD Keyword in the
second column of file with list of genomes to be removed (no
keyword=remove all genomes in list) :sub:`:sub:`:sub:` **Example:**
```\ BASH python3 src/create\_summary\_pdf.py -i input\_dir -c 1,2,3,4
:sub:`~` **Output:** File summary.pdf will be created in directory
input\_dir (or in output\_dir, if provided).

Providing a gold standard binning file is optional. It allows to insert
a row with information of the gold standard in the table of
contamination and completeness.
