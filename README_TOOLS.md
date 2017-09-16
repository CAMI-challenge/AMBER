### precision_recall.py
~~~BASH
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

### precision_recall_per_bin.py
~~~BASH
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
~~~
**Example:**
~~~BASH
./precision_recall_per_bin.py -g test/gsa_mapping.binning \
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

### utils/exclude_genomes.py
~~~BASH
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
~~~
**Example:**

The example computes the table of precision and recall and pipes it to utils/exclude_genomes.py.
~~~BASH
./precision_recall_per_bin.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2 | \
./utils/exclude_genomes.py -r test/unique_common.tsv -k "circular element"
~~~

**Output:**

The output the is the table from precision_recall_per_bin.py without the excluded genomes.

### precision_recall_average.py
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
./precision_recall_per_bin.py -g test/gsa_mapping.binning \
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

### precision_recall_by_bpcount.py
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

### rand_index.py
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

### genome_recovery.py
~~~BASH
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
~~~
**Example:**
~~~BASH
./precision_recall_per_bin.py -g test/gsa_mapping.binning \
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

### plot_by_genome.py
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
./precision_recall_per_bin.py -g test/gsa_mapping.binning \
-f test/CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta.gz \
test/naughty_carson_2 | \
./plot_by_genome.py
~~~
**Output:**
Figure is shown on screen.

### utils/convert_fasta_bins_to_biobox_format.py
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

### utils/add_length_column.py
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
./utils/add_length_column.py -g test/gsa_mapping.binning \
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
Column _LENGTH in the gold standard mapping eliminates the need for a FASTA or FASTQ file in program evaluate.py.

### create_summary_pdf.py
_create_summary_pdf.py_ must be run after tool _evaluate.py_. The input directory of _create_summary_pdf.py_ must be the output directory of _evaluate.py_.
~~~BASH
usage: create_summary_pdf.py [-h] [-c FIGURE_CODES] -i INPUT_DIR
                             [-o OUTPUT_DIR] [-g GOLD_STANDARD_FILE]
                             [-r REMOVE_GENOMES] [-k KEYWORD]

Combine figures and table of completeness and contamination into file
summary.pdf

optional arguments:
  -h, --help            show this help message and exit
  -c FIGURE_CODES, --figure_codes FIGURE_CODES
                        Comma-separated figure codes without spaces (1=average
                        precision vs. average recall, 2=ari vs. %assigned bps,
                        3=weighed precision vs. weighed recall, 4=table of
                        contamination and completeness, 5=precision vs. recall
                        per bin; default 1,2,3,4)
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Directory to read file summary.tsv and figures from
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write pdf and temporary files (if
                        different from input_dir)
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard - ground truth - file
  -r REMOVE_GENOMES, --remove_genomes REMOVE_GENOMES
                        File with list of genomes to be removed
  -k KEYWORD, --keyword KEYWORD
                        Keyword in the second column of file with list of
                        genomes to be removed (no keyword=remove all genomes
                        in list)
~~~
**Example:**
~~~BASH
./create_summary_pdf.py -i input_dir -c 1,2,3,4 
~~~
**Output:**
File summary.pdf will be created in directory input_dir (or in output_dir, if provided).

Providing a gold standard binning file is optional. It allows to insert a row with information of the gold standard in the table of contamination and completeness.