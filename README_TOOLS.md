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

### src/utils/add_length_column.py
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
Column _LENGTH in the gold standard mapping eliminates the need for a FASTA or FASTQ file in program amber.py.

### src/create_summary_pdf.py
_create_summary_pdf.py_ must be run after tool _amber.py_. The input directory of _create_summary_pdf.py_ must be the output directory of _amber.py_.
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
python3 src/create_summary_pdf.py -i input_dir -c 1,2,3,4 
~~~
**Output:**
File summary.pdf will be created in directory input_dir (or in output_dir, if provided).

Providing a gold standard binning file is optional. It allows to insert a row with information of the gold standard in the table of contamination and completeness.