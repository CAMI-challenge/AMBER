# Copyright 2020 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse

HELP_GOLD_STANDARD_FILE = "Gold standard - ground truth - file"
HELP_FILE = "File containing precision and recall for each genome"
HELP_QUERY_FILE = "Binning file"
HELP_QUERY_FILES = "Binning files"
HELP_FASTA_FILE = "FASTA or FASTQ file with sequences of gold standard (required if gold standard file misses column _LENGTH)"
HELP_FILTER = "Filter out [FILTER]%% smallest genome bins (default: 0)"
HELP_GENOMES_FILE = "File with list of genomes to be removed"
HELP_LABEL = "Binning name"
HELP_LABELS = "Comma-separated binning names"
HELP_KEYWORD = "Keyword in the second column of file with list of genomes to be removed (no keyword=remove all genomes in list)"
HELP_MAP_BY_RECALL = "Map genomes to bins by maximizing completeness"
HELP_THRESHOLDS_COMPLETENESS = "Comma-separated list of min. completeness thresholds (default %%: 50,70,90)"
HELP_THRESHOLDS_CONTAMINATION = "Comma-separated list of max. contamination thresholds (default %%: 10,5)"

PARSER_GS = argparse.ArgumentParser(add_help=False)
PARSER_GS.add_argument("bin_file", help=HELP_QUERY_FILE)
PARSER_GS.add_argument("-g", "--gold_standard_file", help=HELP_GOLD_STANDARD_FILE, required=True)
PARSER_GS.add_argument("-f", "--fasta_file", help=HELP_FASTA_FILE)

PARSER_MULTI = argparse.ArgumentParser(add_help=False)
PARSER_MULTI.add_argument('file', nargs='?', type=argparse.FileType('r'), help=HELP_FILE)
PARSER_MULTI.add_argument('-p', '--filter', help=HELP_FILTER)
PARSER_MULTI.add_argument('-l', '--label', help=HELP_LABEL, required=False)

PARSER_MULTI2 = argparse.ArgumentParser(add_help=False)
PARSER_MULTI2.add_argument("bin_files", nargs='+', help=HELP_QUERY_FILES)
PARSER_MULTI2.add_argument("-g", "--gold_standard_file", help=HELP_GOLD_STANDARD_FILE, required=True)
PARSER_MULTI2.add_argument('-l', '--labels', help=HELP_LABELS, required=False)

