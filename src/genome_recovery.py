#!/usr/bin/env python

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
import numpy as np
import collections
import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src.utils import argparse_parents
from src.utils import filter_tail
from src.utils import load_data


def get_defaults():
    return [.5, .7, .9], [.1, .05]


def calc_table(pd_bins_rank, min_completeness, max_contamination):
    if not min_completeness:
        min_completeness = get_defaults()[0]
    if not max_contamination:
        max_contamination = get_defaults()[1]

    genome_recovery = np.zeros((len(max_contamination), len(min_completeness)), dtype=int)
    for row in pd_bins_rank.itertuples():
        if np.isnan(row.precision_bp):
            continue
        contamination = 1 - row.precision_bp
        for i in range(len(max_contamination)):
            for j in range(len(min_completeness)):
                if row.recall_bp > min_completeness[j] and contamination < max_contamination[i]:
                    genome_recovery[i][j] += 1
    return genome_recovery


def calc_dict(pd_bins_rank, min_completeness, max_contamination):
    if not min_completeness:
        min_completeness = get_defaults()[0]
    if not max_contamination:
        max_contamination = get_defaults()[1]

    genome_recovery = calc_table(pd_bins_rank, min_completeness, max_contamination)
    counts_dict = collections.OrderedDict()
    for i, cont in zip(range(len(max_contamination)), max_contamination):
        for j, compl in zip(range(len(min_completeness)), min_completeness):
            counts_dict['>{}compl<{}cont'.format(compl, cont)] = genome_recovery[i][j]
    return counts_dict


def print_table(genome_recovery, label, min_completeness, max_contamination, stream=sys.stdout):
    if not label:
        label = ""
    if not min_completeness:
        min_completeness = get_defaults()[0]
    if not max_contamination:
        max_contamination = get_defaults()[1]

    completeness = "".join(">{}% complete\t".format(str(int(x * 100))) for x in min_completeness).rstrip("\t")
    contamination = ["<{}% contamination".format(str(int(x * 100))) for x in max_contamination]

    stream.write("%s\t%s\n" % (label, completeness))

    for y_label, num_genomes in zip(contamination, genome_recovery):
        contamination_values = "".join("{}\t".format(str(values)) for values in num_genomes).rstrip("\t")
        stream.write("%s\t%s\n" % (y_label, contamination_values))


def main():
    parser = argparse.ArgumentParser(description="Calculate number of genome bins recovered with more than the specified thresholds of completeness and contamination. Default: >50%, >70%, >90% completeness vs. <10%, <5% contamination",
                                     parents=[argparse_parents.PARSER_MULTI])
    parser.add_argument('-x', '--min_completeness', help=argparse_parents.HELP_THRESHOLDS_COMPLETENESS, required=False)
    parser.add_argument('-y', '--max_contamination', help=argparse_parents.HELP_THRESHOLDS_CONTAMINATION, required=False)
    args = parser.parse_args()
    if not args.file and sys.stdin.isatty():
        parser.print_help()
        parser.exit(1)
    metrics = load_data.load_tsv_table(sys.stdin if not sys.stdin.isatty() else args.file)

    min_completeness = None
    max_contamination = None
    if args.min_completeness:
        min_completeness = [int(x.strip())/100.0 for x in args.min_completeness.split(',')]
    if args.max_contamination:
        max_contamination = [int(x.strip())/100.0 for x in args.max_contamination.split(',')]

    if args.filter:
        metrics = filter_tail.filter_tail(metrics, args.filter)
    results = calc_table(metrics, min_completeness, max_contamination)
    print_table(results, args.label, min_completeness, max_contamination)


if __name__ == "__main__":
    main()
