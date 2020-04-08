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

import sys
import os

try:
    import argparse_parents
    import load_data
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import argparse_parents
        import load_data
    finally:
        sys.path.remove(os.path.dirname(__file__))


def filter_data(bin_metrics, genome_to_unique_common, keyword):
    filtered_bin_metrics = []
    for bin in bin_metrics:
        bin_id = bin['most_abundant_genome']
        if bin_id not in genome_to_unique_common or (keyword and genome_to_unique_common[bin_id] != keyword):
            filtered_bin_metrics.append(bin)
    if len(filtered_bin_metrics) == 0:
        sys.exit('All bins have been filtered out due to option --remove_genomes or --min_length.')
    return filtered_bin_metrics
