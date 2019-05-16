#!/usr/bin/env python

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
        bin_id = bin['mapping_id']
        if bin_id not in genome_to_unique_common or (keyword and genome_to_unique_common[bin_id] != keyword):
            filtered_bin_metrics.append(bin)
    if len(filtered_bin_metrics) == 0:
        sys.exit('All bins have been filtered out due to option --remove_genomes or --min_length.')
    return filtered_bin_metrics
