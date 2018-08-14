#!/usr/bin/env python

import argparse
import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src import precision_recall_per_bin
from src import precision_recall_average
from src.utils import argparse_parents
from src.utils import exclude_genomes
from src.utils import load_data


def evaluate_all(gold_standard, queries, labels, filter_tail_percentage, genomes_file, keyword, map_by_recall):
    precision_recall_average.print_precision_recall_table_header()
    labels_iterator = iter(labels)
    for query in queries:

        bin_metrics = precision_recall_per_bin.compute_metrics(query, gold_standard, map_by_recall)

        if genomes_file:
            bin_metrics = exclude_genomes.filter_data(bin_metrics, genomes_file, keyword)

        avg_precision, avg_recall, std_deviation_precision, std_deviation_recall, std_error_precision, std_error_recall = \
            precision_recall_average.compute_precision_and_recall(bin_metrics, filter_tail_percentage)

        precision_recall_average.print_precision_recall(next(labels_iterator) if len(labels) > 0 else query.path.split('/')[-1],
                                                        avg_precision,
                                                        avg_recall,
                                                        std_deviation_precision,
                                                        std_deviation_recall,
                                                        std_error_precision,
                                                        std_error_recall)
