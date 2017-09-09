#!/usr/bin/env python

import collections

TOOL = "tool"
PRECISION = "precision"
RECALL = "recall"
AVG_PRECISION = "avg_precision"
STD_DEV_PRECISION = "std_dev_precision"
SEM_PRECISION = "sem_precision"
AVG_RECALL = "avg_recall"
STD_DEV_RECALL = "std_dev_recall"
SEM_RECALL = "sem_recall"
RI_BY_BP = "rand_index_by_bp"
RI_BY_SEQ = "rand_index_by_seq"
ARI_BY_BP = "a_rand_index_by_bp"
ARI_BY_SEQ = "a_rand_index_by_seq"
PERCENTAGE_ASSIGNED_BPS = "percent_assigned_bps"
ACCURACY = "accuracy"

abbreviations = collections.OrderedDict([('avg_precision', 'precision averaged over genome bins'),
                 ('std_dev_precision', 'standard deviation of precision averaged over genome bins'),
                 ('sem_precision', 'standard error of the mean of precision averaged over genome bins'),
                 ('avg_recall', 'recall averaged over genome bins'),
                 ('std_dev_recall', 'standard deviation of recall averaged over genome bins'),
                 ('sem_recall', 'standard error of the mean of recall averaged over genome bins'),
                 ('precision', 'precision weighed by base pairs'),
                 ('recall', 'recall weighed by base pairs'),
                 ('rand_index_by_bp', 'Rand index weighed by base pairs'),
                 ('rand_index_by_seq', 'Rand index weighed by sequence counts'),
                 ('a_rand_index_by_bp', 'adjusted Rand index weighed by base pairs'),
                 ('a_rand_index_by_seq', 'adjusted Rand index weighed by sequence counts'),
                 ('percent_assigned_bps', 'percentage of base pairs that were assigned to bins'),
                 ('accuracy', 'accuracy'),
                 ('>0.5compl<0.1cont', 'number of bins with more than 50% completeness and less than 10% contamination'),
                 ('>0.7compl<0.1cont', 'number of bins with more than 70% completeness and less than 10% contamination'),
                 ('>0.9compl<0.1cont', 'number of bins with more than 90% completeness and less than 10% contamination'),
                 ('>0.5compl<0.05cont', 'number of bins with more than 50% completeness and less than 5% contamination'),
                 ('>0.7compl<0.05cont', 'number of bins with more than 70% completeness and less than 5% contamination'),
                 ('>0.9compl<0.05cont', 'number of bins with more than 90% completeness and less than 5% contamination')])
