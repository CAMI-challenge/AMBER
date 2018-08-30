#!/usr/bin/env python

import collections

TOOL = "tool"
RANK = "rank"
BINNING_TYPE = "binning type"

AVG_PRECISION = "average purity"
AVG_PRECISION_STD = "standard deviation of average purity"
AVG_PRECISION_SEM = "standard error of average purity"

AVG_RECALL = "average completeness"
AVG_RECALL_STD = "standard deviation of average completeness"
AVG_RECALL_SEM = "standard error of average completeness"

AVG_PRECISION_PER_BP = "average purity per bp"

AVG_RECALL_PER_BP = "average completeness per bp"

RI_BY_BP = "Rand index by bp counts"
RI_BY_SEQ = "Rand index by sequence counts"
ARI_BY_BP = "adjusted Rand index by bp counts"
ARI_BY_SEQ = "adjusted Rand index by sequence counts"

PERCENTAGE_ASSIGNED_BPS = "percentage of assigned bps"
ACCURACY = "accuracy"
MISCLASSIFICATION = "misclassification rate"

# PRECISION = "avg_purity_per_bp"
# RECALL = "avg_completeness_per_bp"
# AVG_PRECISION = "avg_purity"
# STD_DEV_PRECISION = "std_dev_purity"
# SEM_PRECISION = "sem_purity"
# AVG_RECALL = "avg_completeness"
# STD_DEV_RECALL = "std_dev_completeness"
# SEM_RECALL = "sem_completeness"
# RI_BY_BP = "rand_index_by_bp"
# RI_BY_SEQ = "rand_index_by_seq"
# ARI_BY_BP = "a_rand_index_by_bp"
# ARI_BY_SEQ = "a_rand_index_by_seq"
# PERCENTAGE_ASSIGNED_BPS = "percent_assigned_bps"
# ACCURACY = "accuracy"

# abbreviations = collections.OrderedDict([(AVG_PRECISION, 'purity averaged over genome bins'),
#                  (STD_DEV_PRECISION, 'standard deviation of purity averaged over genome bins'),
#                  (SEM_PRECISION, 'standard error of the mean of purity averaged over genome bins'),
#                  (AVG_RECALL, 'completeness averaged over genome bins'),
#                  (STD_DEV_RECALL, 'standard deviation of completeness averaged over genome bins'),
#                  (SEM_RECALL, 'standard error of the mean of completeness averaged over genome bins'),
#                  (PRECISION, 'average purity per base pair'),
#                  (RECALL, 'average completeness per base pair'),
#                  ('rand_index_by_bp', 'Rand index weighed by base pairs'),
#                  ('rand_index_by_seq', 'Rand index weighed by sequence counts'),
#                  ('a_rand_index_by_bp', 'adjusted Rand index weighed by base pairs'),
#                  ('a_rand_index_by_seq', 'adjusted Rand index weighed by sequence counts'),
#                  ('percent_assigned_bps', 'percentage of base pairs that were assigned to bins'),
#                  ('accuracy', 'accuracy'),
#                  ('>0.5compl<0.1cont', 'number of bins with more than 50% completeness and less than 10% contamination'),
#                  ('>0.7compl<0.1cont', 'number of bins with more than 70% completeness and less than 10% contamination'),
#                  ('>0.9compl<0.1cont', 'number of bins with more than 90% completeness and less than 10% contamination'),
#                  ('>0.5compl<0.05cont', 'number of bins with more than 50% completeness and less than 5% contamination'),
#                  ('>0.7compl<0.05cont', 'number of bins with more than 70% completeness and less than 5% contamination'),
#                  ('>0.9compl<0.05cont', 'number of bins with more than 90% completeness and less than 5% contamination')])
