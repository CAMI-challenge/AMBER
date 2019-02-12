#!/usr/bin/env python

TOOL = "tool"
RANK = "rank"
BINNING_TYPE = "binning type"

AVG_PRECISION_BP = "Average purity (bp)"
AVG_PRECISION_BP_SEM = "Std error of average purity (bp)"

AVG_PRECISION_SEQ = "Average purity (frag)"
AVG_PRECISION_SEQ_SEM = "Std error of average purity (frag)"

AVG_RECALL_BP = "Average completeness (bp)"
AVG_RECALL_BP_SEM = "Std error of average completeness (bp)"

AVG_RECALL_SEQ = "Average completeness (frag)"
AVG_RECALL_SEQ_SEM = "Std error of average completeness (frag)"

AVG_PRECISION_PER_BP = "Average purity per bp"
AVG_PRECISION_PER_SEQ = "Average purity per frag"

AVG_RECALL_PER_BP = "Average completeness per bp"
AVG_RECALL_PER_SEQ = "Average completeness per frag"

RI_BY_BP = "Rand index by bp counts"
RI_BY_SEQ = "Rand index by sequence counts"
ARI_BY_BP = "Adjusted Rand index by bp counts"
ARI_BY_SEQ = "Adjusted Rand index by sequence counts"

PERCENTAGE_ASSIGNED_BPS = "% assigned bps of known origin"
PERCENTAGE_ASSIGNED_SEQS = "% assigned sequences of known origin"

PERCENTAGE_ASSIGNED_BPS_UNKNOWN = "% assigned bps of unknown origin"
PERCENTAGE_ASSIGNED_SEQS_UNKNOWN = "% assigned sequences of unknown origin"

ACCURACY_PER_BP = "Accuracy per bp"
ACCURACY_PER_SEQ = "Accuracy per frag"
MISCLASSIFICATION_PER_BP = "Misclassification rate per bp"
MISCLASSIFICATION_PER_SEQ = "Misclassification rate per frag"
